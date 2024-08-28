"""Use PySAGES, ASE and Allegro to learn the FES of molecule with Spectral-ABF Sampling."""
import argparse
from functools import partial
import os
import numpy as np
from ase import Atoms
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase import units
from ase.io import read
from ase.constraints import FixBondLengths
import matplotlib.pyplot as plt
import pysages
import dill as pickle
from pysages.colvars import DihedralAngle
from pysages.methods import SpectralABF
from nequip.ase.nequip_calculator import NequIPCalculator


# Temperature and timestep
T = 300
dt = 1.0 * units.fs
PI = np.pi

def generate_simulation(model: str='adpmodelexpu1k.pth', T=T, restart=False, fixbond=False):
    """Generate a simulation object with ase backend for PySAGES."""

    # Read in the simulation system
    if restart:
        dimer= read('md.traj', '-1')
    else:
        dimer: Atoms = read('ala-ase.xyz')
        dimer.set_cell([32.168, 32.168, 32.168])
        dimer.set_pbc([1, 1, 1])

    if fixbond:
        constr = FixBondLengths([[1, 0], [1, 2], [1, 3],
                                 [8, 9], [10, 11], [10, 12],
                                 [10, 13], [16, 17],[18, 19],
                                 [18, 20],[18, 21], [6,7]])
        dimer.set_constraint(constr)


    # Set up MLFF as ase calculator
    calc = NequIPCalculator.from_deployed_model(model, device='cuda')
    dimer.set_calculator(calc)

    # Set up the NVT simulation object
    MaxwellBoltzmannDistribution(dimer, temperature_K=T)
    dyn = Langevin(
        dimer,
        dt,
        temperature_K=T,
        friction=0.01,
        trajectory='md.traj',
        logfile='md.log',
        loginterval=5000
    )
    return dyn


def get_args():
    """Get args from input."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--savedir", type=str, required=True)
    parser.add_argument("--expname", type=str, default="5ns")
    parser.add_argument("--restart", action="store_true")
    parser.add_argument("--fixbond", action="store_true")

    return parser.parse_args()


def get_model(savedir):
    """Get the path of the model in the directory."""
    for filename in os.listdir(savedir):
        if filename.endswith(".pth"):
            return os.path.join(savedir, filename)
    raise ValueError(f"Model not found in {savedir}!")

def plot_helper(surface: np.ndarray, filename: str, title: str):
    fig, ax = plt.subplots()
    im = ax.imshow(
        surface, 
        interpolation="bicubic", 
        origin="lower", 
        extent=[-PI, PI, -PI, PI], 
        aspect=1
    )
    ax.contour(surface, 
               levels=15, 
               linewidths=0.75, 
               colors="k", 
               extent=[-PI, PI, -PI, PI]
               )
    plt.colorbar(im)
    plt.title(title)
    fig.savefig(filename)

def plot_data(result, data_key: str, filename: str, title: str, shape=(50, 50)):
    data = result[data_key]
    A = data.max() - data
    surface = A.reshape(shape)
    plot_helper(surface, filename, title)

def main(args):
    """Main function to perform simulation with Spectral-ABF."""
    # Setting up the collective variables of interest
    cvs = [DihedralAngle([4, 6, 8, 14]), DihedralAngle([6, 8, 14, 16])]
    shape = (32, 32)
    grid = pysages.Grid(
        lower=(-np.pi, -np.pi),
        upper=(np.pi, np.pi),
        shape=shape,
        periodic=True
    )
    method = SpectralABF(cvs, grid)
    timesteps = 5000000

    model = get_model(args.savedir)
    generate_simulation_func = partial(
        generate_simulation,
        model=model,
        restart=args.restart
    )

    # Potentially restart from file
    if args.restart:
        with open("restart.pickle", "rb") as f:
            raw_result = pickle.load(f)
        raw_result = pysages.run(raw_result, generate_simulation_func, timesteps)
    else:
        raw_result = pysages.run(method, generate_simulation_func, timesteps)

    # Save the result
    with open("restart.pickle", "wb") as f:
        pickle.dump(raw_result, f)

    # Extract free energy from result
    result = pysages.analyze(raw_result)
    plot_data(result, "free_energy", "energy.png", shape=shape)


if __name__ == "__main__":
    args = get_args()
    main(args)
