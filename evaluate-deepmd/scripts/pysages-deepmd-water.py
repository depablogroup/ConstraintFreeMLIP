from ase.io import read
from deepmd.calculator import DP
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.langevin import Langevin
from ase import units
from ase import Atoms
import numpy as np
import pysages
import dill as pickle
from pysages.colvars import DihedralAngle
from pysages.methods import SpectralABF
import matplotlib.pyplot as plt
import os
import argparse
from functools import partial


T = 300.0
dt = 1.0 * units.fs

pi = np.pi

def generate_simulation(model: str='frozen_model.pb', T=T, restart=False):
    """Generate a simulation object with ase backend for PySAGES."""

    # Read in the simulation system
    if restart:
        dimer= read('md.traj', '-1')
    else:
        dimer: Atoms = read('mixture.xyz')
        dimer.set_cell([15.000, 15.000, 15.000])
        dimer.set_pbc([1, 1, 1])

    # Set up MLFF as ase calculator
    calc= DP(model)
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
    return parser.parse_args()

def get_model(savedir):
    """Get the path of the model in the directory."""
    for filename in os.listdir(savedir):
        if filename.endswith(".pb"):
            return os.path.join(savedir, filename)
    raise ValueError(f"Model not found in {savedir}!")


# Functions for ploting and storing data
def plot_energy(result):
    energy = result["free_energy"]
    mesh = np.asarray(result["mesh"])
    A = energy
    A = (A.max() - A).reshape((32,32))
    surface = A.reshape(32,32)
    np.savetxt("FES-5ns.csv", np.hstack([mesh, A.reshape(-1, 1)]))

    fig, ax = plt.subplots()
    im = ax.imshow(
        surface, 
        interpolation="bicubic", 
        origin="lower", 
        extent=[-pi, pi, -pi, pi], 
        aspect=1
    )
    ax.contour(surface, 
               levels=15, 
               linewidths=0.75, 
               colors="k", 
               extent=[-pi, pi, -pi, pi]
               )
    plt.colorbar(im)
    fig.savefig("energy.png")


def plot_histogram(result):
    surface = np.asarray(result["histogram"])
    surface=surface.reshape(32,32)

    fig, ax = plt.subplots()
    im = ax.imshow(
        surface, 
        interpolation="bicubic", 
        origin="lower", 
        extent=[-pi, pi, -pi, pi], 
        aspect=1
    )
    ax.contour(surface,
               levels=15,
               linewidths=0.75,
               colors="k",
               extent=[-pi, pi, -pi, pi]
               )
    plt.colorbar(im)
    fig.savefig("histogram.png")


def main(args):
    """Main function to perform simulation with Spectral-ABF."""
    # Setting up the collective variables of interest
    cvs = [DihedralAngle([4, 6, 8, 14]), DihedralAngle([6, 8, 14, 16])]
    grid = pysages.Grid(
        lower=(-np.pi, -np.pi),
        upper=(np.pi, np.pi),
        shape=(32, 32),
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

    # Extract free energy from result
    result = pysages.analyze(raw_result)

    plot_energy(result)
    plot_histogram(result)


if __name__ == "__main__":
    args = get_args()
    main(args)
