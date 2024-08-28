import numpy as np
import pysages
from pysages.colvars import Distance, DihedralAngle
from pysages.methods import ABF, SpectralABF, Unbiased, HistogramLogger
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.langevin import Langevin
from ase.io.trajectory import Trajectory
from ase import units
from ase.io import read
from ase.calculators.amber import Amber
import matplotlib.pyplot as plt
import argparse
import dill as pickle

# Temperature and timestep
T = 300
dt = 1.0 * units.fs
PI = np.pi

# Simulation
def generate_simulation(T=T, dt=dt):
    atoms = read('ala-ase.xyz')
    calc = Amber(amber_exe='sander -O ',
                infile='mm.in',
                outfile='mm.out',
                topologyfile='strip.parm7',
                incoordfile='strip.rst7'
                )
    atoms.calc = calc

    MaxwellBoltzmannDistribution(atoms, temperature_K=T)
    dyn = Langevin(atoms, dt, temperature_K=T, friction=0.02)
    def _printenergy(a=atoms):
        """Function to print the potential, kinetic and total energy"""
        epot = a.get_potential_energy() / len(a)
        ekin = a.get_kinetic_energy() / len(a)
        with open("energy.txt", "a") as f:
          print('Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  '
            'Etot = %.3feV' % (epot, ekin, ekin / (1.5 * units.kB), epot + ekin), file=f)

    dyn.attach(_printenergy, interval=3000)
    traj = Trajectory('moldyn3.traj', 'w', atoms)
    dyn.attach(traj.write, interval=5000)

    return dyn

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--unbiased", action="store_true")
    return parser.parse_args()

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

def plot_histogram(result, filename="histogram.png"):
    surface = np.asarray(result["histogram"])
    plot_helper(surface, filename, "Histogram")

def save_energy_forces(result):
    energy = np.asarray(result["free_energy"])
    forces = np.asarray(result["mean_force"])
    grid = np.asarray(result["mesh"])
    hist = np.asarray(result["histogram"])
    
    np.savetxt("FES.csv", np.hstack([grid, energy.reshape(-1,1)]))
    np.savetxt("Histogram.csv", np.hstack([grid, hist.reshape(-1,1)]))


def main(args):
    cvs = [
        DihedralAngle([4, 6, 8, 14]),DihedralAngle([6, 8, 14, 16])
    ]
    shape = (32, 32)
    grid = pysages.Grid(
        lower=(-PI,-PI), 
        upper=(PI, PI), 
        shape=shape, 
        periodic=True
        )
    
    if args.unbiased:
        callback = HistogramLogger(5000)
        method = Unbiased(cvs)
        state = pysages.run(method, generate_simulation, 5000000, callback=callback)
    else:
        method = SpectralABF(cvs, grid)
        state = pysages.run(method, generate_simulation, 5000000)
    
    # Save the result
    with open("restart.pickle", "wb") as f:
        pickle.dump(state, f)

    if not args.unbiased:
        result = pysages.analyze(state)
        plot_data(result, "free_energy", "energy.png", shape=shape)
        plot_histogram(result)
        save_energy_forces(result)
    else:
        bins = 50
        hist, xedges,yedges = np.histogram2d(
                                callback.data[:,0], 
                                callback.data[:,1],
                                bins=bins,range=[[-PI, PI],[-PI, PI]]
                              )
        xedges = (xedges[1:] + xedges[:-1]) / 2
        yedges = (yedges[1:] + yedges[:-1]) / 2
        mesh = np.dstack(np.meshgrid(xedges, yedges)).reshape(-1, 2)
        np.savetxt("hist-fromlogger.csv", np.hstack([mesh, hist.reshape(-1, 1)]))


if __name__ == "__main__":
    args = get_args()
    main(args)
