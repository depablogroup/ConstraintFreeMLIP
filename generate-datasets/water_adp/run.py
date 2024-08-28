# Importing necessary modules
import argparse
from functools import partial
from pysages.methods import Unbiased, SpectralABF
from pysages.colvars import DihedralAngle
from pysages.backends import SamplingContext
from pysages.utils import try_import
import numpy as np
import pysages
import dill as pickle
import matplotlib.pyplot as plt


openmm = try_import("openmm", "simtk.openmm")
unit = try_import("openmm.unit", "simtk.unit")
app = try_import("openmm.app", "simtk.openmm.app")
from openmm import *
from openmm.unit import *

from parmed import load_file, unit as u
from parmed.openmm import StateDataReporter

PI = np.pi
T = 600

def generate_simulation(T=600.0 * kelvin, dt=1.0 * femtoseconds):
    print("Loading AMBER files...")
    pdb = app.PDBFile('mixture.pdb')
    ff = app.ForceField("amber99sb.xml", "tip3p.xml")
    topology = pdb.topology
    system = ff.createSystem(topology,
        nonbondedMethod=app.PME,
        rigidWater=False,
        switchDistance=0.6 * nanometer,
        nonbondedCutoff=0.7 * nanometer,
        constraints=None,
    )
    # Create the integrator to do Langevin dynamics
    integrator = openmm.LangevinIntegrator(
        T,  # Temperature of heat bath
        1.0 / picoseconds,  # Friction coefficient
        dt,  # Time step
    )
    platform = Platform.getPlatformByName("CUDA")
    # Create the Simulation object
    class ForceReporter(object):
        def __init__(self, file, reportInterval):
            self._out = open(file, 'w')
            self._reportInterval = reportInterval

        def __del__(self):
            self._out.close()

        def describeNextReport(self, simulation):
            steps = self._reportInterval - simulation.currentStep%self._reportInterval
            return (steps, False, False, True, False, None)

        def report(self, simulation, state):
            forces = state.getForces().value_in_unit(kilojoules/mole/nanometer)
            for f in forces:
                self._out.write('%g %g %g\n' % (f[0], f[1], f[2]))
    sim = app.Simulation(topology, system, integrator, platform)

    # Set the particle positions
    sim.context.setPositions(pdb.getPositions(frame=-1))
    sim.reporters.append(app.PDBReporter("output.pdb", 100000))
    sim.reporters.append(app.DCDReporter("output.dcd", 2000))
    sim.reporters.append(StateDataReporter("data.txt", 2000, step=True, separator=" "))
    sim.reporters.append(ForceReporter('forces.txt', 2000))

    return sim


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


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--unbiased", action="store_true")
    parser.add_argument("--restart", action="store_true")
    return parser.parse_args()


def main(args):
    cvs = [DihedralAngle([4, 6, 8, 14]),DihedralAngle([6, 8, 14, 16])]
    timesteps = 80000000
    if args.unbiased:
        method = Unbiased(cvs)
    else:
        shape = (50, 50)
        grid = pysages.Grid(
            lower=(-PI, -PI), 
            upper=(PI, PI), 
            shape=shape, 
            periodic=True
        )
        method = SpectralABF(cvs, grid=grid)
    sampling_context = SamplingContext(method, generate_simulation)

    # Potentially restart from file
    if args.restart:
        with open("restart.pickle", "rb") as f:
            state = pickle.load(f)
        state = pysages.run(state, generate_simulation, timesteps)
    else:
        state = pysages.run(sampling_context, timesteps)

    # Save the result
    with open("restart.pickle", "wb") as f:
        pickle.dump(state, f)

    if not args.unbiased:
        result = pysages.analyze(state)
        plot_data(result, "free_energy", "energy.png", shape=shape)
        plot_histogram(result)
        save_energy_forces(result)

if __name__ == "__main__":
    args = get_args()
    main(args)
