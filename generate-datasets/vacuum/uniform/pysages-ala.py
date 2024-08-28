
import numpy as np
import pysages
from pysages.colvars import Distance, DihedralAngle
from pysages.methods import ABF, SpectralABF
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase.md.langevin import Langevin
from ase.io.trajectory import Trajectory
from ase import units
from ase.io import read
from ase.calculators.amber import Amber
from ase.build import molecule
from ase import Atoms
import matplotlib.pyplot as plt


# Temperature and timestep
T = 500
dt = 1.0 * units.fs
PI = np.pi

# Simulation
def generate_simulation(T=T, dt=dt):

    atoms = read('ala-ase.xyz')
    calc = Amber(amber_exe='sander -O ',
                infile='mm.in',
                outfile='mm.out',
                topologyfile='strip.parm7',
                incoordfile='strip.rst7')

    atoms.calc = calc
    MaxwellBoltzmannDistribution(atoms, temperature_K=T)
    dyn = Langevin(atoms, dt,temperature_K=T, friction=0.02)
    def _printenergy(a=atoms):
        """Function to print the potential, kinetic and total energy"""
        epot = a.get_potential_energy() / len(a)
        ekin = a.get_kinetic_energy() / len(a)
        with open("energy.txt", "a") as f:
          print('Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  '
            'Etot = %.3feV' % (epot, ekin, ekin / (1.5 * units.kB), epot + ekin), file=f)

    dyn.attach(_printenergy, interval=2000)

    traj = Trajectory('moldyn3.traj', 'w', atoms)
    dyn.attach(traj.write, interval=2000)

    return dyn


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
    surface = np.asarray(result["histogram"]) / 300
    plot_helper(surface, filename, "Histogram")

def save_energy_forces(result):
    energy = np.asarray(result["free_energy"])
    forces = np.asarray(result["mean_force"])
    grid = np.asarray(result["mesh"])
    hist = np.asarray(result["histogram"])
    
    np.savetxt("FES.csv", np.hstack([grid, energy.reshape(-1,1)]))
    np.savetxt("Histogram.csv", np.hstack([grid, hist.reshape(-1,1)]))


def main():
    cvs = [
        DihedralAngle([4,6,8,14]),DihedralAngle([6,8,14,16])
    ]
    shape = (50, 50)
    grid = pysages.Grid(
        lower=(-PI, -PI), 
        upper=(PI, PI), 
        shape=(50, 50), 
        periodic=True
        )
    method = SpectralABF(cvs, grid)

    state = pysages.run(method, generate_simulation, 2000000)
    result = pysages.analyze(state)
    plot_data(result, "free_energy", "energy.png", shape=shape)
    plot_histogram(result)
    save_energy_forces(result)

if __name__ == "__main__":
    main()
