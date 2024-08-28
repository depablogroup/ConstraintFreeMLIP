#!/usr/bin/env python3

from pysages.utils import try_import
openmm = try_import("openmm", "simtk.openmm")
unit = try_import("openmm.unit", "simtk.unit")
app = try_import("openmm.app", "simtk.openmm.app")
from openmm import *
from openmm.unit import *
from parmed import load_file, unit as u
from parmed.openmm import StateDataReporter


print("Loading input files...")
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
    300 * kelvin,  # Temperature of heat bath
    1.0 / picoseconds,  # Friction coefficient
    1.0 * femtoseconds,  # Time step
)
platform = Platform.getPlatformByName("CUDA")
    
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


sim.reporters.append(app.DCDReporter("output2.dcd", 1))
sim.reporters.append(StateDataReporter("data.txt", 1, step=True, separator=" "))
sim.reporters.append(ForceReporter('forces.txt', 1))


for i in range(1000):
    pdb = app.PDBFile(str(i)+'.pdb')
    sim.context.setPositions(pdb.getPositions())
    sim.step(1)
