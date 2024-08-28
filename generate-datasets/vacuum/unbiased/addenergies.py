import numpy as np
from ase.units import kcal,mol,Ang
from ase.io import write, read
from ase.calculators.singlepoint import SinglePointCalculator


out_filename = 'test.extxyz'
x, y, epot1, epot2, w, m = np.genfromtxt("energies.txt", unpack=True)
for i in range(1001):
    dimer = read('noenergy.extxyz', index=i)
    forces = dimer.get_forces()
    calculator = SinglePointCalculator(
        dimer, 
        energy=(epot1[i] + epot2[i])*kcal/mol, 
        forces=forces*kcal/(mol*Ang))
    dimer.calc = calculator
    write(out_filename, dimer, format='extxyz', append=True)
exit()
