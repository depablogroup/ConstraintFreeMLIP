
import numpy as np
from ase.io import read
from ase import units
from nequip.ase.nequip_calculator import NequIPCalculator


# Temperature and timestep
T = 300
dt = 1.0 * units.fs

# Initialization
binsx, binsy = 20, 20
hist = np.zeros((binsx, binsy))
numbers = np.zeros((binsx, binsy))
cv1min, cv1max = -np.pi, np.pi
cv2min, cv2max = -np.pi, np.pi
dcv1 = (cv1max - cv1min) / binsx
dcv2 = (cv2max - cv2min) / binsy
x = np.linspace(cv2min, cv2max, binsx)
y = np.linspace(cv1min, cv1max, binsx)
X,Y = np.meshgrid(x, y)
XY = np.array([X.flatten(), Y.flatten()]).T
grid = XY


dimer = read('ala-ase.xyz')
dimer.set_pbc([0, 0, 0])
calc= NequIPCalculator.from_deployed_model('model_corrected.pth',device='cuda')
dimer.set_calculator(calc)

for i in range(80000):
    atoms = read('test80.extxyz', index=i)
    positions = atoms.get_positions()
    dimer.set_positions(positions)
    e = dimer.get_potential_energy()
    phi = (dimer.get_dihedral(4, 6, 8, 14) / 180.0) * np.pi - 2 * np.pi
    if phi < -np.pi:
        phi = phi + 2 * np.pi
    psi = (dimer.get_dihedral(6, 8, 14, 16) / 180.0) * np.pi - 2 * np.pi
    
    if psi < -np.pi:
        psi = psi + 2 * np.pi
    index_x = int((phi - cv1min) / dcv1)
    index_y = int((psi - cv2min) / dcv2)

    if 0 <= index_x < binsx and 0 <= index_y < binsy:
        hist[index_x][index_y] += e
        numbers[index_x][index_y] += 1
        
    print(i)

for i in range(binsx):
    for j in range(binsy):
        if numbers[i,j] > 0:
            hist[i][j] = hist[i][j] / numbers[i][j]

np.savetxt('enthalpycorr.txt', np.column_stack([grid, hist.flatten()]))
np.savetxt('histcorr.txt', np.column_stack([grid, numbers.flatten()]))
