from ase.io.trajectory import Trajectory
from ase.io import write


traj = Trajectory('moldyn3.traj')
for i in range(len(traj)):
    atoms=traj[i]
    write('test.extxyz', atoms, append=True)
exit()
