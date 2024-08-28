import numpy as np


data_file='../alanine_dipeptide.npy'
#ranges = [train, val, test]
all_data = np.load(data_file, allow_pickle=True).item()
all_data['force'] = all_data['force'] / 96.0
data = all_data
data['pbc'] = np.array([True]*3)
data['energy'] = np.zeros(50001)[:, None]
data['force'] = data['force']
data['lattices'] = data['lengths'][:, None] * np.eye(3)
np.savez('data_corr', **data)
