from ase.io import read
from kaldo.conductivity import Conductivity
from kaldo.forceconstants import ForceConstants
from kaldo.phonons import Phonons
import numpy as np

### Set up force constant objects via interface to LAMMPS ####

# Replicate the unit cell 'nrep'=3 times
nrep = 3
supercell = np.array([nrep, nrep, nrep])

forceconstants = ForceConstants.from_folder(folder='fc_Si46',supercell=supercell,format='lammps')


# Configure phonon object
# 'k_points': number of k-points
# 'is_classic': specify if the system is classic, True for classical and False for quantum
# 'temperature: temperature (Kelvin) at which simulation is performed
# 'folder': name of folder containing phonon property and thermal conductivity calculations
# 'storage': Format to storage phonon properties ('formatted' for ASCII format data, 'numpy' 
#            for python numpy array and 'memory' for quick calculations, no data stored)

# Define the k-point mesh using 'kpts' parameter
k_points = 3#'k_points'=5 k points in each direction
phonons_config = {'kpts': [k_points, k_points, k_points],
                  'is_classic': False, 
                  'temperature': 300, #'temperature'=300K
                  'folder': 'Si46_ald',
		   'storage': 'formatted'}
# Set up phonon object by passing in configuration details and the forceconstants object computed above
phonons = Phonons(forceconstants=forceconstants, **phonons_config)

### Set up the Conductivity object and thermal conductivity calculations ####

# Compute thermal conductivity (t.c.) by solving Boltzmann Transport
# Equation (BTE) with various of methods.

# 'phonons': phonon object obtained from the above calculations
# 'method': specify methods to solve for BTE  
# ('qhgk' for Quasi Harmonic Green Kubo (QHGK),'rta' for RTA)


print('\n')
qhgk_cond_matrix = Conductivity(phonons=phonons, method='qhgk').conductivity.sum(axis=0)
print('Conductivity from Quasi QHGK (W/m-K): %.3f' % (-1*np.mean(np.diag(qhgk_cond_matrix))))
print(-1*qhgk_cond_matrix)

print('\n')
rta_cond_matrix = Conductivity(phonons=phonons, method='rta').conductivity.sum(axis=0)
print('Conductivity from RTA (W/m-K): %.3f' % (np.mean(np.diag(rta_cond_matrix))))
print(rta_cond_matrix)
