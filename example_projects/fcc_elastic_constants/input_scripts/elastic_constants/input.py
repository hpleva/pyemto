import numpy as np
import ase
from pyemto.utilities import distort, rotation_matrix
from pyemto.examples.emto_input_generator import *
from pymatgen import Lattice
from ase.visualize import view
from ase.build import cut, make_supercell
import sys

find_primitive = False
make_supercell = None
coords_are_cartesian = True

nkx = 13
nky = 13
nkz = 13

# Primitive fcc
prims0 = np.array([
    [0.5,0.5,0],
    [0.5,0,0.5],
    [0,0.5,0.5]])

basis0 = np.array([
    [0.0,0.0,0.0]
])

species = [['Co','Co','Cr','Cr','Fe','Fe','Mn','Mn','Ni','Ni']]
splts = [[1,-1]*5]

ncpu = 1

folder = os.getcwd()
emtopath = folder
latpath = emtopath

slurm_options=['#SBATCH -A snic2018-3-205',
               '#SBATCH -n {0}'.format(ncpu),
               #'#SBATCH --mem-per-cpu=4000',
               #'module load emto/5.8.1/build16-QNA'
]


deltas = np.linspace(0,0.05,6)
# We need to use a non-zero value for the first delta to break the symmetry of the structure.
deltas[0] = 0.001

# Only two distortions for cubic (third one is bulk modulus EOS fit)
distortions = ['Cprime', 'C44']

for i, distortion in enumerate(distortions):
    print('#'*100)
    print('distortion = ',distortion)
    print('#'*100)
    for delta in deltas:
        print('#'*100)
        print('delta = ',delta)
        print('#'*100)

        # These distortion matrices are from the EMTO book.
        
        if distortion == 'Cprime':
            dist_matrix = np.array([
                    [1+delta,0,0],
                    [0,1-delta,0],
                    [0,0,1/(1-delta**2)]
                    ])

        elif distortion == 'C44':
            dist_matrix = np.array([
                    [1, delta, 0],
                    [delta ,1, 0],
                    [0, 0, 1/(1-delta**2)]
                    ])

        # Calculate new lattice vectors and atomic positions
        prims = distort(dist_matrix, prims0)
        basis = distort(dist_matrix, basis0)

        # Each different distortion might need different set of nkx, nky, nkz
        if distortion == 'Cprime':
            nkx = 21; nky = 21; nkz = 21
        elif distortion == 'C44':
            nkx = 20; nky = 20; nkz = 25

        input_creator = EMTO(folder=emtopath, EMTOdir='/home/x_henle/emto/KED_ELF/emto-dev/build')
        input_creator.prepare_input_files(latpath=latpath,
                                          jobname='CoCrFeMnNi_d{0}_{1:4.2f}'.format(i+1,delta),
                                          species=species,
                                          splts=splts,
                                          #concs=concs,
                                          prims=prims,
                                          basis=basis,
                                          find_primitive=find_primitive,
                                          coords_are_cartesian=coords_are_cartesian,
                                          latname='d{0}_{1:4.2f}'.format(i+1, delta),
                                          #nz1=32,
                                          ncpa=15,
                                          sofc='Y',
                                          nkx=nkx,
                                          nky=nky,
                                          nkz=nkz,
                                          ncpu=ncpu,
                                          parallel=False,
                                          alpcpa=0.6,
                                          runtime='24:00:00',
                                          KGRN_file_type='scf',
                                          KFCD_file_type='fcd',
                                          amix=0.01,
                                          #efgs=-1.0,
                                          #depth=2.0,
                                          tole=1e-6,
                                          tolef=1e-6,
                                          iex=4,
                                          niter=200,
                                          kgrn_nfi=91,
                                          #strt='B',
                                          make_supercell=make_supercell,
                                          slurm_options=slurm_options)
        
        sws_range = np.array([2.602])
        
        input_creator.write_bmdl_kstr_shape_input()
        input_creator.write_kgrn_kfcd_swsrange(sws=sws_range)
