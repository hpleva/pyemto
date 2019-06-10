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
coa_array = np.linspace(1.3,1.8,7)

ncpu = 16

folder = os.getcwd()
emtopath = folder
latpath = emtopath

slurm_options=['#SBATCH -n {0}'.format(ncpu),
               '#SBATCH --ntasks-per-node=4',
               '##SBATCH --qos=short',
               'ulimit -s unlimited',
               'module load compiler/intel/16 intel-mkl/16 openmpi/2.0',
               'export OMP_NUM_THREADS=${SLURM_NTASKS_PER_NODE}',
               'export OMP_STACKSIZE=400m']


# Calculate equilibrium volume

for i in range(len(coa_array)):

    coa = coa_array[i]

    # Primitive hcp
    prims0 = np.array([
    [1.0000000000,0.0000000000,0.0000000000],
    [-0.5000000000,np.sqrt(3.)/2,0.0000000000],
    [0.0000000000,0.0000000000,coa]
    ])

    basis0 = np.array([
    [0.0,0.0,0.0],
    [0.0,0.57735027,coa/2.0]
    ])

    species = [['Ti','Ti'],['Ti','Ti']]
    splts = [[1,-1],[1,-1]]

    concs=[[0.5,0.5],[0.5,0.5]]

    input_creator = EMTO(folder=emtopath, EMTOdir='/home/minds/jianwang/vasptool/EMTO/5.8.1/build')
    input_creator.prepare_input_files(latpath=latpath,
                                      jobname='Ti_hcp_ca{}'.format(i+1),
                                      species=species,
                                      splts=splts,
                                      concs=concs,
                                      prims=prims0,
                                      basis=basis0,
                                      find_primitive=find_primitive,
                                      coords_are_cartesian=coords_are_cartesian,
                                      latname='hcp_ca{}'.format(i+1),
                                      #nz1=32,
                                      ncpa=20,
                                      sofc='Y',
                                      nkx=nkx,
                                      nky=nky,
                                      nkz=nkz,
                                      ncpu=ncpu,
                                      parallel=False,
                                      alpcpa=0.6,
                                      runtime='06:00:00',
                                      KGRN_file_type='scf',
                                      KFCD_file_type='fcd',
                                      amix=0.01,
                                      #efgs=-1.0,
                                      #depth=2.0,
                                      tole=1e-6,
                                      iex=4,
                                      niter=200,
                                      kgrn_nfi=91,
                                      #strt='B',
                                      make_supercell=make_supercell,
                                      slurm_options=slurm_options)

    sws_range = np.linspace(2.5, 2.75, 15)

    input_creator.write_bmdl_kstr_shape_input()
    input_creator.write_kgrn_kfcd_swsrange(sws=sws_range)
