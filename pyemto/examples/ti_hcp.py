#!/usr/bin/env python

import pyemto
import os
import numpy as np

# 1. Calculate the equilibrium volume and bulk modulus of hcp Ti.
# 2. Calculate the elastic constants of hcp Ti.

which_calculation = 'volume'
#which_calculation = 'elastic_constants'

which_mode = 'generate'
#which_mode = 'analyze'

folder = os.getcwd()                # Get current working directory.
emtopath = folder + "/ti_hcp"         # Folder where the calculation take place.
latpath = emtopath                  # Folder where the structure files are located.

ncpu = 16
slurm_options=['#SBATCH -n {0}'.format(ncpu),
               '#SBATCH --ntasks-per-node=4',
               '##SBATCH --qos=short',
               'ulimit -s unlimited',
               'module load compiler/intel/16 intel-mkl/16 openmpi/2.0',
               'export OMP_NUM_THREADS=${SLURM_NTASKS_PER_NODE}',
               'export OMP_STACKSIZE=400m']

ti_hcp=pyemto.System(folder=emtopath, EMTOdir='/home/minds/jianwang/vasptool/EMTO/5.8.1/build')

# Initialize the bulk system using the bulk() function:
ti_hcp.bulk(lat='hcp',
            latpath=latpath,
            atoms=['Ti'],
            sws=3.0,
            amix=0.02,
            efmix=0.9,
            #expan='M',
            sofc='Y',
            xc='P07',           # Use PBEsol
            nky=31,             # k-points
            nkz=19,             # k-points
            runtime='24:00:00', # Allow large enough timelimit for SLURM
            slurm_options=slurm_options,
            KGRN_file_type='scf',
            KFCD_file_type='fcd',
            ncpu=ncpu,
            )

sws = np.linspace(2.9,3.1,7) # A list of 7 different volumes from 2.9 to 3.1

if which_calculation == 'volume' and which_mode == 'generate':

    # Unfortunately the "lattice_constants_batch_generate" does not
    # automatically generate the lattice input files, so we generate them here
    # manually:

    # Create hcp structure files
    cas = np.linspace(1.50,1.7,7) # 7 different equally spaced c/a values from 1.5 to 1.7
    latnames = ['hcp_ca1','hcp_ca2','hcp_ca3','hcp_ca4','hcp_ca5','hcp_ca6','hcp_ca7']

    structure = pyemto.System(folder=emtopath)
    for i in range(len(latnames)):
        structure.lattice.set_values(jobname_lat=latnames[i],
                                     latpath=latpath,
                                     lat='hcp',
                                     #kappaw=[0.0,-20.0],
                                     msgl=0,
                                     ca=cas[i])

        structure.lattice.bmdl.write_input_file(folder=latpath)
        structure.lattice.kstr.write_input_file(folder=latpath)
        structure.lattice.shape.write_input_file(folder=latpath)
        structure.lattice.batch.write_input_file(folder=latpath)

    ti_hcp.lattice_constants_batch_generate(sws=sws, auto_ca=True)

elif which_calculation == 'volume' and which_mode == 'analyze':
    ti_hcp.lattice_constants_analyze(sws=sws)

# Results printed by the lattice_constants_analyze function.
# They are inputed to the elastic_constants functions.
# hcp Ti should get results close to what is below.
# Check Vitos book pages 108-109 for explanation what the R and cs parameters
# are.

sws0 =      3.002260
ca0  =      1.610122
B0   =    115.952318
E0   =  -1705.738844
R0   =      0.019532
cs0  =    498.360422

if which_calculation == 'elastic_constants' and which_mode == 'generate':
    ti_hcp.elastic_constants_batch_generate(sws=sws0, ca=ca0)
if which_calculation == 'elastic_constants' and which_mode == 'analyze':
    ti_hcp.elastic_constants_analyze(sws=sws0,bmod=B0, ca=ca0, R=R0, cs=cs0)

