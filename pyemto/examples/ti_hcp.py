#!/usr/bin/python

# This script will automatically accomplish:

#  1. Calculate the equilibrium volume and bulk modulus of hcp Ti.
#  2. Calculate the elastic constants of hcp Ti.

import pyemto
import os
import numpy as np

# It is recommended to always use absolute paths
folder = os.getcwd()                # Get current working directory.
latpath = "/wrk/hpleva/structures"  # Folder where the structure files are located.
emtopath = folder+"/ti_hcp"         # Folder where the calculation take place.

ti_hcp=pyemto.System(folder=emtopath)

# Initialize the bulk system using the bulk() function:
ti_hcp.bulk(lat='hcp',
            latpath=latpath,
            atoms=['Ti'],
            sws=3.0,
            amix=0.02,
            efmix=0.9,
            expan='M',
            sofc='Y',
            xc='P07',           # Use PBEsol
            nky=31,             # k-points
            nkz=19,             # k-points
            runtime='24:00:00') # Allow large enough timelimit for SLURM   

sws = np.linspace(2.9,3.1,7) # A list of 7 different volumes from 2.9 to 3.1 


sws0,ca0,B0,e0,R0,cs0 = ti_hcp.lattice_constants_batch_calculate(sws=sws)
ti_hcp.elastic_constants_batch_calculate(sws=sws0,bmod=B0,ca=ca0,R=R0,cs=cs0)

# If the batch jobs are submitted by hand use these functions.

# To evaluate the results, comment out the _generate functions
# and uncomment the _analyze functions.

ti_hcp.lattice_constants_batch_generate(sws=sws)
#ti_hcp.lattice_constants_analyze(sws=sws)

# Results. These are inputed to the elastic_constants functions.
#sws0 =      3.002260
#ca0  =      1.610122
#B0   =    115.952318
#E0   =  -1705.738844
#R0   =      0.019532
#cs0  =    498.360422

#ti_hcp.elastic_constants_batch_generate(sws=sws0,ca=ca0)
#ti_hcp.elastic_constants_analyze(sws=sws0,bmod=B0,ca=ca0,R=R0,cs=cs0)


