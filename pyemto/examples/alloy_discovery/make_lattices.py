#!/usr/bin/python

# Make inputs for  the functionality of the Bmdl, Kstr and Shape programs
# for alloy discovery script

import pyemto
import os
import numpy as np

# It is recommended to use absolute paths
folder = os.getcwd()
latpath = folder
emtopath = folder



# Create BMDL, KSTR and SHAPE files
"""
#for HCP
structure=pyemto.System(folder=emtopath)
#structure.lattice.set_values(latpath=latpath,kappaw=[0.0,-20.0]) # Double Taylor-expansion

# Create hcp structure files
cas = np.linspace(1.50,1.7,7) # 7 different equally spaced c/a values from 1.5 to 1.7
latnames = ['hcp_ca1','hcp_ca2','hcp_ca3','hcp_ca4','hcp_ca5','hcp_ca6','hcp_ca7']
dmaxs = [2.39196429,2.41107143,2.43017857,2.44928571,2.46839286,2.4875,2.50660714]

for i in range(len(latnames)):
    structure.lattice.set_values(jobname=latnames[i],latpath=latpath,
                                 lat='hcp',kappaw=[0.0,-20.0],msgl=0,ca=cas[i],
                                 dmax=dmaxs[i])

    structure.lattice.bmdl.write_input_file(folder=latpath)
    structure.lattice.kstr.write_input_file(folder=latpath)
    structure.lattice.shape.write_input_file(folder=latpath)
    structure.lattice.batch.write_input_file(folder=latpath)

"""
# for BCC
structure=pyemto.System(folder=emtopath)
structure.lattice.set_values(jobname='bcc',latpath=latpath,lat='bcc',
                             kappaw=[0.0,-20.0],msgl=0,dmax=2.2)
    
structure.lattice.bmdl.write_input_file(folder=latpath)
structure.lattice.kstr.write_input_file(folder=latpath)
structure.lattice.shape.write_input_file(folder=latpath)
structure.lattice.batch.write_input_file(folder=latpath)

# for FCC
structure=pyemto.System(folder=emtopath)
structure.lattice.set_values(jobname='fcc',latpath=latpath,lat='fcc',
                             kappaw=[0.0,-20.0],msgl=0,dmax=1.7)
    
structure.lattice.bmdl.write_input_file(folder=latpath)
structure.lattice.kstr.write_input_file(folder=latpath)
structure.lattice.shape.write_input_file(folder=latpath)
structure.lattice.batch.write_input_file(folder=latpath)
