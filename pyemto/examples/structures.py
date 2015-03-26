#!/usr/bin/python

# Test the functionality of the Bmdl, Kstr and Shape classes

import pyemto
import os
import numpy as np

# It is recommended to use absolute paths
folder = os.getcwd()
latpath = folder
emtopath = folder

structure=pyemto.System(folder=emtopath)

# Create BMDL, KSTR and SHAPE files and run them

#deltas = np.linspace(0.0,0.05,6)

structure.lattice.set_values(latpath=latpath,kappaw=[0.0,-20.0]) # Double Taylor-expansion

#dmaxs = [2.52,2.49,2.455,2.43,2.4,2.4] # Orthorhombic distortion for c/a = 1.603
#dmaxs = [2.51,2.51,2.51,2.49,2.51,2.49] # Monoclinic distortion for c/a = 1.603

#jobnames = ['bcco0','bcco1','bcco2','bcco3','bcco4','bcco5'] # bcc ortho
#jobnames = ['bccm0','bccm1','bccm2','bccm3','bccm4','bccm5'] # bcc mono

#jobnames = ['fcco0','fcco1','fcco2','fcco3','fcco4','fcco5'] # fcc ortho
#jobnames = ['fccm0','fccm1','fccm2','fccm3','fccm4','fccm5'] # fcc mono

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


# bcc ortho
"""
for i in range(6):
    structure.lattice.distortion(lat='bcc',dist='ortho',index=i)
    
    structure.lattice.bmdl.write_input_file(folder=latpath)
    structure.lattice.kstr.write_input_file(folder=latpath)
    structure.lattice.shape.write_input_file(folder=latpath)
    structure.lattice.batch.write_input_file(folder=latpath)


# bcc mono

for i in range(6):
    structure.lattice.distortion(lat='bcc',dist='mono',index=i)
    
    structure.lattice.bmdl.write_input_file(folder=latpath)
    structure.lattice.kstr.write_input_file(folder=latpath)
    structure.lattice.shape.write_input_file(folder=latpath)
    structure.lattice.batch.write_input_file(folder=latpath)
"""

# fcc ortho
"""
for i in range(6):
    structure.lattice.distortion(lat='fcc',dist='ortho',index=i)
    
    structure.lattice.bmdl.write_input_file(folder=latpath)
    structure.lattice.kstr.write_input_file(folder=latpath)
    structure.lattice.shape.write_input_file(folder=latpath)
    structure.lattice.batch.write_input_file(folder=latpath)


# fcc mono

for i in range(6):
    structure.lattice.distortion(lat='fcc',dist='mono',index=i)
    
    structure.lattice.bmdl.write_input_file(folder=latpath)
    structure.lattice.kstr.write_input_file(folder=latpath)
    structure.lattice.shape.write_input_file(folder=latpath)
    structure.lattice.batch.write_input_file(folder=latpath)
"""
