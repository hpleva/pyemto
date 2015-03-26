#!/usr/bin/python

# Calculate the equilibrium volume and bulk modulus of
# fcc CoCrFeMnNi high-entropy alloy using DLM.

# !!! NOTE !!!
# This script DOES NOT automatically start running
# the batch scripts. It only generates the input files
# and batch scripts which the user runs by themselves
# !!! NOTE !!!

import pyemto
import os
import numpy as np

# It is recommended to always use absolute paths
folder = os.getcwd()                 # Get current working directory.
latpath = "/wrk/hpleva/structures"   # Folder where the structure output files are.
emtopath = folder+"/cocrfemnni_fcc"  # Folder where the calculations will be performed.

cocrfemnni=pyemto.System(folder=emtopath)

sws = np.linspace(2.50,2.70,11) # 11 different volumes from 2.5 Bohr to 2.7 Bohr

# Set KGRN and KFCD values using a for loop.
# Use write_input_file functions to write input files to disk:

for i in range(len(sws)):
    cocrfemnni.bulk(lat='fcc',
                    jobname='cocrfemnni',
                    latpath=latpath,
                    atoms=['Co','Co','Cr','Cr','Fe','Fe','Mn','Mn','Ni','Ni'],
                    splts=[-1.0,1.0,-1.0,1.0,-1.0,1.0,-1.0,1.0,-1.0,1.0],
                    sws=sws[i],
                    amix=0.02,
                    efmix=0.9,
                    expan='M',
                    sofc='Y',
                    afm='M',         # Fixed-spin DLM calculation.
                    iex=7,           # We want to use self-consistent GGA (PBE).
                    nz2=16,
                    tole=1.0E-8,
                    ncpa=10,
                    nky=21,
                    tfermi=5000,
                    dx=0.015,        # Dirac equation parameters
                    dirac_np=1001,   # Dirac equation parameters
                    nes=50,          # Dirac equation parameters
                    dirac_niter=500) # Dirac equation parameters

    cocrfemnni.emto.kgrn.write_input_file(folder=emtopath)
    cocrfemnni.emto.kfcd.write_input_file(folder=emtopath)
    cocrfemnni.emto.batch.write_input_file(folder=emtopath)
