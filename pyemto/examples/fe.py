#!/usr/bin/python

# An example script of how to automatically
# calculate the equilibrium volume and
# bulk modulus of bulk bcc Fe using pyEMTO.

import os
import pyemto

emtopath = os.getcwd() + "/fe" # Define a folder for our KGRN and KFCD input and output files.
latpath  = os.getcwd() + "/fe" # Define a folder where the BMDL, KSTR and KFCD input and output files
                               # will be located.

fe = pyemto.System(folder=emtopath)  # Create a new instance of the pyemto System-class.

# Let's calculate the equilibrium volume and bulk modulus of Fe.

# Initialize the system using the bulk.() function:

fe.bulk(lat     = 'bcc',   # We want to use the bcc structure.
        latpath = latpath, # Path to the folder where the structure files are located.
        afm     = 'F',     # We want to do a ferromagnetic calculation.
        atoms   = ['Fe'],  # A list of atoms.
        splts   = [2.0],   # A list of magnetic splittings.
        expan   = 'M',     # We want to use the double-Taylor expansion.
        sofc    = 'Y',     # We want to use soft-core approximation.
        xc      = 'PBE',   # We want to use the PBE functional.
        amix    = 0.05,    # Density mixing.
        nky     = 21)      # Number of k-points.

sws = [2.6,2.62,2.64,2.66,2.68,2.70] # A list of Wigner-Seitz radii

# Generate all the necessary input files with this function:
fe.lattice_constants_batch_generate(sws=sws)

#The newly created batch scripts are then submitted by hand.

# Analyze the results using this function once all the calculations
# have finished:
#fe.lattice_constants_analyze(sws=sws)

# This function combines the features of the previous two:
#fe.lattice_constants_batch_calculate(sws=sws)
