#!/usr/bin/python

# An example script showing how to automatically calculate
# the elastic constants of bulk bcc Fe using pyEMTO.

import pyemto

emtopath = "/wrk/hpleva/pyEMTO_examples/fe_elastic_constants" # Define a folder for our KGRN and KFCD input and output files.
latpath  = "/wrk/hpleva/structures"                           # Define a folder where the BMDL, KSTR and KFCD input and output files
                                                              # will be located.

fe = pyemto.System(folder=emtopath)  # Create a new instance of the pyemto System-class.

# Let's calculate the elastic constants of Fe.

fe.bulk(lat     = 'bcc',   # We want to use the bcc structure.
        latpath = latpath, # Path to the folder where the structure files are located.
        afm     = 'F',     # We want to do a ferromagnetic calculation.
        atoms   = ['Fe'],  # A list of atoms.
        splts   = [2.0],   # A list of magnetic splittings.
        expan   = 'M',     # We want to use the double-Taylor expansion.
        sofc    = 'Y',     # We want to use soft-core approximation.
        xc      = 'PBE',   # We want to use the PBE functional.
        amix    = 0.05)    # Mixing parameter.  

sws0 = 2.64 # Eq. WS-radius that we computed previously.
B0   = 195  # Eq. bulk modulus.

# Generate all the necessary input files with this function:
#fe.elastic_constants_batch_generate(sws=sws0)

#The newly created batch scripts are then submitted by hand.

# Analyze the results using this function once all the calculations
# have finished:
#fe.elastic_constants_analyze(sws=sws,bmod=B0)

# This function combines the features of the previous two:
fe.elastic_constants_batch_calculate(sws=sws0,bmod=B0)
