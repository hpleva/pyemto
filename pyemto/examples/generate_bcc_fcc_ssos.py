from emto_input_generator import *
import numpy as np

folder = os.getcwd()               # Get current working directory.
emtopath = folder+"/bcc_fcc_ssos"  # Folder where the calculations will be performed.
latpath = emtopath

# bcc SSOS-1
prims = np.array([[0.0,2.0,3.0],
                  [-0.5,2.5,2.5],
                  [-0.5,-0.5,0.5]])

basis = np.array([[0.0,0.0,0.0],
                  [0.6,0.6,0.4],
                  [0.4,0.4,0.6],
                  [0.2,0.2,0.8],
                  [0.8,0.8,0.2]])

#species = ['A1+','A2+','A3+','A4+','A5+']
species_cpa = ["Co","Cr","Fe","Mn","Ni"]

input_creator = EMTO(folder=emtopath)

input_creator.init_structure(latpath=latpath,
                             prims=prims,
                             basis=basis,
                             latname='ssos_bcc_1')

"""
input_creator.init_bulk(atoms_cpa=species_cpa)

sws_range = np.linspace(2,3,6)

input_creator.write_bmdl_kstr_shape_input()
#input_creator.write_kgrn_kfcd_input()
input_creator.write_kgrn_kfcd_swsrange(sws=sws_range)

input_creator.draw_structure('output')
"""
