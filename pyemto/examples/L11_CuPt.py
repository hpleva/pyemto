from emto_input_generator import *
import numpy as np

folder = os.getcwd()           # Get current working directory.
emtopath = folder+"/L11_CuPt"  # Folder where the calculations will be performed.
latpath = emtopath

# L11 CuPt
prims = np.array([[1.0,0.5,0.5],
                   [0.5,1.0,0.5],
                   [0.5,0.5,1.0]])

basis = np.array([[0.0,0.0,0.0],
                  [0.5,0.5,0.5]])

species = ["Cu","Pt"]
species_cpa = [["cu","pt"],["cu","pt"]]

input_creator = EMTO(folder=emtopath)

input_creator.init_structure(latpath=latpath,
                             prims=prims,
                             basis=basis,
                             atoms=species,
                             latname='L11')

input_creator.init_bulk(atoms_cpa=species_cpa)

sws_range = np.linspace(2,3,6)

#input_creator.write_bmdl_kstr_shape_input()
#input_creator.write_kgrn_kfcd_input()

input_creator.write_kgrn_kfcd_swsrange(sws=sws_range)

#input_creator.draw_structure('standard_conv')
