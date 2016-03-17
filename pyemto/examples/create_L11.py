from emto_input_generator import *

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

input_creator = Create_EMTO_input(folder=emtopath,
                                  latpath=latpath,
                                  prims=prims,
                                  basis=basis,
                                  species=species,
                                  latname='L11')

input_creator.create_structure_input()
input_creator.create_KGRN_KFCD_input(sws=3.0)
#input_creator.draw_structure('standard_conv')
