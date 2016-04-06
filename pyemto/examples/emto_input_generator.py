#!/usr/bin/python
import numpy as np
from pymatgen import Lattice, Structure
from pymatgen.vis.structure_vtk import StructureVis
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.structure_matcher import StructureMatcher

import os

# Help python to find the pyemto folder
import sys
sys.path.insert(0, "/home/henrik/Desktop/pyemto_temp/pyemto/")

import pyemto
import pyemto.common.common as common

class EMTO:
    """This class can be used to create EMTO input files from
       an arbitrary structure. What is needed as input:
       -primitive lattice vectors,
       -basis vectors,
       -list of atomic species that occupy the basis sites.
    """

    def __init__(self,folder=None,EMTOdir=None):
        """ """
        if folder is None:
            self.folder = os.getcwd()
        else:
            self.folder = folder
        if EMTOdir is None:
            self.EMTOdir = '/home/henrik/local_emto_stuff/openmp-stable-cmake'    
        else:
            self.EMTOdir = EMTOdir

        self.sg2ibz = {1:14, 2:14, 3:12, 4:12, 5:13, 6:12, 7:12, 8:13, 9:13, 10:12,
                       11:12, 12:13, 13:12, 14:12, 15:13, 16:8, 17:8, 18:8, 19:8, 20:9,
                       21:9, 22:11, 23:10, 24:10, 25:8, 26:8, 27:8, 28:8, 29:8, 30:8,
                       31:8, 32:8, 33:8, 34:8, 35:9, 36:9, 37:9, 38:9, 39:9, 40:9,
                       41:9, 42:11, 43:11, 44:10, 45:10, 46:10, 47:8, 48:8, 49:8, 50:8,
                       51:8, 52:8, 53:8, 54:8, 55:8, 56:8, 57:8, 58:8, 59:8, 60:8,
                       61:8, 62:8, 63:9, 64:9, 65:9, 66:9, 67:9, 68:9, 69:11, 70:11,
                       71:10, 72:10, 73:10, 74:10, 75:5, 76:5, 77:5, 78:5, 79:6, 80:6,
                       81:5, 82:6, 83:5, 84:5, 85:5, 86:5, 87:6, 88:6, 89:5, 90:5,
                       91:5, 92:5, 93:5, 94:5, 95:5, 96:5, 97:6, 98:6, 99:5, 100:5,
                       101:5, 102:5, 103:5, 104:5, 105:5, 106:5, 107:6, 108:6, 109:6, 110:6,
                       111:5, 112:5, 113:5, 114:5, 115:5, 116:5, 117:5, 118:5, 119:6, 120:6,
                       121:6, 122:6, 123:5, 124:5, 125:5, 126:5, 127:5, 128:5, 129:5, 130:5,
                       131:5, 132:5, 133:5, 134:5, 135:5, 136:5, 137:5, 138:5, 139:6, 140:6,
                       141:6, 142:6, 143:4, 144:4, 145:4, 146:7, 147:4, 148:7, 149:4, 150:4,
                       151:4, 152:4, 153:4, 154:4, 155:7, 156:4, 157:4, 158:4, 159:4, 160:7,
                       161:7, 162:4, 163:4, 164:4, 165:4, 166:7, 167:7, 168:4, 169:4, 170:4,
                       171:4, 172:4, 173:4, 174:4, 175:4, 176:4, 177:4, 178:4, 179:4, 180:4,
                       181:4, 182:4, 183:4, 184:4, 185:4, 186:4, 187:4, 188:4, 189:4, 190:4,
                       191:4, 192:4, 193:4, 194:4, 195:1, 196:2, 197:3, 198:1, 199:3, 200:1,
                       201:1, 202:2, 203:2, 204:3, 205:1, 206:3, 207:1, 208:1, 209:2, 210:2,
                       211:3, 212:1, 213:1, 214:3, 215:1, 216:2, 217:3, 218:1, 219:2, 220:3,
                       221:1, 222:1, 223:1, 224:1, 225:2, 226:2, 227:2, 228:2, 229:3, 230:3}

        self.sg2bl = {1:'simple triclinic', 2:'simple triclinic', 3:'simple monoclinic', 4:'simple monoclinic', 5:'base-centered monoclinic',
                      6:'simple monoclinic', 7:'simple monoclinic', 8:'base-centered monoclinic', 9:'base-centered monoclinic', 10:'simple monoclinic',
                      11:'simple monoclinic', 12:'base-centered monoclinic', 13:'simple monoclinic', 14:'simple monoclinic', 15:'base-centered monoclinic',
                      16:'simple orthorhombic', 17:'simple orthorhombic', 18:'simple orthorhombic', 19:'simple orthorhombic', 20:'base-centered orthorhombic',
                      21:'base-centered orthorhombic', 22:'face-centered orthorhombic', 23:'body-centered orthorhombic', 24:'body-centered orthorhombic', 25:'simple orthorhombic',
                      26:'simple orthorhombic', 27:'simple orthorhombic', 28:'simple orthorhombic', 29:'simple orthorhombic', 30:'simple orthorhombic',
                      31:'simple orthorhombic', 32:'simple orthorhombic', 33:'simple orthorhombic', 34:'simple orthorhombic', 35:'base-centered orthorhombic',
                      36:'base-centered orthorhombic', 37:'base-centered orthorhombic', 38:'base-centered orthorhombic', 39:'base-centered orthorhombic', 40:'base-centered orthorhombic',
                      41:'base-centered orthorhombic', 42:'face-centered orthorhombic', 43:'face-centered orthorhombic', 44:'body-centered orthorhombic', 45:'body-centered orthorhombic',
                      46:'body-centered orthorhombic', 47:'simple orthorhombic', 48:'simple orthorhombic', 49:'simple orthorhombic', 50:'simple orthorhombic',
                      51:'simple orthorhombic', 52:'simple orthorhombic', 53:'simple orthorhombic', 54:'simple orthorhombic', 55:'simple orthorhombic',
                      56:'simple orthorhombic', 57:'simple orthorhombic', 58:'simple orthorhombic', 59:'simple orthorhombic', 60:'simple orthorhombic',
                      61:'simple orthorhombic', 62:'simple orthorhombic', 63:'base-centered orthorhombic', 64:'base-centered orthorhombic', 65:'base-centered orthorhombic',
                      66:'base-centered orthorhombic', 67:'base-centered orthorhombic', 68:'base-centered orthorhombic', 69:'face-centered orthorhombic', 70:'face-centered orthorhombic',
                      71:'body-centered orthorhombic', 72:'body-centered orthorhombic', 73:'body-centered orthorhombic', 74:'body-centered orthorhombic', 75:'simple tetragonal',
                      76:'simple tetragonal', 77:'simple tetragonal', 78:'simple tetragonal', 79:'body-centered tetragonal', 80:'body-centered tetragonal',
                      81:'simple tetragonal', 82:'body-centered tetragonal', 83:'simple tetragonal', 84:'simple tetragonal', 85:'simple tetragonal',
                      86:'simple tetragonal', 87:'body-centered tetragonal', 88:'body-centered tetragonal', 89:'simple tetragonal', 90:'simple tetragonal',
                      91:'simple tetragonal', 92:'simple tetragonal', 93:'simple tetragonal', 94:'simple tetragonal', 95:'simple tetragonal',
                      96:'simple tetragonal', 97:'body-centered tetragonal', 98:'body-centered tetragonal', 99:'simple tetragonal', 100:'simple tetragonal',
                      101:'simple tetragonal', 102:'simple tetragonal', 103:'simple tetragonal', 104:'simple tetragonal', 105:'simple tetragonal',
                      106:'simple tetragonal', 107:'body-centered tetragonal', 108:'body-centered tetragonal', 109:'body-centered tetragonal', 110:'body-centered tetragonal',
                      111:'simple tetragonal', 112:'simple tetragonal', 113:'simple tetragonal', 114:'simple tetragonal', 115:'simple tetragonal',
                      116:'simple tetragonal', 117:'simple tetragonal', 118:'simple tetragonal', 119:'body-centered tetragonal', 120:'body-centered tetragonal',
                      121:'body-centered tetragonal', 122:'body-centered tetragonal', 123:'simple tetragonal', 124:'simple tetragonal', 125:'simple tetragonal',
                      126:'simple tetragonal', 127:'simple tetragonal', 128:'simple tetragonal', 129:'simple tetragonal', 130:'simple tetragonal',
                      131:'simple tetragonal', 132:'simple tetragonal', 133:'simple tetragonal', 134:'simple tetragonal', 135:'simple tetragonal',
                      136:'simple tetragonal', 137:'simple tetragonal', 138:'simple tetragonal', 139:'body-centered tetragonal', 140:'body-centered tetragonal',
                      141:'body-centered tetragonal', 142:'body-centered tetragonal', 143:'hexagonal', 144:'hexagonal', 145:'hexagonal',
                      146:'rhombohedral', 147:'hexagonal', 148:'rhombohedral', 149:'hexagonal', 150:'hexagonal',
                      151:'hexagonal', 152:'hexagonal', 153:'hexagonal', 154:'hexagonal', 155:'rhombohedral',
                      156:'hexagonal', 157:'hexagonal', 158:'hexagonal', 159:'hexagonal', 160:'rhombohedral',
                      161:'rhombohedral', 162:'hexagonal', 163:'hexagonal', 164:'hexagonal', 165:'hexagonal',
                      166:'rhombohedral', 167:'rhombohedral', 168:'hexagonal', 169:'hexagonal', 170:'hexagonal',
                      171:'hexagonal', 172:'hexagonal', 173:'hexagonal', 174:'hexagonal', 175:'hexagonal',
                      176:'hexagonal', 177:'hexagonal', 178:'hexagonal', 179:'hexagonal', 180:'hexagonal',
                      181:'hexagonal', 182:'hexagonal', 183:'hexagonal', 184:'hexagonal', 185:'hexagonal',
                      186:'hexagonal', 187:'hexagonal', 188:'hexagonal', 189:'hexagonal', 190:'hexagonal',
                      191:'hexagonal', 192:'hexagonal', 193:'hexagonal', 194:'hexagonal', 195:'simple cubic',
                      196:'face-centered cubic', 197:'body-centered cubic', 198:'simple cubic', 199:'body-centered cubic', 200:'simple cubic',
                      201:'simple cubic', 202:'face-centered cubic', 203:'face-centered cubic', 204:'body-centered cubic', 205:'simple cubic',
                      206:'body-centered cubic', 207:'simple cubic', 208:'simple cubic', 209:'face-centered cubic', 210:'face-centered cubic',
                      211:'body-centered cubic', 212:'simple cubic', 213:'simple cubic', 214:'body-centered cubic', 215:'simple cubic',
                      216:'face-centered cubic', 217:'body-centered cubic', 218:'simple cubic', 219:'face-centered cubic', 220:'body-centered cubic',
                      221:'simple cubic', 222:'simple cubic', 223:'simple cubic', 224:'simple cubic', 225:'face-centered cubic',
                      226:'face-centered cubic', 227:'face-centered cubic', 228:'face-centered cubic', 229:'body-centered cubic', 230:'body-centered cubic'}

        # BMDL, KSTR, SHAPE, KGRN and KFCD class instances
        self.input_system = pyemto.System(folder=self.folder,EMTOdir=self.EMTOdir)
        return

    def calc_ws_radius(self,struct):
        bohr2angst = 0.52917721
        vol_unit = struct.volume/struct.num_sites
        sws = (3*vol_unit/4.0/np.pi)**(1.0/3)/bohr2angst
        return sws

    def make_basis_array(self,struct):
        """Returns a 2D numpy array of the basis atom coordinates
           in !!Cartesian!! coordinates.
        """
        len_basis = struct.num_sites
        emto_basis = np.zeros((len_basis,3))
        for i in range(len_basis):
            emto_basis[i,:] = struct.sites[i].coords
        return emto_basis

    def make_sites_array(self,struct):
        len_basis = struct.num_sites
        emto_sites = []
        for i in range(len_basis):
            emto_sites.append(struct.sites[i].species_string)
        return emto_sites

    def make_sites_index(self,species):
        len_array = len(species)
        index_array = np.zeros(len_array,dtype='int32')
        for i in range(len_array):
            if species[i] == 'A+':
                index_array[i] = 1
            else:
                index_array[i] = int(species[i][1:-1])
        return index_array

    def init_structure(self,prims=None,basis=None,latpath=None,
                       coords_are_cartesian=False,latname=None,kappaw=None):
        if prims is None:
            sys.exit('EMTO.init_structure(): \'prims\' has to be given!')
        if basis is None:
            sys.exit('EMTO.init_structure(): \'basis\' has to be given!')
        #if atoms is None:
        #    sys.exit('EMTO.init_structure(): \'atoms\' has to be given!')
        if latpath is None:
            self.latpath = os.getcwd()
        else:
            self.latpath = latpath
        if latname is None:
            self.latname = 'structure'
        else:
            self.latname = latname

        self.prims   = np.asarray(prims)
        self.basis   = np.asarray(basis)
        self.len_basis = len(self.basis[:,0])
        # Construct species array. We use A+, A2+ etc.
        # so that the order of basis vectors does not
        # change when standard structure is created.
        atoms = []
        for i in range(self.len_basis):
            atoms.append('A{0}+'.format(i+1))
        self.species = np.asarray(atoms)
        self.coords_are_cartesian = coords_are_cartesian
        self.ibz = None
        if kappaw == None:
            self.kappaw = [0.0]
        else:
            self.kappaw = kappaw
        #
        self.pmg_input_lattice = Lattice(self.prims)
        self.pmg_input_struct  = Structure(self.pmg_input_lattice, self.species, self.basis, 
                                           coords_are_cartesian=self.coords_are_cartesian)
        self.sws = self.calc_ws_radius(self.pmg_input_struct)

        #
        self.finder = SpacegroupAnalyzer(self.pmg_input_struct)
        self.stm = StructureMatcher(ltol=0.001,stol=0.001,angle_tol=1.0)
        #
        print("Input structure information:")
        print(self.pmg_input_struct)
        print("Volume: ",self.pmg_input_struct.volume)
        print("Lattice vectors:")
        print(self.pmg_input_struct.lattice.matrix)
        print("")
        #
        self.conv_struct = self.finder.get_conventional_standard_structure(international_monoclinic=False)
        self.prim_struct = self.finder.get_primitive_standard_structure(international_monoclinic=False)
        #print('self.prim_struct:')
        #print(self.prim_struct)
        #print("")
        # Sort the sequence of basis vectors to be the same as the input structure.
        # I want to do this because when the conventional (and EMTO) cell is created, sites end up in random order
        # and the user might not realize that and end up putting species in wrong sites.
        #species_input_index  = self.make_sites_index(self.make_sites_array(self.pmg_input_struct))
        #species_output_index = self.make_sites_index(self.make_sites_array(self.prim_struct))
        #sorting_array = np.argsort(species_output_index)
        #print(species_input_index)
        #print(species_output_index)
        #print(sorting_array)
        #prim_struct_temp = self.prim_struct.copy()
        #print('prim_struct_temp:')
        #print(prim_struct_temp)
        #print("")
        #prim_struct_temp[0] = prim_struct_temp[2]
        #print('self.prim_struct:')
        #print(self.prim_struct)
        #print("")
        #print('prim_struct_temp:')
        #print(prim_struct_temp)
        #print("")
        #for i in range(len(sorting_array)):
        #    self.prim_struct[i] = prim_struct_temp[sorting_array[i]]
        #
        self.finder_prim = SpacegroupAnalyzer(self.prim_struct)
        self.finder_space = self.finder_prim.get_spacegroup_number()
        self.ibz_string = self.sg2bl[self.finder_space]
        self.ibz = self.sg2ibz[self.finder_space]
        #
        print("Detected standard conventional structure:")
        print(self.conv_struct)
        print("Volume: ",self.conv_struct.volume)
        print("Lattice vectors:")
        print(self.conv_struct.lattice.matrix)
        print("")
        print("Detected standard primitive structure:")
        print(self.prim_struct)
        print("Volume: ",self.prim_struct.volume)
        print("Lattice vectors:")
        print(self.prim_struct.lattice.matrix)
        print("")
        print("The spacegroup symbol of your structure: {}".format(self.finder_prim.get_spacegroup_symbol()))
        print("The spacegroup number of your structure: {}".format(self.finder_prim.get_spacegroup_number()))
        print("The crystal system of your structure   : {}".format(self.finder_prim.get_crystal_system()))
        print("The lattice type of your structure     : {}".format(self.finder_prim.get_lattice_type()))
        print("The point group of your structure      : {}".format(self.finder_prim.get_point_group()))
        print("The Bravais lattice of your structure  : {}".format(self.ibz_string))
        print("Number of basis atoms                  : {}".format(self.prim_struct.num_sites))
        print("")
        #
        if self.sg2ibz[self.finder_space] == 1:
            self.primaa = self.prim_struct.lattice.matrix[0,:]
            self.primbb = self.prim_struct.lattice.matrix[1,:]
            self.primcc = self.prim_struct.lattice.matrix[2,:]
            norm_tmp = np.linalg.norm(self.primaa)
            self.output_prima = self.primaa/norm_tmp
            self.output_primb = self.primbb/norm_tmp
            self.output_primc = self.primcc/norm_tmp
            self.output_basis = self.make_basis_array(self.prim_struct)
            # Apply transformation on the basis atoms
            self.output_basis = self.output_basis/norm_tmp
            self.output_boa = 0.0
            self.output_coa = 0.0
            self.output_alpha = 0.0
            self.output_beta = 0.0
            self.output_gamma = 0.0
        elif self.sg2ibz[self.finder_space] == 2:
            self.primaa = self.prim_struct.lattice.matrix[0,:]
            self.primbb = self.prim_struct.lattice.matrix[1,:]
            self.primcc = self.prim_struct.lattice.matrix[2,:]
            norm_tmp = 2*self.primaa[1]
            self.output_prima = self.primcc/norm_tmp
            self.output_primb = self.primaa/norm_tmp
            self.output_primc = self.primbb/norm_tmp
            self.output_basis = self.make_basis_array(self.prim_struct)
            # Apply transformation on the basis atoms 
            self.output_basis = self.output_basis/norm_tmp
            self.output_boa = 0.0
            self.output_coa = 0.0
            self.output_alpha = 0.0
            self.output_beta = 0.0
            self.output_gamma = 0.0
        elif self.sg2ibz[self.finder_space] == 3:
            self.primaa = self.prim_struct.lattice.matrix[0,:]
            self.primbb = self.prim_struct.lattice.matrix[1,:]
            self.primcc = self.prim_struct.lattice.matrix[2,:]
            norm_tmp = 2*self.primaa[1]
            self.output_prima = self.primcc/norm_tmp
            self.output_primb = self.primaa/norm_tmp
            self.output_primc = self.primbb/norm_tmp
            self.output_basis = self.make_basis_array(self.prim_struct)
            # Apply transformation on the basis atoms
            self.output_basis = self.output_basis/norm_tmp
            self.output_boa = 0.0
            self.output_coa = 0.0
            self.output_alpha = 0.0
            self.output_beta = 0.0
            self.output_gamma = 0.0
        elif self.sg2ibz[self.finder_space] == 4:
            # Hexagonal, have to rotate the vectors...
            pass
        elif self.sg2ibz[self.finder_space] == 5:
            self.primaa = self.prim_struct.lattice.matrix[0,:]
            self.primbb = self.prim_struct.lattice.matrix[1,:]
            self.primcc = self.prim_struct.lattice.matrix[2,:]
            norm_tmp = self.primaa[0]
            self.output_prima = self.primaa/norm_tmp
            self.output_primb = self.primbb/norm_tmp
            self.output_primc = self.primcc/norm_tmp
            self.output_basis = self.make_basis_array(self.prim_struct)
            # Apply transformation on the basis atoms
            self.output_basis = self.output_basis/norm_tmp
            self.output_boa = 0.0
            self.output_coa = self.output_primc[2]
            self.output_alpha = 0.0
            self.output_beta = 0.0
            self.output_gamma = 0.0
        elif self.sg2ibz[self.finder_space] == 6:
            self.primaa = self.prim_struct.lattice.matrix[0,:]
            self.primbb = self.prim_struct.lattice.matrix[1,:]
            self.primcc = self.prim_struct.lattice.matrix[2,:]
            norm_tmp = 2*self.primbb[0]
            self.output_prima = self.primaa/norm_tmp
            self.output_primb = self.primbb/norm_tmp
            self.output_primc = self.primcc/norm_tmp
            self.output_basis = self.make_basis_array(self.prim_struct)
            # Apply transformation on the basis atoms
            self.output_basis = self.output_basis/norm_tmp
            self.output_boa = 0.0
            self.output_coa = 2*self.output_prima[2]
            self.output_alpha = 0.0
            self.output_beta = 0.0
            self.output_gamma = 0.0            
        elif self.sg2ibz[self.finder_space] == 7:
            from pyemto.utilities.utils import rotation_matrix
            self.primaa = self.prim_struct.lattice.matrix[0,:]
            self.primbb = self.prim_struct.lattice.matrix[1,:]
            self.primcc = self.prim_struct.lattice.matrix[2,:]
            self.output_basis = self.make_basis_array(self.prim_struct)
            alpha = self.prim_struct.lattice.alpha
            kulma = np.arctan((self.primaa[0]+self.primbb[0]+self.primcc[0])/
                              (self.primaa[2]+self.primbb[2]+self.primcc[2]))
            rot1 = rotation_matrix([0.0,-1.0,0.0],kulma)
            rot2 = np.array([[-np.sqrt(3.0)/2,-0.5,0.0],
                             [0.5,-np.sqrt(3.0)/2,0.0],
                             [0.0,0.0,1.0]])
            self.output_prima = np.dot(rot2,np.dot(rot1,self.primaa))
            self.output_primb = np.dot(rot2,np.dot(rot1,self.primbb))
            self.output_primc = np.dot(rot2,np.dot(rot1,self.primcc))
            scale_a = self.output_prima[1]
            self.output_prima = self.output_prima/scale_a
            self.output_primb = self.output_primb/scale_a
            self.output_primc = self.output_primc/scale_a
            # Apply transformation on the basis atoms
            for i in range(len(self.output_basis[:,0])):
                self.output_basis[i,:] = np.dot(rot2,np.dot(rot1,self.output_basis[i,:]))/scale_a
            self.output_boa = 1.0
            self.output_coa = self.output_prima[2]
            self.output_alpha = 0.0
            self.output_beta = 0.0
            self.output_gamma = 0.0
            #print('emto_coa = ',emto_coa)
        elif self.sg2ibz[self.finder_space] == 8:
            self.primaa = self.prim_struct.lattice.matrix[0,:]
            self.primbb = self.prim_struct.lattice.matrix[1,:]
            self.primcc = self.prim_struct.lattice.matrix[2,:]
            norm_tmp = self.primaa[0]
            self.output_prima = self.primaa/norm_tmp
            self.output_primb = self.primbb/norm_tmp
            self.output_primc = self.primcc/norm_tmp
            self.output_basis = self.make_basis_array(self.prim_struct)
            # Apply transformation on the basis atoms
            self.output_basis = self.output_basis/norm_tmp
            self.output_boa = self.output_primb[1]
            self.output_coa = self.output_primc[2]
            self.output_alpha = 0.0
            self.output_beta = 0.0
            self.output_gamma = 0.0
        elif self.sg2ibz[self.finder_space] == 9:
            self.primaa = self.prim_struct.lattice.matrix[0,:]
            self.primbb = self.prim_struct.lattice.matrix[1,:]
            self.primcc = self.prim_struct.lattice.matrix[2,:]
            self.output_basis = self.make_basis_array(self.prim_struct)
            norm_tmp = 2*self.primaa[0]
            self.output_prima = self.primaa/norm_tmp
            self.output_primb = self.primbb/norm_tmp
            self.output_primc = self.primcc/norm_tmp
            # Apply transformation on the basis atoms
            self.output_basis = self.output_basis/norm_tmp
            self.output_boa = 2*self.output_primb[1]
            self.output_coa = self.output_primc[2]
            self.output_alpha = 0.0
            self.output_beta = 0.0
            self.output_gamma = 0.0
        elif self.sg2ibz[self.finder_space] == 10:
            self.primaa = self.prim_struct.lattice.matrix[0,:]
            self.primbb = self.prim_struct.lattice.matrix[1,:]
            self.primcc = self.prim_struct.lattice.matrix[2,:]
            self.output_basis = self.make_basis_array(self.prim_struct)
            norm_tmp = 2*self.primbb[0]
            self.output_prima = self.primbb/norm_tmp
            self.output_primb = self.primcc/norm_tmp
            self.output_primc = self.primaa/norm_tmp
            # Apply transformation on the basis atoms
            self.output_basis = self.output_basis/norm_tmp
            self.output_boa = 2*self.output_primb[1]
            self.output_coa = 2*self.output_primc[2]
            self.output_alpha = 0.0
            self.output_beta = 0.0
            self.output_gamma = 0.0
        elif self.sg2ibz[self.finder_space] == 11:
            self.primaa = self.prim_struct.lattice.matrix[0,:]
            self.primbb = self.prim_struct.lattice.matrix[1,:]
            self.primcc = self.prim_struct.lattice.matrix[2,:]
            self.output_basis = self.make_basis_array(self.prim_struct)
            norm_tmp = 2*self.primbb[0]
            self.output_prima = self.primbb/norm_tmp
            self.output_primb = self.primcc/norm_tmp
            self.output_primc = self.primaa/norm_tmp
            # Apply transformation on the basis atoms
            self.output_basis = self.output_basis/norm_tmp
            self.output_boa = 2*self.output_primc[1]
            self.output_coa = 2*self.output_primc[2]
            self.output_alpha = 0.0
            self.output_beta = 0.0
            self.output_gamma = 0.0
        elif self.sg2ibz[self.finder_space] == 12:
            # Simple monoclinic, need to rotate lattice vectors...
            pass
        elif self.sg2ibz[self.finder_space] == 13:
            from pymatgen.util.coord_utils import get_angle
            self.primaa = self.prim_struct.lattice.matrix[0,:]
            self.primbb = self.prim_struct.lattice.matrix[1,:]
            self.primcc = self.prim_struct.lattice.matrix[2,:]
            gamma = get_angle(self.primcc,self.primaa+self.primbb)
            rot1 = np.array([[1.0,0.0,0.0],
                             [0.0,np.cos(np.radians(180-gamma)),-np.sin(np.radians(180-gamma))],
                             [0.0,np.sin(np.radians(180-gamma)),np.cos(np.radians(180-gamma))]])
            rot2 = np.array([[0.0,0.0,1.0],
                             [0.0,1.0,0.0],
                             [-1.0,0.0,0.0]])
            bc_norm = np.linalg.norm(self.primaa+self.primbb)
            self.output_prima = np.dot(rot2,np.dot(rot1,self.primcc))/bc_norm
            self.output_primb = np.dot(rot2,np.dot(rot1,self.primaa))/bc_norm
            self.output_primc = np.dot(rot2,np.dot(rot1,self.primbb))/bc_norm
            self.output_basis = self.make_basis_array(self.prim_struct)
            # Apply transformation on the basis atoms
            for i in range(len(self.output_basis[:,0])):
                self.output_basis[i,:] = np.dot(rot2,np.dot(rot1,self.output_basis[i,:]))/bc_norm
            self.output_boa = np.abs(self.output_prima[1])
            self.output_coa = np.abs(2*self.output_primc[2])
            self.output_alpha = 0.0
            self.output_beta = 0.0
            self.output_gamma = gamma
        elif self.sg2ibz[self.finder_space] == 14:
            self.primaa = self.prim_struct.lattice.matrix[0,:]
            self.primbb = self.prim_struct.lattice.matrix[1,:]
            self.primcc = self.prim_struct.lattice.matrix[2,:]
            self.output_basis = self.make_basis_array(self.prim_struct)
            norm_tmp = self.primaa[0]
            self.output_prima = self.primaa/norm_tmp
            self.output_primb = self.primbb/norm_tmp
            self.output_primc = self.primcc/norm_tmp
            # Apply transformation on the basis atoms
            self.output_basis = self.output_basis/norm_tmp
            # This could be tested, should be OK:
            #self.output_boa = np.sqrt(self.output_primb[0]**2+self.output_primb[1]**2)
            self.output_boa = self.prim_struct.lattice.b/self.prim_struct.lattice.a
            self.output_coa = self.prim_struct.lattice.c/self.prim_struct.lattice.a
            self.output_alpha = self.prim_struct.lattice.alpha
            self.output_beta = self.prim_struct.lattice.beta
            self.output_gamma = self.prim_struct.lattice.gamma
        #
        self.output_sites = self.make_sites_array(self.prim_struct)
        self.output_lattice = Lattice(np.array([self.output_prima,self.output_primb,self.output_primc]))
        self.output_struct = Structure(self.output_lattice, self.output_sites, 
                                       self.output_basis, coords_are_cartesian=True)
        #
        # Print EMTO structure information
        print("")
        print("Generated EMTO structure:")
        print(self.output_struct)
        print("Volume: ",self.output_struct.volume)
        print("WS-rad: ",self.sws)
        print("Lattice vectors:")
        print(self.output_struct.lattice.matrix)
        print("Basis vectors:")
        for i in range(len(self.output_struct.sites)):
            print(self.output_struct.sites[i].coords)
        print("")
        # Sanity check
        print('')
        print('Sanity check:')
        print('Same structure (sites only)     ?: ',self.stm.fit_anonymous(self.pmg_input_struct,self.output_struct))
        print('Same structure (sites+chemistry)?: ',self.stm.fit(self.pmg_input_struct,self.output_struct))
        print("")
        # Generate EMTO structure input files
        self.input_system.lattice.set_values(jobname=self.latname,
                                             latpath=self.latpath,
                                             lat=common.ibz_to_lat(self.ibz),
                                             latparams=[1.0,self.output_boa,self.output_coa],
                                             latvectors=[self.output_alpha,self.output_beta,self.output_gamma],
                                             basis=self.output_basis,
                                             kappaw=self.kappaw,
                                             EMTOdir=self.EMTOdir)
        return

    def write_bmdl_kstr_shape_input(self):
        self.input_system.lattice.write_structure_input_files(folder=self.folder,jobname=self.latname)
        return

    def init_bulk(self,atoms_cpa=None,splts=None,concs=None,its=None,sws=None,**kwargs):
        """Generates and writes down to disk KGRN and KFCD input files, as well as the
           SLURM job script that is used to run the calculations.
        """
        if self.ibz == None:
            sys.exit('self.ibz == None! Run create_structure_input() to generate IBZ \n'+
                     'for your structure, if you do not know it.')
        if sws == None:
            pass
        else:
            self.sws = sws
        # Construct an index array to keep track of the number of atoms in each site.
        if atoms_cpa == None:
            sys.exit('EMTO.init_bulk(): \'atoms_cpa\' has to be given!')
            #self.KGRN_atoms = np.asarray(self.output_sites)
            #index_array = np.ones(len(self.KGRN_atoms),dtype='int32')
            #index_len = np.sum(index_array)
        else:
            index_array = np.ones(len(atoms_cpa),dtype='int32')
            for i in range(len(atoms_cpa)):
                if isinstance(atoms_cpa[i],list):
                    index_array[i] = len(atoms_cpa[i])
                else:
                    index_array[i] = 1
            index_len = np.sum(index_array)
            atoms_flat = []
            for i in range(len(atoms_cpa)):
                if isinstance(atoms_cpa[i],list):
                    for j in range(len(atoms_cpa[i])):
                        atoms_flat.append(atoms_cpa[i][j])
                else:
                    atoms_flat.append(atoms_cpa[i])
            self.KGRN_atoms = np.array(atoms_flat)
        if splts == None:
            # Assume a zero moments array
            self.KGRN_splts = np.zeros(index_len)
        else:
            splts_flat = []
            for i in range(len(splts)):
                if isinstance(splts[i],list):
                    for j in range(len(splts[i])):
                        splts_flat.append(splts[i][j])
                else:
                    splts_flat.append(splts[i])
            self.KGRN_splts = np.array(splts_flat)
        if concs == None:
            # Assume 1.0 concentration, if site is occupied by a single atom,
            # and 1.0/N(i) if site(i) is occupied by N(i) atoms (CPA).
            self.KGRN_concs = np.zeros(index_len)
            running_index = 0
            for i in range(len(index_array)):
                for j in range(index_array[i]):
                    self.KGRN_concs[running_index] = 1.0/index_array[i]
                    running_index += 1
        else:
            concs_flat = []
            for i in range(len(concs)):
                if isinstance(concs[i],list):
                    for j in range(len(concs[i])):
                        atoms_flat.append(concs[i][j])
                else:
                    atoms_flat.append(concs[i])
            self.KGRN_concs = np.array(concs_flat)
        ###
        # Construct iqs, its, and itas arrays (for the KGRN atomblock).
        self.KGRN_iqs = np.zeros(index_len, dtype='int32')
        running_index = 0
        for i in range(len(index_array)):
            for j in range(index_array[i]):
                self.KGRN_iqs[running_index] = i+1
                running_index += 1
        #
        self.KGRN_itas = np.zeros(index_len, dtype='int32')
        running_index = 0
        for i in range(len(index_array)):
            for j in range(index_array[i]):
                self.KGRN_itas[running_index] = j+1
                running_index += 1
        #
        # its = array of the indexes of sublattices (IT in input file). 
        # The concept of sublattice can be used to reduce computational load;
        # if two sites are occupied by two identical elements with
        # identical magnetic moment (different spatial rotation in terms of
        # lattice symmetry is allowed), we only need to compute one
        # instance and then copy the information to the other lattice site.
        self.KGRN_its = np.zeros(index_len, dtype='int32')
        #
        if its == None:
            # Assume distinct sublattices (for safety, not speed!??)
            self.KGRN_its = np.zeros(index_len, dtype='int32')
            running_index = 0
            for i in range(len(index_array)):
                for j in range(index_array[i]):
                    self.KGRN_its[running_index] = i+1
                    running_index += 1
        elif its == 'all_same':
            self.KGRN_its = np.ones(index_len, dtype='int32')
        else:
            its_flat = []
            for i in range(len(its)):
                if isinstance(its[i],list):
                    for j in range(len(its[i])):
                        its_flat.append(its[i][j])
                else:
                    atoms_flat.append(its[i])
            self.KGRN_its = np.array(its_flat, dtype='int32')
        #
        self.input_system.bulk_new(lat=common.ibz_to_lat(self.ibz),
                                   ibz=self.ibz,
                                   latname=self.latname,
                                   latpath=self.latpath,
                                   atoms=self.KGRN_atoms,
                                   concs=self.KGRN_concs,
                                   splts=self.KGRN_splts,
                                   iqs=self.KGRN_iqs,
                                   its=self.KGRN_its,
                                   itas=self.KGRN_itas,
                                   sws=self.sws,
                                   **kwargs)
        #
    def write_kgrn_kfcd_input(self):
        self.input_system.emto.kgrn.write_input_file(folder=self.folder)
        self.input_system.emto.kfcd.write_input_file(folder=self.folder)
        self.input_system.emto.batch.write_input_file(folder=self.folder)
        return

    def write_kgrn_kfcd_swsrange(self,sws=None):
        if sws is None:
            sys.exit('EMTO.write_KGRN_KFCD_swsrange(): An array of' +
                     ' WS-radii \'sws\' has to be given!')
        jobnames = self.input_system.lattice_constants_batch_generate(sws)
        return jobnames

    def draw_structure(self,which='input'):
        self.vis = StructureVis()
        if which == 'input':
            self.vis.set_structure(self.pmg_input_struct)
        elif which == 'output':
            self.vis.set_structure(self.output_struct)
        elif which == 'standard_conv':
            self.vis.set_structure(self.conv_struct)
        elif which == 'standard_prim':
            self.vis.set_structure(self.prim_struct)
        self.vis.show()
        return
