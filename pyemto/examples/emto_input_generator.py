import sys
import numpy as np
from itertools import combinations
from pyemto.utilities.utils import rotation_matrix
import spglib as spg

try:
    from pymatgen import Lattice, Structure
    from pymatgen.vis.structure_vtk import StructureVis
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    from pymatgen.analysis.structure_matcher import StructureMatcher
    from pymatgen.util.coord import get_angle
except ImportError:
    # pymatgen has not been installed
    raise ImportError('emto_input_generator requires pymatgen>=4.4.0 to be installed!')

import os

import pyemto
import pyemto.common.common as common

class EMTO:
    """This class can be used to create EMTO input files from
       an arbitrary structure. What is needed as input:
       -primitive lattice vectors,
       -basis vectors,
       -list of atomic species that occupy the basis sites.
    """

    def __init__(self, folder=None, EMTOdir=None):
        """ """
        if folder is None:
            self.folder = os.getcwd()
        else:
            self.folder = folder
        if EMTOdir is None:
            self.EMTOdir = '/home/EMTO'
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

        self.sg2bl = {1:'simple triclinic', 2:'simple triclinic',
                      3:'simple monoclinic', 4:'simple monoclinic',
                      5:'base-centered monoclinic', 6:'simple monoclinic',
                      7:'simple monoclinic', 8:'base-centered monoclinic',
                      9:'base-centered monoclinic', 10:'simple monoclinic',
                      11:'simple monoclinic', 12:'base-centered monoclinic',
                      13:'simple monoclinic', 14:'simple monoclinic',
                      15:'base-centered monoclinic', 16:'simple orthorhombic',
                      17:'simple orthorhombic', 18:'simple orthorhombic',
                      19:'simple orthorhombic', 20:'base-centered orthorhombic',
                      21:'base-centered orthorhombic', 22:'face-centered orthorhombic',
                      23:'body-centered orthorhombic', 24:'body-centered orthorhombic',
                      25:'simple orthorhombic', 26:'simple orthorhombic',
                      27:'simple orthorhombic', 28:'simple orthorhombic',
                      29:'simple orthorhombic', 30:'simple orthorhombic',
                      31:'simple orthorhombic', 32:'simple orthorhombic',
                      33:'simple orthorhombic', 34:'simple orthorhombic',
                      35:'base-centered orthorhombic', 36:'base-centered orthorhombic',
                      37:'base-centered orthorhombic', 38:'base-centered orthorhombic',
                      39:'base-centered orthorhombic', 40:'base-centered orthorhombic',
                      41:'base-centered orthorhombic', 42:'face-centered orthorhombic',
                      43:'face-centered orthorhombic', 44:'body-centered orthorhombic',
                      45:'body-centered orthorhombic', 46:'body-centered orthorhombic',
                      47:'simple orthorhombic', 48:'simple orthorhombic',
                      49:'simple orthorhombic', 50:'simple orthorhombic',
                      51:'simple orthorhombic', 52:'simple orthorhombic',
                      53:'simple orthorhombic', 54:'simple orthorhombic',
                      55:'simple orthorhombic', 56:'simple orthorhombic',
                      57:'simple orthorhombic', 58:'simple orthorhombic',
                      59:'simple orthorhombic', 60:'simple orthorhombic',
                      61:'simple orthorhombic', 62:'simple orthorhombic',
                      63:'base-centered orthorhombic', 64:'base-centered orthorhombic',
                      65:'base-centered orthorhombic', 66:'base-centered orthorhombic',
                      67:'base-centered orthorhombic', 68:'base-centered orthorhombic',
                      69:'face-centered orthorhombic', 70:'face-centered orthorhombic',
                      71:'body-centered orthorhombic', 72:'body-centered orthorhombic',
                      73:'body-centered orthorhombic', 74:'body-centered orthorhombic',
                      75:'simple tetragonal', 76:'simple tetragonal',
                      77:'simple tetragonal', 78:'simple tetragonal',
                      79:'body-centered tetragonal', 80:'body-centered tetragonal',
                      81:'simple tetragonal', 82:'body-centered tetragonal',
                      83:'simple tetragonal', 84:'simple tetragonal',
                      85:'simple tetragonal', 86:'simple tetragonal',
                      87:'body-centered tetragonal', 88:'body-centered tetragonal',
                      89:'simple tetragonal', 90:'simple tetragonal',
                      91:'simple tetragonal', 92:'simple tetragonal',
                      93:'simple tetragonal', 94:'simple tetragonal',
                      95:'simple tetragonal', 96:'simple tetragonal',
                      97:'body-centered tetragonal', 98:'body-centered tetragonal',
                      99:'simple tetragonal', 100:'simple tetragonal',
                      101:'simple tetragonal', 102:'simple tetragonal',
                      103:'simple tetragonal', 104:'simple tetragonal',
                      105:'simple tetragonal', 106:'simple tetragonal',
                      107:'body-centered tetragonal', 108:'body-centered tetragonal',
                      109:'body-centered tetragonal', 110:'body-centered tetragonal',
                      111:'simple tetragonal', 112:'simple tetragonal',
                      113:'simple tetragonal', 114:'simple tetragonal',
                      115:'simple tetragonal', 116:'simple tetragonal',
                      117:'simple tetragonal', 118:'simple tetragonal',
                      119:'body-centered tetragonal', 120:'body-centered tetragonal',
                      121:'body-centered tetragonal', 122:'body-centered tetragonal',
                      123:'simple tetragonal', 124:'simple tetragonal',
                      125:'simple tetragonal', 126:'simple tetragonal',
                      127:'simple tetragonal', 128:'simple tetragonal',
                      129:'simple tetragonal', 130:'simple tetragonal',
                      131:'simple tetragonal', 132:'simple tetragonal',
                      133:'simple tetragonal', 134:'simple tetragonal',
                      135:'simple tetragonal', 136:'simple tetragonal',
                      137:'simple tetragonal', 138:'simple tetragonal',
                      139:'body-centered tetragonal', 140:'body-centered tetragonal',
                      141:'body-centered tetragonal', 142:'body-centered tetragonal',
                      143:'hexagonal', 144:'hexagonal',
                      145:'hexagonal', 146:'rhombohedral',
                      147:'hexagonal', 148:'rhombohedral',
                      149:'hexagonal', 150:'hexagonal',
                      151:'hexagonal', 152:'hexagonal',
                      153:'hexagonal', 154:'hexagonal',
                      155:'rhombohedral', 156:'hexagonal',
                      157:'hexagonal', 158:'hexagonal',
                      159:'hexagonal', 160:'rhombohedral',
                      161:'rhombohedral', 162:'hexagonal',
                      163:'hexagonal', 164:'hexagonal',
                      165:'hexagonal', 166:'rhombohedral',
                      167:'rhombohedral', 168:'hexagonal',
                      169:'hexagonal', 170:'hexagonal',
                      171:'hexagonal', 172:'hexagonal',
                      173:'hexagonal', 174:'hexagonal',
                      175:'hexagonal', 176:'hexagonal',
                      177:'hexagonal', 178:'hexagonal',
                      179:'hexagonal', 180:'hexagonal',
                      181:'hexagonal', 182:'hexagonal',
                      183:'hexagonal', 184:'hexagonal',
                      185:'hexagonal', 186:'hexagonal',
                      187:'hexagonal', 188:'hexagonal',
                      189:'hexagonal', 190:'hexagonal',
                      191:'hexagonal', 192:'hexagonal',
                      193:'hexagonal', 194:'hexagonal',
                      195:'simple cubic', 196:'face-centered cubic',
                      197:'body-centered cubic', 198:'simple cubic',
                      199:'body-centered cubic', 200:'simple cubic',
                      201:'simple cubic', 202:'face-centered cubic',
                      203:'face-centered cubic', 204:'body-centered cubic',
                      205:'simple cubic', 206:'body-centered cubic',
                      207:'simple cubic', 208:'simple cubic',
                      209:'face-centered cubic', 210:'face-centered cubic',
                      211:'body-centered cubic', 212:'simple cubic',
                      213:'simple cubic', 214:'body-centered cubic',
                      215:'simple cubic', 216:'face-centered cubic',
                      217:'body-centered cubic', 218:'simple cubic',
                      219:'face-centered cubic', 220:'body-centered cubic',
                      221:'simple cubic', 222:'simple cubic',
                      223:'simple cubic', 224:'simple cubic',
                      225:'face-centered cubic', 226:'face-centered cubic',
                      227:'face-centered cubic', 228:'face-centered cubic',
                      229:'body-centered cubic', 230:'body-centered cubic'}

        # BMDL, KSTR, SHAPE, KGRN and KFCD class instances
        self.input_system = pyemto.System(folder=self.folder, EMTOdir=self.EMTOdir)
        #
        self.fit_angle_tol = 5e-6
        self.fit_norm_ratio_tol = 5e-6
        return

    def calc_ws_radius(self, struct):
        bohr2angst = 0.52917721
        vol_unit = struct.volume/struct.num_sites
        sws = (3*vol_unit/4.0/np.pi)**(1.0/3)/bohr2angst
        return sws

    def make_basis_array(self, struct):
        """Returns a 2D numpy array of the basis atom coordinates
           in !!Cartesian!! coordinates.
        """
        len_basis = struct.num_sites
        emto_basis = np.zeros((len_basis, 3))
        for i in range(len_basis):
            emto_basis[i, :] = struct.sites[i].coords
        return emto_basis

    def make_sites_array(self, struct):
        len_basis = struct.num_sites
        emto_sites = []
        for i in range(len_basis):
            emto_sites.append(struct.sites[i].specie.number)
        return emto_sites

    def make_cpa_sites_array(self, struct):
        len_basis = struct.num_sites
        self.atoms_cpa = []
        self.concs_cpa = []
        self.m_splits_cpa = []
        self.fxms_cpa = []
        self.s_wss_cpa = []
        self.ws_wsts_cpa = []
        self.qtrs_cpa = []
        self.nrms_cpa = []
        self.a_scrs_cpa = []
        self.b_scrs_cpa = []
        self.tetas_cpa = []
        self.phis_cpa = []
        self.asrs_cpa = []
        for i in range(len_basis):
            atom_number = struct.sites[i].specie.number
            for j in range(len(self.pmg_species)):
                if atom_number == self.pmg_species[j]:
                    self.atoms_cpa.append(self.species[j])
                    self.concs_cpa.append(self.concs[j])
                    self.m_splits_cpa.append(self.m_splits[j])
                    self.fxms_cpa.append(self.fxms[j])
                    self.s_wss_cpa.append(self.s_wss[j])
                    self.ws_wsts_cpa.append(self.ws_wsts[j])
                    self.qtrs_cpa.append(self.qtrs[j])
                    self.nrms_cpa.append(self.nrms[j])
                    self.a_scrs_cpa.append(self.a_scrs[j])
                    self.b_scrs_cpa.append(self.b_scrs[j])
                    self.tetas_cpa.append(self.tetas[j])
                    self.phis_cpa.append(self.phis[j])
                    self.asrs_cpa.append(self.asrs[j])
                    break

    def get_equivalent_sites(self):
        """Find all the sites that have exactly the same species,
        concentrations, and magnetic moments"""
        splt_tol = 1e-6
        conc_tol = 1e-6
        species_sorted = []
        m_splits_sorted = []
        concs_sorted = []
        for i in range(len(self.species)):
            tmp1 = []
            tmp2 = []
            tmp3 = []
            ind_sorted = np.argsort(self.species[i])
            for ind in ind_sorted:
                tmp1.append(self.species[i][ind])
                tmp2.append(self.m_splits[i][ind])
                tmp3.append(self.concs[i][ind])
            species_sorted.append(tmp1)
            m_splits_sorted.append(tmp2)
            concs_sorted.append(tmp3)

        eqv_sites = np.zeros((len(species_sorted), len(species_sorted)), dtype=np.int) + 9999
        for i in range(len(species_sorted)-1):
            for j in range(i+1, len(species_sorted)):
                eqv_sites[i,j] = 1
                if len(species_sorted[i]) != len(species_sorted[j]):
                    # Sites i and j contain different amound of atoms.
                    # For now, take them to be non-equivalent, although
                    # they could still be equivalent in the case that
                    # some element has been split into two or more parts
                    # concentration-wise (whole and the parts should have
                    # identical magnetic moments).
                    eqv_sites[i, j] = 0
                else:
                    for a1, a2, splt1, splt2, conc1, conc2 in zip(species_sorted[i], species_sorted[j],
                            m_splits_sorted[i], m_splits_sorted[j], concs_sorted[i], concs_sorted[j]):
                        if a1 != a2 or np.abs(splt1 - splt2) > splt_tol or np.abs(conc1 - conc2) > conc_tol:
                            # Some pair of atoms (in the sorted lists) were not
                            # the same => sites i and j are not equivalent.
                            eqv_sites[i, j] = 0
                            break
        output_sites = np.ones(len(species_sorted), dtype=np.int) * 9999
        next_available = 1
        for i in range(len(species_sorted)-1):
            if output_sites[i] == 9999:
                output_sites[i] = next_available
                next_available += 1
            for j in range(i+1, len(species_sorted)):
                if eqv_sites[i, j] == 1:
                    output_sites[j] = output_sites[i]
        if output_sites[-1] == 9999:
            output_sites[-1] = next_available
        return output_sites

    def get_crystal_system(self, n):
        """n = space group number"""
        def f(i, j):
            return i <= n <= j

        cs = {"triclinic": (1, 2), "monoclinic": (3, 15),
              "orthorhombic": (16, 74), "tetragonal": (75, 142),
              "trigonal": (143, 167), "hexagonal": (168, 194),
              "cubic": (195, 230)}

        crystal_system = None
        for k, v in cs.items():
            if f(*v):
                crystal_system = k
                break
        return crystal_system

    def get_space_group_number(self, sg):
        return int(sg.split()[-1].lstrip('(').rstrip(')'))

    def get_ibz_from_crystal_system(self, cs, lat):
        spg_cell = (
            lat,
            [[0.0, 0.0, 0.0]],
            [1]
        )
        space_group = spg.get_spacegroup(spg_cell)
        space_group_number = self.get_space_group_number(space_group)
        ibz = self.sg2ibz[space_group_number]

        ibz_out = None
        if cs == 'cubic':
            if ibz in (1,2,3):
                ibz_out = ibz
        elif cs == 'hexagonal':
            if ibz == 4:
                ibz_out = ibz
        elif cs == 'tetragonal':
            if ibz in (5,6):
                ibz_out = ibz
            elif ibz == 1:
                ibz_out = 5
            elif ibz == 3:
                ibz_out = 6
        elif  cs == 'trigonal':
            if ibz in (1,7):
                ibz_out = ibz
        elif cs == 'orthorhombic':
            if ibz in (8,9,10,11):
                ibz_out = ibz
            elif ibz in (1,5):
                # Base-centered orthorhombic without basis can
                # look like simple cubic or simple tetragonal.
                # Check the angles of the input lattice to
                # see whether the lattice vectors are those of
                # simple or base-centered.
                angles = np.array([get_angle(lat[0,:], lat[1,:]),
                                   get_angle(lat[0,:], lat[2,:]),
                                   get_angle(lat[1,:], lat[2,:])])
                if not all(np.abs(angles-90) < 1e-4):
                    ibz_out = 9
                else:
                    ibz_out = 8
            elif ibz == 2:
                ibz_out = 11
            elif ibz in (3,6):
                ibz_out = 10
        elif cs == 'monoclinic':
            if ibz in (12,13):
                ibz_out = ibz
            elif ibz in (1,5,8):
                num_angle_is_90 = 0
                angles = np.array([get_angle(lat[0,:], lat[1,:]),
                                   get_angle(lat[0,:], lat[2,:]),
                                   get_angle(lat[1,:], lat[2,:])])
                for angle in angles:
                    if np.abs(angle-90) < 1e-4:
                        num_angle_is_90 += 1
                if num_angle_is_90 < 2:
                    ibz_out = 13
                else:
                    ibz_out = 12
        elif cs == 'triclinic':
            if ibz == 14:
                ibz_out = 14
            elif ibz in (1,5,8,12):
                ibz_out = 14

        if ibz_out is None:
            raise ValueError(f"Cannot fit crystal system \'{cs}\'" +
                f" and IBZ of the lattice vectors \'{ibz}\' together!")
        else:
            return ibz_out

    def prepare_input_files(self, prims=None, basis=None, latpath=None,
                            coords_are_cartesian=False, latname=None,
                            species=None, find_primitive=True,
                            use_sublattices=True,
                            concs=None, m_splits=None, its=None, ws_wsts=None,
                            s_wss=None, make_supercell=None, fxms=None,
                            qtrs=None, nrms=None, a_scrs=None, b_scrs=None,
                            tetas=None, phis=None, asrs=None,
                            **kwargs):
        if prims is None:
            sys.exit('EMTO.init_structure(): \'prims\' has to be given!')
        if basis is None:
            sys.exit('EMTO.init_structure(): \'basis\' has to be given!')
        if latpath is None:
            self.latpath = os.getcwd()
        else:
            self.latpath = latpath
        if latname is None:
            self.latname = 'structure'
        else:
            self.latname = latname

        self.prims = np.array(prims)
        self.basis = np.array(basis)
        self.len_basis = len(self.basis[:, 0])

        if species is None:
            sys.exit('EMTO.init_structure(): \'species\' has to be given!')
        else:
            self.species = []
            for i in range(len(species)):
                if isinstance(species[i], list):
                    tmp = []
                    for j in range(len(species[i])):
                        tmp.append(species[i][j])
                    self.species.append(tmp)
                else:
                    self.species.append([species[i]])

        if m_splits is None:
            # Assume a zero moments array
            self.m_splits = []
            for i in range(len(self.species)):
                if isinstance(self.species[i], list):
                    tmp = []
                    for j in range(len(self.species[i])):
                        tmp.append(0.0)
                    self.m_splits.append(tmp)
                else:
                    self.m_splits.append([0.0])
        else:
            self.m_splits = []
            for i in range(len(m_splits)):
                if isinstance(m_splits[i], list):
                    tmp = []
                    for j in range(len(m_splits[i])):
                        tmp.append(m_splits[i][j])
                    self.m_splits.append(tmp)
                else:
                    self.m_splits.append([m_splits[i]])

        if fxms is None:
            # Assume a zero moments array
            self.fxms = []
            for i in range(len(self.species)):
                if isinstance(self.species[i], list):
                    tmp = []
                    for j in range(len(self.species[i])):
                        tmp.append('N')
                    self.fxms.append(tmp)
                else:
                    self.fxms.append(['N'])
        else:
            self.fxms = []
            for i in range(len(fxms)):
                if isinstance(fxms[i], list):
                    tmp = []
                    for j in range(len(fxms[i])):
                        tmp.append(fxms[i][j])
                    self.fxms.append(tmp)
                else:
                    self.fxms.append([fxms[i]])

        if nrms is None:
            # Assume a zero moments array
            self.nrms = []
            for i in range(len(self.species)):
                if isinstance(self.species[i], list):
                    tmp = []
                    for j in range(len(self.species[i])):
                        tmp.append(1)
                    self.nrms.append(tmp)
                else:
                    self.nrms.append([1])
        else:
            self.nrms = []
            for i in range(len(nrms)):
                if isinstance(nrms[i], list):
                    tmp = []
                    for j in range(len(nrms[i])):
                        tmp.append(nrms[i][j])
                    self.nrms.append(tmp)
                else:
                    self.nrms.append([nrms[i]])

        if a_scrs is None:
            # Assume a zero moments array
            self.a_scrs = []
            for i in range(len(self.species)):
                if isinstance(self.species[i], list):
                    tmp = []
                    for j in range(len(self.species[i])):
                        tmp.append(0.67)
                    self.a_scrs.append(tmp)
                else:
                    self.a_scrs.append([0.67])
        else:
            self.a_scrs = []
            for i in range(len(a_scrs)):
                if isinstance(a_scrs[i], list):
                    tmp = []
                    for j in range(len(a_scrs[i])):
                        tmp.append(a_scrs[i][j])
                    self.a_scrs.append(tmp)
                else:
                    self.a_scrs.append([a_scrs[i]])

        if b_scrs is None:
            # Assume a zero moments array
            self.b_scrs = []
            for i in range(len(self.species)):
                if isinstance(self.species[i], list):
                    tmp = []
                    for j in range(len(self.species[i])):
                        tmp.append(1.05)
                    self.b_scrs.append(tmp)
                else:
                    self.b_scrs.append([1.05])
        else:
            self.b_scrs = []
            for i in range(len(b_scrs)):
                if isinstance(b_scrs[i], list):
                    tmp = []
                    for j in range(len(b_scrs[i])):
                        tmp.append(b_scrs[i][j])
                    self.b_scrs.append(tmp)
                else:
                    self.b_scrs.append([b_scrs[i]])

        if tetas is None:
            # Assume a zero moments array
            self.tetas = []
            for i in range(len(self.species)):
                if isinstance(self.species[i], list):
                    tmp = []
                    for j in range(len(self.species[i])):
                        tmp.append(0.0)
                    self.tetas.append(tmp)
                else:
                    self.tetas.append([0.0])
        else:
            self.tetas = []
            for i in range(len(tetas)):
                if isinstance(tetas[i], list):
                    tmp = []
                    for j in range(len(tetas[i])):
                        tmp.append(tetas[i][j])
                    self.tetas.append(tmp)
                else:
                    self.tetas.append([tetas[i]])

        if phis is None:
            # Assume a zero moments array
            self.phis = []
            for i in range(len(self.species)):
                if isinstance(self.species[i], list):
                    tmp = []
                    for j in range(len(self.species[i])):
                        tmp.append(0.0)
                    self.phis.append(tmp)
                else:
                    self.phis.append([0.0])
        else:
            self.phis = []
            for i in range(len(phis)):
                if isinstance(phis[i], list):
                    tmp = []
                    for j in range(len(phis[i])):
                        tmp.append(phis[i][j])
                    self.phis.append(tmp)
                else:
                    self.phis.append([phis[i]])

        if asrs is None:
            # Assume a zero moments array
            self.asrs = []
            for i in range(len(self.species)):
                if isinstance(self.species[i], list):
                    tmp = []
                    for j in range(len(self.species[i])):
                        tmp.append(1.0)
                    self.asrs.append(tmp)
                else:
                    self.asrs.append([1.0])
        else:
            self.asrs = []
            for i in range(len(asrs)):
                if isinstance(asrs[i], list):
                    tmp = []
                    for j in range(len(asrs[i])):
                        tmp.append(asrs[i][j])
                    self.asrs.append(tmp)
                else:
                    self.asrs.append([asrs[i]])

        # EMTOx does not use s_wss, ws_wsts, qtrs
        if s_wss is None:
            self.s_wss = []
            for i in range(len(self.species)):
                if isinstance(self.species[i], list):
                    tmp = []
                    for j in range(len(self.species[i])):
                        tmp.append(1.0)
                    self.s_wss.append(tmp)
                else:
                    self.s_wss.append([1.0])
        else:
            self.s_wss = []
            for i in range(len(s_wss)):
                if isinstance(s_wss[i], list):
                    tmp = []
                    for j in range(len(s_wss[i])):
                        tmp.append(s_wss[i][j])
                    self.s_wss.append(tmp)
                else:
                    tmp = []
                    for j in range(len(self.species[i])):
                        tmp.append(s_wss[i])
                    self.s_wss.append(tmp)

        if ws_wsts is None:
            self.ws_wsts = []
            for i in range(len(self.species)):
                if isinstance(self.species[i], list):
                    tmp = []
                    for j in range(len(self.species[i])):
                        tmp.append(1.0)
                    self.ws_wsts.append(tmp)
                else:
                    self.ws_wsts.append([1.0])
        else:
            self.ws_wsts = []
            for i in range(len(ws_wsts)):
                if isinstance(ws_wsts[i], list):
                    tmp = []
                    for j in range(len(ws_wsts[i])):
                        tmp.append(ws_wsts[i][j])
                    self.ws_wsts.append(tmp)
                else:
                    tmp = []
                    for j in range(len(self.species[i])):
                        tmp.append(ws_wsts[i])
                    self.ws_wsts.append(tmp)

        if qtrs is None:
            self.qtrs = []
            for i in range(len(self.species)):
                if isinstance(self.species[i], list):
                    tmp = []
                    for j in range(len(self.species[i])):
                        tmp.append(0.0)
                    self.qtrs.append(tmp)
                else:
                    self.qtrs.append([0.0])
        else:
            self.qtrs = []
            for i in range(len(qtrs)):
                if isinstance(qtrs[i], list):
                    tmp = []
                    for j in range(len(qtrs[i])):
                        tmp.append(qtrs[i][j])
                    self.qtrs.append(tmp)
                else:
                    tmp = []
                    for j in range(len(self.species[i])):
                        tmp.append(qtrs[i])
                    self.qtrs.append(tmp)

        if concs is None:
            # Assume a zero moments array
            self.concs = []
            for i in range(len(self.species)):
                if isinstance(self.species[i], list):
                    tmp = []
                    for j in range(len(self.species[i])):
                        tmp.append(1.0/len(self.species[i]))
                    self.concs.append(tmp)
                else:
                    self.concs.append([1.0])
        else:
            self.concs = []
            for i in range(len(concs)):
                if isinstance(concs[i], list):
                    tmp = []
                    tmp_sum = 0.0
                    for j in range(len(concs[i])):
                        tmp.append(concs[i][j])
                        tmp_sum += concs[i][j]
                    if tmp_sum < 1.1:
                        if np.abs(tmp_sum - 1.0) > 1.e-6:
                            sys.exit('Concentrations {0} for site {1} do not add up to 1.0!!!'.format(concs[i], i+1))
                    else:
                        if np.abs(tmp_sum - 100.0) > 1.e-3:
                            sys.exit('Concentrations {0} for site {1} do not add up to 100!!!'.format(concs[i], i+1))                        
                    self.concs.append(tmp)
                else:
                    self.concs.append([concs[i]])

        # Check that all species, concs, and m_splits arrays have the same dimensions
        for a, b in combinations([self.basis, self.species, self.concs,
            self.m_splits, self.fxms, self.nrms, self.a_scrs, self.b_scrs,
            self.tetas, self.phis, self.asrs], 2):
            if len(a) != len(b):
                print(a, 'len = ', len(a))
                print(b, 'len = ', len(b))
                sys.exit('The above input arrays have inconsistent lengths!!!')
        for a, b in combinations([self.species, self.concs, self.m_splits,
            self.fxms, self.nrms, self.a_scrs, self.b_scrs, self.tetas,
            self.phis, self.asrs], 2):
            for sublist1, sublist2 in zip(a, b):
                if len(sublist1) != len(sublist2):
                    print(sublist1, 'len = ', len(sublist1))
                    print(sublist2, 'len = ', len(sublist2))
                    sys.exit('The above input array elements have inconsistent lengths!!!')

        self.find_primitive = find_primitive
        self.pmg_species = self.get_equivalent_sites()
        #
        self.coords_are_cartesian = coords_are_cartesian
        self.ibz = None
        self.make_supercell = make_supercell
        #
        self.pmg_input_lattice = Lattice(self.prims)
        self.pmg_input_struct = Structure(self.pmg_input_lattice, self.pmg_species, self.basis,
                                          coords_are_cartesian=self.coords_are_cartesian)
        #
        # if self.make_supercell is not None:
            # self.pmg_input_struct.make_supercell(self.make_supercell)
        #
        self.sws = self.calc_ws_radius(self.pmg_input_struct)
        #
        self.finder = SpacegroupAnalyzer(self.pmg_input_struct, symprec=0.0001,
            angle_tolerance=0.0001)
        self.stm = StructureMatcher(ltol=0.001, stol=0.001, angle_tol=0.001,
            attempt_supercell=True)
        #
        print("Input structure information:")
        print(self.pmg_input_struct)
        print("Volume: ", self.pmg_input_struct.volume)
        print("Lattice vectors:")
        print(self.pmg_input_struct.lattice.matrix)
        print("")
        #
        # spglib
        spg_cell = (
            self.pmg_input_lattice.matrix,
            self.pmg_input_struct.frac_coords,
            self.pmg_species
        )
    
        if not self.find_primitive:
            space_group = spg.get_spacegroup(spg_cell)
            space_group_number = self.get_space_group_number(space_group)
            crystal_system = self.get_crystal_system(space_group_number)
            self.spg_ibz = self.get_ibz_from_crystal_system(crystal_system,
                self.pmg_input_lattice.matrix)
            self.ibz = self.spg_ibz

            self.equivalent_atoms = spg.get_symmetry(spg_cell)['equivalent_atoms']

            # Leaving only unique symmetrically equivalent atoms in the
            # structure does not work.
            # unique_atoms = []
            # include_atoms = []
            # for i in range(len(eq_atoms)):
                # if eq_atoms[i] not in unique_atoms:
                    # include_atoms.append(i)
                    # unique_atoms.append(eq_atoms[i])
            # pmg_eq_unique_species = []
            # pmg_eq_unique_basis = []
            # for ind in include_atoms:
                # pmg_eq_unique_species.append(self.pmg_species[ind])
                # pmg_eq_unique_basis.append(self.pmg_input_struct.frac_coords[ind])
            # spg_cell = (
                # self.pmg_input_lattice.matrix,
                # pmg_eq_unique_basis,
                # pmg_eq_unique_species
                # )

            self.spg_space_group = spg.get_spacegroup(spg_cell)
            self.spg_space_group_number = self.get_space_group_number(self.spg_space_group)
            self.spg_space_group_symbol = self.spg_space_group
            self.spg_prim_lat = self.pmg_input_lattice.matrix
            self.spg_prim_pos = self.pmg_input_struct.frac_coords
            self.spg_prim_species = self.pmg_species
        else:
            self.spg_prim_lat, self.spg_prim_pos, self.spg_prim_species = \
            spg.standardize_cell(spg_cell, to_primitive=True)
            self.spg_space_group = spg.get_spacegroup((
                self.spg_prim_lat, self.spg_prim_pos, self.spg_prim_species
                ))
            tmp = spg.get_symmetry((
                self.spg_prim_lat, self.spg_prim_pos, self.spg_prim_species
                ))
            self.equivalent_atoms = tmp['equivalent_atoms']

            self.spg_space_group_number = self.get_space_group_number(self.spg_space_group)
            self.spg_space_group_symbol = self.spg_space_group
            self.spg_ibz = self.sg2ibz[self.spg_space_group_number]
            self.ibz = self.spg_ibz

        self.prim_struct = Structure(Lattice(self.spg_prim_lat),
            self.spg_prim_species, self.spg_prim_pos)

        print("Output structure information:")
        print(self.prim_struct)
        print("Volume: ", self.prim_struct.volume)
        print("Lattice vectors:")
        print(self.prim_struct.lattice.matrix)
        print("")
        #
        self.primaa = self.prim_struct.lattice.matrix[0, :]
        self.primbb = self.prim_struct.lattice.matrix[1, :]
        self.primcc = self.prim_struct.lattice.matrix[2, :]
        self.output_basis = self.make_basis_array(self.prim_struct)
        # Below we calculate the transformation that maps
        # self.primaX to lattice vectors used by EMTO.
        # This transform depends on the type of the Bravais lattice,
        # so each case must be treated separately.
        if self.spg_ibz == 1:
            norm_tmp = np.linalg.norm(self.primaa)
            self.output_prima = self.primaa/norm_tmp
            self.output_primb = self.primbb/norm_tmp
            self.output_primc = self.primcc/norm_tmp
            # Apply transformation on the basis atoms
            self.output_basis = self.output_basis/norm_tmp
            self.output_boa = 0.0
            self.output_coa = 0.0
            self.output_alpha = 0.0
            self.output_beta = 0.0
            self.output_gamma = 0.0
            self.emto_prima = np.array([1, 0, 0])
            self.emto_primb = np.array([0, 1, 0])
            self.emto_primc = np.array([0, 0, 1])
            self.emto_basis = self.output_basis

        elif self.spg_ibz == 2:
            norm_tmp = 2*self.primaa[1]
            self.output_prima = self.primcc/norm_tmp
            self.output_primb = self.primaa/norm_tmp
            self.output_primc = self.primbb/norm_tmp
            # Apply transformation on the basis atoms
            self.output_basis = self.output_basis/norm_tmp
            self.output_boa = 0.0
            self.output_coa = 0.0
            self.output_alpha = 0.0
            self.output_beta = 0.0
            self.output_gamma = 0.0
            self.emto_prima = np.array([0.5, 0.5, 0])
            self.emto_primb = np.array([0, 0.5, 0.5])
            self.emto_primc = np.array([0.5, 0, 0.5])
            self.emto_basis = self.output_basis

        elif self.spg_ibz == 3:
            norm_tmp = 2*self.primaa[1]
            self.output_prima = self.primcc/norm_tmp
            self.output_primb = self.primaa/norm_tmp
            self.output_primc = self.primbb/norm_tmp
            # Apply transformation on the basis atoms
            self.output_basis = self.output_basis/norm_tmp
            self.output_boa = 0.0
            self.output_coa = 0.0
            self.output_alpha = 0.0
            self.output_beta = 0.0
            self.output_gamma = 0.0
            self.emto_prima = np.array([0.5, 0.5, -0.5])
            self.emto_primb = np.array([-0.5, 0.5, 0.5])
            self.emto_primc = np.array([0.5, -0.5, 0.5])
            self.emto_basis = self.output_basis

        elif self.spg_ibz == 4:
            rot1 = rotation_matrix([0.0, 0.0, 1.0], 0./180*np.pi)
            self.output_prima = np.dot(rot1, self.primaa)
            self.output_primb = np.dot(rot1, self.primbb)
            self.output_primc = np.dot(rot1, self.primcc)
            for i in range(len(self.output_basis[:, 0])):
                self.output_basis[i, :] = np.dot(rot1, self.output_basis[i, :])
            self.output_boa = 0.0
            self.output_coa = self.output_primc[2]
            self.output_alpha = 0.0
            self.output_beta = 0.0
            self.output_gamma = 0.0
            # EMTO convention:
            self.emto_prima = np.array([1., 0, 0])
            self.emto_primb = np.array([-0.5, np.sqrt(3.)/2, 0])
            self.emto_primc = np.array([0., 0, self.output_coa])
            self.emto_basis = self.output_basis

        elif self.spg_ibz == 5:
            norm_tmp = self.primaa[0]
            self.output_prima = self.primaa/norm_tmp
            self.output_primb = self.primbb/norm_tmp
            self.output_primc = self.primcc/norm_tmp
            # Apply transformation on the basis atoms
            self.output_basis = self.output_basis/norm_tmp
            self.output_boa = 0.0
            self.output_coa = self.output_primc[2]
            self.output_alpha = 0.0
            self.output_beta = 0.0
            self.output_gamma = 0.0
            self.emto_prima = np.array([1.0, 0.0, 0.0])
            self.emto_primb = np.array([0.0, 1.0, 0.0])
            self.emto_primc = np.array([0.0, 0.0, self.output_coa])
            self.emto_basis = self.output_basis

        elif self.spg_ibz == 6:
            self.output_prima = self.primbb
            self.output_primb = self.primcc
            self.output_primc = self.primaa
            # Apply transformation on the basis atoms
            self.output_basis = self.output_basis
            self.output_boa = 0.0
            self.output_coa = 2*self.output_prima[2]
            self.output_alpha = 0.0
            self.output_beta = 0.0
            self.output_gamma = 0.0
            self.emto_prima = np.array([0.5, -0.5, self.output_coa/2])
            self.emto_primb = np.array([0.5, 0.5, -self.output_coa/2])
            self.emto_primc = np.array([-0.5, 0.5, self.output_coa/2])
            self.emto_basis = self.output_basis

        elif self.spg_ibz == 7:
            alpha = self.prim_struct.lattice.alpha
            kulma = np.arctan((self.primaa[0]+self.primbb[0]+self.primcc[0])/
                              (self.primaa[2]+self.primbb[2]+self.primcc[2]))
            rot1 = rotation_matrix([0.0, -1.0, 0.0], kulma)
            rot2 = np.array([[-np.sqrt(3.0)/2, -0.5, 0.0],
                             [0.5, -np.sqrt(3.0)/2, 0.0],
                             [0.0, 0.0, 1.0]])
            self.output_prima = np.dot(rot2, np.dot(rot1, self.primaa))
            self.output_primb = np.dot(rot2, np.dot(rot1, self.primbb))
            self.output_primc = np.dot(rot2, np.dot(rot1, self.primcc))
            scale_a = self.output_prima[1]
            print('scale_a = ',scale_a)
            self.output_prima = self.output_prima/scale_a
            self.output_primb = self.output_primb/scale_a
            self.output_primc = self.output_primc/scale_a
            # Apply transformation on the basis atoms
            for i in range(len(self.output_basis[:, 0])):
                self.output_basis[i,:] = np.dot(rot2, np.dot(rot1, self.output_basis[i, :]))/scale_a
            self.output_boa = 1.0
            self.output_coa = self.output_prima[2]
            self.output_alpha = 0.0
            self.output_beta = 0.0
            self.output_gamma = 0.0
            self.emto_prima = np.array([0.0, 1.0, self.output_coa])
            self.emto_primb = np.array([-np.sqrt(3.)/2, -0.5, self.output_coa])
            self.emto_primc = np.array([np.sqrt(3.)/2, -0.5, self.output_coa])
            self.emto_basis = self.output_basis

        elif self.spg_ibz == 8:
            if (np.abs(self.primaa[0]) < np.abs(self.primbb[1])) and \
               (np.abs(self.primbb[1]) < np.abs(self.primcc[2])):
                norm_tmp = self.primaa[0]
                self.output_prima = self.primaa/norm_tmp
                self.output_primb = self.primbb/norm_tmp
                self.output_primc = self.primcc/norm_tmp
                # Apply transformation on the basis atoms
                self.output_basis = self.output_basis/norm_tmp

            elif np.abs(np.abs(self.primaa[0]) - np.abs(self.primbb[1])) < 1.e-6 and \
               np.abs(self.primbb[1]) < np.abs(self.primcc[2]):
                norm_tmp = self.primaa[0]
                self.output_prima = self.primaa/norm_tmp
                self.output_primb = self.primbb/norm_tmp
                self.output_primc = self.primcc/norm_tmp
                # Apply transformation on the basis atoms
                self.output_basis = self.output_basis/norm_tmp
                
            elif np.abs(self.primaa[0]) < np.abs(self.primcc[2]):
                norm_tmp = self.primcc[2]
                rot1 = rotation_matrix([0.0, 0.0, 1.0], -90./180*np.pi)
                rot2 = rotation_matrix([-1.0, 0.0, 0.0], 90./180*np.pi)
                self.output_prima = np.dot(rot2, np.dot(rot1, self.primbb))/norm_tmp
                self.output_primb = np.dot(rot2, np.dot(rot1, self.primcc))/norm_tmp
                self.output_primc = np.dot(rot2, np.dot(rot1, self.primaa))/norm_tmp
                print(self.output_prima)
                print(self.output_primb)
                print(self.output_primc)
                # Apply transformation on the basis atoms
                for i in range(len(self.output_basis[:, 0])):
                    self.output_basis[i,:] = np.dot(rot2, np.dot(rot1, self.output_basis[i, :]))/norm_tmp
            else:
                norm_tmp = self.primaa[0]
                self.output_prima = self.primaa/norm_tmp
                self.output_primb = self.primbb/norm_tmp
                self.output_primc = self.primcc/norm_tmp
                # Apply transformation on the basis atoms
                self.output_basis = self.output_basis/norm_tmp
            #
            self.output_boa = self.output_primb[1]
            self.output_coa = self.output_primc[2]
            self.output_alpha = 0.0
            self.output_beta = 0.0
            self.output_gamma = 0.0
            self.emto_prima = np.array([1.0, 0.0, 0.0])
            self.emto_primb = np.array([0.0, self.output_boa, 0.0])
            self.emto_primc = np.array([0.0, 0.0 ,self.output_coa])
            self.emto_basis = self.output_basis

        elif self.spg_ibz == 9:
            if np.abs(self.primbb[1] - 0.5) < 1e-12 and \
               np.abs(self.primcc[1] + 0.5) < 1e-12:
                rot1 = rotation_matrix([0.0, 1.0, 0.0], 90./180*np.pi)
                self.output_prima = np.dot(rot1, self.primaa)
                self.output_primb = np.dot(rot1, self.primbb)
                self.output_primc = np.dot(rot1, self.primcc)
                # Redefine lattice vectors
                tmp = np.copy(self.output_prima)
                self.output_prima[:] = self.output_primc[:]
                self.output_primc[:] = tmp
                # Mirror along the xy-plane
                self.output_primc *= -1
                # Scale lattice vectors so that a1 and a2 x-components are 0.5
                norm_tmp = 2*self.output_prima[0]
                self.output_prima /= norm_tmp
                self.output_primb /= norm_tmp
                self.output_primc /= norm_tmp
                # Apply transformation on the basis atoms
                for i in range(len(self.output_basis[:, 0])):
                    self.output_basis[i,:] = np.dot(rot1, self.output_basis[i, :])
                for i in range(len(self.output_basis[:, 0])):
                    self.output_basis[i,2] *= -1
                self.output_basis /= norm_tmp
                #print(self.output_prima)
                #print(self.output_primb)
                #print(self.output_primc)
                
            else:
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
            # EMTO convention:
            self.emto_prima = np.array([0.5, -self.output_boa/2, 0])
            self.emto_primb = np.array([0.5, self.output_boa/2, 0])
            self.emto_primc = np.array([0, 0, self.output_coa])
            self.emto_basis = self.output_basis

        elif self.spg_ibz == 10:
            self.output_prima = np.zeros_like(self.primaa)
            self.output_primb = np.zeros_like(self.primbb)
            self.output_primc = np.zeros_like(self.primcc)
            self.output_prima[0] = self.primaa[1]
            self.output_prima[1] = self.primaa[0]
            self.output_prima[2] = self.primaa[2]
            self.output_primb[0] = self.primcc[1]
            self.output_primb[1] = self.primcc[0]
            self.output_primb[2] = self.primcc[2]
            self.output_primc[0] = self.primbb[1]
            self.output_primc[1] = self.primbb[0]
            self.output_primc[2] = self.primbb[2]
            norm_tmp = 2*self.output_prima[0]
            self.output_prima /= norm_tmp
            self.output_primb /= norm_tmp
            self.output_primc /= norm_tmp
            # Apply transformation on the basis atoms
            basis_tmp = np.copy(self.output_basis)
            self.output_basis[:, 0] = basis_tmp[:, 1]
            self.output_basis[:, 1] = basis_tmp[:, 0]
            self.output_basis = self.output_basis/norm_tmp
            self.output_boa = 2*self.output_primc[1]
            self.output_coa = 2*self.output_primc[2]
            self.output_alpha = 0.0
            self.output_beta = 0.0
            self.output_gamma = 0.0
            self.emto_prima = np.array([0.5, -self.output_boa/2, self.output_coa/2])
            self.emto_primb = np.array([0.5, self.output_boa/2, -self.output_coa/2])
            self.emto_primc = np.array([-0.5, self.output_boa/2, self.output_coa/2])
            self.emto_basis = self.output_basis

        elif self.spg_ibz == 11:
            rot1 = rotation_matrix([1, 1, 1], 120./180*np.pi)
            self.output_prima = np.dot(rot1, self.primaa)
            self.output_primb = np.dot(rot1, self.primbb)
            self.output_primc = np.dot(rot1, self.primcc)
            norm_tmp = 2*self.output_prima[0]
            self.output_prima /= norm_tmp
            self.output_primb /= norm_tmp
            self.output_primc /= norm_tmp
            # Apply transformation on the basis atoms
            for i in range(len(self.output_basis[:, 0])):
                self.output_basis[i, :] = np.dot(rot1, self.output_basis[i, :])
            self.output_basis /= norm_tmp
            self.output_boa = 2*self.output_primc[1]
            self.output_coa = 2*self.output_primc[2]
            self.output_alpha = 0.0
            self.output_beta = 0.0
            self.output_gamma = 0.0
            # EMTO convention:
            self.emto_prima = np.array([0.5, 0, self.output_coa/2])
            self.emto_primb = np.array([0.5, self.output_boa/2, 0])
            self.emto_primc = np.array([0, self.output_boa/2, self.output_coa/2])
            self.emto_basis = self.output_basis

        elif self.spg_ibz == 12:
            bc_norm = np.linalg.norm(self.primaa)
            # Rotate 90 degreen counter clockwise around the x-axis
            rot1 = rotation_matrix([1, 0, 0], -90./180*np.pi)
            self.output_prima = np.dot(rot1, self.primaa/bc_norm)
            self.output_primb = np.dot(rot1, self.primcc/bc_norm)
            self.output_primc = np.dot(rot1, self.primbb/bc_norm)
            # Mirror a3 from negative z-axis to positive side
            self.output_primc *= -1.0
            # spg uses gamma > 90, so we redefine the a3 lattice vector so that
            # gamma < 90:
            self.output_primb[0] *= -1.0
            gamma = get_angle(self.output_prima, self.output_primb)
            y_fac = self.output_primb[1]
            shift = np.abs(2*self.output_primb[0])
            #
            # Apply transformation on the basis atoms
            for i in range(len(self.output_basis[:, 0])):
                self.output_basis[i, :] = np.dot(rot1, self.output_basis[i, :])/bc_norm
            # Transform basis because self.output_primc was mirrored:
            for i in range(len(self.output_basis[:, 0])):
                self.output_basis[i, 2] *= -1.0
            # Transform basis because gamma was changed above:
            for i in range(len(self.output_basis[:, 0])):
                #self.output_basis[i, :] = np.dot(shift_mat, self.output_basis[i, :])
                if self.output_basis[i, 1] > 0:
                    self.output_basis[i, 0] += shift * np.abs(self.output_basis[i, 1] / y_fac)
                else:
                    self.output_basis[i, 0] -= shift * np.abs(self.output_basis[i, 1] / y_fac)
            self.output_boa = np.linalg.norm(self.output_primb)
            self.output_coa = np.linalg.norm(self.output_primc)
            self.output_alpha = 0.0
            self.output_beta = 0.0
            self.output_gamma = gamma
            self.emto_prima = np.array([1.0, 0, 0])
            self.emto_primb = np.array([self.output_boa*np.cos(np.radians(self.output_gamma)),
                                        self.output_boa*np.sin(np.radians(self.output_gamma)), 0])
            self.emto_primc = np.array([0, 0, self.output_coa])
            self.emto_basis = self.output_basis

        elif self.spg_ibz == 13:
            gamma = get_angle(self.primcc, self.primaa+self.primbb)
            switch_x_y = np.array([[0, -1, 0],
                                       [1,  0, 0],
                                       [0,  0, 1]])
            rot1 = np.array([[1.0,0.0,0.0],
                             [0.0,np.cos(np.radians(180-gamma)),-np.sin(np.radians(180-gamma))],
                             [0.0,np.sin(np.radians(180-gamma)),np.cos(np.radians(180-gamma))]])
            rot2 = np.array([[0.0,0.0,1.0],
                             [0.0,1.0,0.0],
                             [-1.0,0.0,0.0]])
            bc_norm = np.linalg.norm(self.primaa+self.primbb)
            self.output_prima = np.dot(rot2, np.dot(rot1, np.dot(switch_x_y, self.primcc)))/bc_norm
            self.output_primb = np.dot(rot2, np.dot(rot1, np.dot(switch_x_y, self.primaa)))/bc_norm
            self.output_primc = np.dot(rot2, np.dot(rot1, np.dot(switch_x_y, self.primbb)))/bc_norm
            # Apply transformation on the basis atoms
            for i in range(len(self.output_basis[:, 0])):
                self.output_basis[i, :] = np.dot(rot2, np.dot(rot1, np.dot(switch_x_y, self.output_basis[i, :])))/bc_norm
            self.output_boa = np.abs(self.output_prima[1])
            self.output_coa = np.abs(2*self.output_primc[2])
            self.output_alpha = 0.0
            self.output_beta = 0.0
            self.output_gamma = gamma
            self.emto_prima = np.array([0.0, -self.output_boa, 0])
            self.emto_primb = np.array([0.5*np.sin(np.radians(self.output_gamma)),
                                        -0.5*np.cos(np.radians(self.output_gamma)),
                                        -self.output_coa/2])
            self.emto_primc = np.array([0.5*np.sin(np.radians(self.output_gamma)),
                                        -0.5*np.cos(np.radians(self.output_gamma)),
                                        self.output_coa/2])
            self.emto_basis = self.output_basis

        elif self.spg_ibz == 14:
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
            self.emto_prima = np.array([1.0, 0, 0])
            self.emto_primb = np.array([self.output_boa*np.cos(np.radians(self.output_gamma)),
                                        self.output_boa*np.sin(np.radians(self.output_gamma)),
                                        0])
            self.emto_primc = np.array([self.output_coa*np.cos(np.radians(self.output_beta)),
                                        self.output_coa*(np.cos(np.radians(self.output_alpha)) -
                                        np.cos(np.radians(self.output_beta)) *
                                        np.cos(np.radians(self.output_gamma))) / np.sin(np.radians(self.output_gamma)),
                                        self.output_coa*np.sqrt(1 - np.cos(np.radians(self.output_gamma))**2 -
                                        np.cos(np.radians(self.output_alpha))**2 -
                                        np.cos(np.radians(self.output_beta))**2 +
                                        2*np.cos(np.radians(self.output_alpha))*
                                        np.cos(np.radians(self.output_beta))*
                                        np.cos(np.radians(self.output_gamma)))/np.sin(np.radians(self.output_gamma))])
            self.emto_basis = self.output_basis

        self.output_sites = self.make_sites_array(self.prim_struct)
        self.output_lattice = Lattice(np.array([self.emto_prima, self.emto_primb, self.emto_primc]))
        self.output_struct = Structure(self.output_lattice, self.output_sites,
                                       self.emto_basis, coords_are_cartesian=True)
        #
        # Print EMTO structure information
        print("")
        print("Generated EMTO structure:")
        print(self.output_struct)
        print("Volume: ", self.output_struct.volume)
        print("WS-rad: ", self.sws)
        print("Lattice vectors:")
        print(self.output_struct.lattice.matrix)
        print("Basis vectors:")
        for i in range(len(self.output_struct.sites)):
            print(self.output_struct.sites[i].coords)
        print("")
        # Print symmetry information
        print("spglib reports the following information:")
        print("The spacegroup symbol of input structure: {}".format(self.spg_space_group))
        print("The spacegroup number of input structure: {}".format(self.spg_space_group_number))
        print("The Bravais lattice of input structure  : {}".format(self.sg2bl[self.spg_space_group_number]))
        print("Number of basis atoms                   : {}".format(self.prim_struct.num_sites))
        print("EMTO IBZ                                : {}".format(self.spg_ibz))
        print("")
        emto_cell = (
        self.output_lattice.matrix,
        self.output_struct.frac_coords,
        self.output_sites
        )
        self.emto_space_group = spg.get_spacegroup(emto_cell)
        self.emto_space_group_number = int(self.emto_space_group.split()[-1].lstrip('(').rstrip(')'))
        self.emto_space_group_symbol = self.emto_space_group
        self.emto_prim_lat, self.emto_prim_pos, self.emto_prim_species = spg.standardize_cell(emto_cell, to_primitive=True)
        self.emto_struct = Structure(Lattice(self.emto_prim_lat), self.emto_prim_species, self.emto_prim_pos)
        self.emto_ibz = self.sg2ibz[self.emto_space_group_number]
        print("spglib reports the following information:")
        print("The spacegroup symbol of EMTO structure : {}".format(self.emto_space_group))
        print("The spacegroup number of EMTO structure : {}".format(self.emto_space_group_number))
        print("The Bravais lattice of EMTO structure   : {}".format(self.sg2bl[self.emto_space_group_number]))
        print("Number of basis atoms                   : {}".format(self.output_struct.num_sites))
        print("EMTO IBZ                                : {}".format(self.emto_ibz))
        print("")

        # print(self.prim_struct)
        # print(self.emto_struct)
        
        # should_exit = False
        # if (self.spg_space_group != self.emto_space_group):
            # print("Input and output spacegroups are different!!!")
            # should_exit = True
        # if (self.spg_ibz != self.emto_ibz):
            # print("Input and output IBZ are different!!!")
            # should_exit = True
        # if should_exit:
            # sys.exit("Structure conversion went wrong! Check the symmetry information above.")
        #
        fitted_angles = [get_angle(self.output_prima, self.emto_prima),
            get_angle(self.output_primb, self.emto_primb),
            get_angle(self.output_primc, self.emto_primc)]
        for i, angle in enumerate(fitted_angles):
            #print(angle)
            if angle > self.fit_angle_tol:
                sys.exit('Error: Angle between lattice vectors {0} is {1} > {2}!!!'.format(i+1, angle, self.fit_angle_tol))
        fitted_ratios = [np.linalg.norm(self.output_prima) / np.linalg.norm(self.emto_prima),
            np.linalg.norm(self.output_primb) / np.linalg.norm(self.emto_primb),
            np.linalg.norm(self.output_primc) / np.linalg.norm(self.emto_primc)]
        for i, ratio in enumerate(fitted_ratios):
            #print(ratio)
            if np.abs(ratio - 1.0) > self.fit_norm_ratio_tol:
                sys.exit('Error: Ratio between lattice vector {0} norms is {1} > {2}!!!'.format(i+1, ratio, self.fit_norm_ratio_tol))
        print('Structure similarity check (input vs. output for EMTO):')
        fit1 = self.stm.fit_anonymous(self.pmg_input_struct, self.prim_struct)
        fit2 = self.stm.fit(self.pmg_input_struct, self.prim_struct)
        fit3 = self.stm.fit_anonymous(self.prim_struct, self.output_struct)
        fit4 = self.stm.fit(self.prim_struct, self.output_struct)
        fit5 = self.stm.fit_anonymous(self.pmg_input_struct, self.output_struct)
        fit6 = self.stm.fit(self.pmg_input_struct, self.output_struct)
        print('Input  -> spglib (sites only)     ?: ', fit1)
        print('Input  -> spglib (sites+chemistry)?: ', fit2)
        print('spglib -> EMTO   (sites only)     ?: ', fit3)
        print('spglib -> EMTO   (sites+chemistry)?: ', fit4)
        print('Input  -> EMTO   (sites only)     ?: ', fit5)
        print('Input  -> EMTO   (sites+chemistry)?: ', fit6)
        print("")
        
        if not all([fit1, fit2, fit3, fit4, fit5, fit6]):
            sys.exit('Some structures are not identical (check for False above) !!!')


        # Finally, save atoms, m_splits, and concs of the output structure
        # to be read by init_bulk function.
        # EMTOx needs ASR in the shape.
        self.make_cpa_sites_array(self.output_struct)
        #
        # Prepare KGRN, KFCD, and SLURM input files next.
        if self.ibz is None:
            sys.exit('self.ibz == None! Run create_structure_input() to generate IBZ \n'+
                     'for your structure.')
        # Construct an index array to keep track of the number of atoms in each site.
        if self.atoms_cpa is None:
            sys.exit('EMTO.init_bulk(): \'self.atoms_cpa\' does not exist!!! (Did you run init_structure?)')
        else:
            index_array = np.ones(len(self.atoms_cpa), dtype='int32')
            for i in range(len(self.atoms_cpa)):
                if isinstance(self.atoms_cpa[i], list):
                    index_array[i] = len(self.atoms_cpa[i])
                else:
                    index_array[i] = 1
            index_len = np.sum(index_array)
            atoms_flat = []
            for i in range(len(self.atoms_cpa)):
                if isinstance(self.atoms_cpa[i], list):
                    for j in range(len(self.atoms_cpa[i])):
                        atoms_flat.append(self.atoms_cpa[i][j])
                else:
                    atoms_flat.append(self.atoms_cpa[i])
            self.KGRN_atoms = np.array(atoms_flat)

        if self.m_splits_cpa is None:
            sys.exit('EMTO.init_bulk(): \'self.m_splits_cpa\' does not exist!!! (Did you run init_structure?)')
        else:
            m_splits_flat = []
            for i in range(len(self.m_splits_cpa)):
                if isinstance(self.m_splits_cpa[i], list):
                    for j in range(len(self.m_splits_cpa[i])):
                        m_splits_flat.append(self.m_splits_cpa[i][j])
                else:
                    m_splits_flat.append(self.m_splits_cpa[i])
            self.KGRN_m_splits = np.array(m_splits_flat)

        if self.fxms_cpa is None:
            sys.exit('EMTO.init_bulk(): \'self.fxms_cpa\' does not exist!!! (Did you run init_structure?)')
        else:
            fxms_flat = []
            for i in range(len(self.fxms_cpa)):
                if isinstance(self.fxms_cpa[i], list):
                    for j in range(len(self.fxms_cpa[i])):
                        fxms_flat.append(self.fxms_cpa[i][j])
                else:
                    fxms_flat.append(self.fxms_cpa[i])
            self.KGRN_fxms = np.array(fxms_flat)

        if self.nrms_cpa is None:
            sys.exit('EMTO.init_bulk(): \'self.nrms_cpa\' does not exist!!! (Did you run init_structure?)')
        else:
            nrms_flat = []
            for i in range(len(self.nrms_cpa)):
                if isinstance(self.nrms_cpa[i], list):
                    for j in range(len(self.nrms_cpa[i])):
                        nrms_flat.append(self.nrms_cpa[i][j])
                else:
                    nrms_flat.append(self.nrms_cpa[i])
            self.KGRN_nrms = np.array(nrms_flat)

        if self.a_scrs_cpa is None:
            sys.exit('EMTO.init_bulk(): \'self.a_scrs_cpa\' does not exist!!! (Did you run init_structure?)')
        else:
            a_scrs_flat = []
            for i in range(len(self.a_scrs_cpa)):
                if isinstance(self.a_scrs_cpa[i], list):
                    for j in range(len(self.a_scrs_cpa[i])):
                        a_scrs_flat.append(self.a_scrs_cpa[i][j])
                else:
                    a_scrs_flat.append(self.a_scrs_cpa[i])
            self.KGRN_a_scrs = np.array(a_scrs_flat)

        if self.b_scrs_cpa is None:
            sys.exit('EMTO.init_bulk(): \'self.b_scrs_cpa\' does not exist!!! (Did you run init_structure?)')
        else:
            b_scrs_flat = []
            for i in range(len(self.b_scrs_cpa)):
                if isinstance(self.b_scrs_cpa[i], list):
                    for j in range(len(self.b_scrs_cpa[i])):
                        b_scrs_flat.append(self.b_scrs_cpa[i][j])
                else:
                    b_scrs_flat.append(self.b_scrs_cpa[i])
            self.KGRN_b_scrs = np.array(b_scrs_flat)

        if self.tetas_cpa is None:
            sys.exit('EMTO.init_bulk(): \'self.tetas_cpa\' does not exist!!! (Did you run init_structure?)')
        else:
            tetas_flat = []
            for i in range(len(self.tetas_cpa)):
                if isinstance(self.tetas_cpa[i], list):
                    for j in range(len(self.tetas_cpa[i])):
                        tetas_flat.append(self.tetas_cpa[i][j])
                else:
                    tetas_flat.append(self.tetas_cpa[i])
            self.KGRN_tetas = np.array(tetas_flat)

        if self.phis_cpa is None:
            sys.exit('EMTO.init_bulk(): \'self.phis_cpa\' does not exist!!! (Did you run init_structure?)')
        else:
            phis_flat = []
            for i in range(len(self.phis_cpa)):
                if isinstance(self.phis_cpa[i], list):
                    for j in range(len(self.phis_cpa[i])):
                        phis_flat.append(self.phis_cpa[i][j])
                else:
                    phis_flat.append(self.phis_cpa[i])
            self.KGRN_phis = np.array(phis_flat)

        if self.asrs_cpa is None:
            sys.exit('EMTO.init_bulk(): \'self.asrs_cpa\' does not exist!!! (Did you run init_structure?)')
        else:
            asrs_flat = []
            for i in range(len(self.asrs_cpa)):
                if isinstance(self.asrs_cpa[i], list):
                    for j in range(len(self.asrs_cpa[i])):
                        asrs_flat.append(self.asrs_cpa[i][j])
                else:
                    asrs_flat.append(self.asrs_cpa[i])
            self.SHAPE_asrs = np.array(asrs_flat)

        if self.concs_cpa is None:
            sys.exit('EMTO.init_bulk(): \'self.concs_cpa\' does not exist!!! (Did you run init_structure?)')
        else:
            concs_flat = []
            for i in range(len(self.concs_cpa)):
                if isinstance(self.concs_cpa[i], list):
                    for j in range(len(self.concs_cpa[i])):
                        concs_flat.append(self.concs_cpa[i][j])
                else:
                    concs_flat.append(self.concs_cpa[i])
            self.KGRN_concs = np.array(concs_flat)

        # S_wss
        if self.s_wss_cpa is None:
            sys.exit('EMTO.init_bulk(): \'self.s_wss_cpa\' does not exist!!! (Did you run init_structure?)')
        else:
            s_wss_flat = []
            for i in range(len(self.s_wss_cpa)):
                if isinstance(self.s_wss_cpa[i], list):
                    for j in range(len(self.s_wss_cpa[i])):
                        s_wss_flat.append(self.s_wss_cpa[i][j])
                else:
                    s_wss_flat.append(self.s_wss_cpa[i])
            self.KGRN_s_wss = np.array(s_wss_flat)

        # ws_wsts
        if self.ws_wsts_cpa is None:
            sys.exit('EMTO.init_bulk(): \'self.ws_wsts_cpa\' does not exist!!! (Did you run init_structure?)')
        else:
            ws_wsts_flat = []
            for i in range(len(self.ws_wsts_cpa)):
                if isinstance(self.ws_wsts_cpa[i], list):
                    for j in range(len(self.ws_wsts_cpa[i])):
                        ws_wsts_flat.append(self.ws_wsts_cpa[i][j])
                else:
                    ws_wsts_flat.append(self.ws_wsts_cpa[i])
            self.KGRN_ws_wsts = np.array(ws_wsts_flat)

        if self.qtrs_cpa is None:
            sys.exit('EMTO.init_bulk(): \'self.qtrs_cpa\' does not exist!!! (Did you run init_structure?)')
        else:
            qtrs_flat = []
            for i in range(len(self.qtrs_cpa)):
                if isinstance(self.qtrs_cpa[i], list):
                    for j in range(len(self.qtrs_cpa[i])):
                        qtrs_flat.append(self.qtrs_cpa[i][j])
                else:
                    qtrs_flat.append(self.qtrs_cpa[i])
            self.KGRN_qtrs = np.array(qtrs_flat)

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
        self.use_sublattices = use_sublattices
        if its is None:
            if self.use_sublattices:
                # Use spglib calculated equivalent_atoms to assign sites
                # to sublattices.
                self.equivalent_its = np.array(self.equivalent_atoms) + 9000
                eq_dict = {}
                running_index = 0
                for i in range(len(self.equivalent_its)):
                    tmp = self.equivalent_its[i]
                    if tmp not in eq_dict.keys():
                        eq_dict[tmp] = running_index+1
                        running_index += 1

                for i in range(len(self.equivalent_its)):
                    self.equivalent_its[i] = eq_dict[self.equivalent_its[i]]

                # print(self.equivalent_its)

                self.KGRN_its = np.zeros(index_len, dtype='int32')
                running_index = 0
                for i in range(len(index_array)):
                    for j in range(index_array[i]):
                        self.KGRN_its[running_index] = self.equivalent_its[i]
                        running_index += 1
            else:
                running_index = 0
                for i in range(len(index_array)):
                    for j in range(index_array[i]):
                        self.KGRN_its[running_index] = i+1
                        running_index += 1
        else:
            its_flat = []
            for i in range(len(its)):
                if isinstance(its[i], list):
                    for j in range(len(its[i])):
                        its_flat.append(its[i][j])
                else:
                    its_flat.append(its[i])
            self.KGRN_its = np.array(its_flat, dtype='int32')

        # Pass structure related input settings to KSTR and SHAPE
        self.input_system.lattice.set_values(
            jobname_lat=self.latname,
            latpath=self.latpath,
            lat=common.ibz_to_lat(self.ibz),
            latparams=[1.0, self.output_boa, self.output_coa],
            latvectors=[self.output_alpha, self.output_beta, self.output_gamma],
            basis=self.output_basis,
            EMTOdir=self.EMTOdir,
            asrs=self.SHAPE_asrs,
            **kwargs)

        # Pass input settings to KGRN and KFCD
        self.input_system.bulk_new(lat=common.ibz_to_lat(self.ibz),
                                   ibz=self.ibz,
                                   latname=self.latname,
                                   latpath=self.latpath,
                                   atoms=self.KGRN_atoms,
                                   concs=self.KGRN_concs,
                                   m_splits=self.KGRN_m_splits,
                                   fxms=self.KGRN_fxms,
                                   iqs=self.KGRN_iqs,
                                   its=self.KGRN_its,
                                   itas=self.KGRN_itas,
                                   sws=self.sws,
                                   s_wss=self.KGRN_s_wss,
                                   ws_wsts=self.KGRN_ws_wsts,
                                   qtrs=self.KGRN_qtrs,
                                   nrms=self.KGRN_nrms,
                                   a_scrs=self.KGRN_a_scrs,
                                   b_scrs=self.KGRN_b_scrs,
                                   tetas=self.KGRN_tetas,
                                   phis=self.KGRN_phis,
                                   **kwargs)

    def write_bmdl_kstr_shape_input(self):
        self.input_system.lattice.write_structure_input_files(folder=self.folder, jobname_lat=self.latname)
        return

    def write_kgrn_kfcd_input(self):
        self.input_system.emto.kgrn.write_input_file(folder=self.folder)
        self.input_system.emto.kfcd.write_input_file(folder=self.folder)
        self.input_system.emto.batch.write_input_file(folder=self.folder)
        return

    def write_kgrn_kfcd_swsrange(self, sws=None):
        if sws is None:
            sys.exit('EMTO.write_KGRN_KFCD_swsrange(): An array of' +
                     ' WS-radii \'sws\' has to be given!')
        jobnames = self.input_system.lattice_constants_batch_generate(sws)
        return jobnames

    def draw_structure(self, which='input'):
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
