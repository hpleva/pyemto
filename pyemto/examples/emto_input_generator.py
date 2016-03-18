#!/usr/bin/python
import numpy as np
from pymatgen import Lattice, Structure
from pymatgen.vis.structure_vtk import StructureVis
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.structure_matcher import StructureMatcher

import os

# Help python to find the pyemto folder
import sys
sys.path.insert(0, "/home/henrik/local_emto_stuff/pyemto")

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

        self.sg2ibz = {8:13,12:13,42:11,160:7,166:7,225:2}
        self.sg2bl  = {8:'base-centered monoclinic',12:'base-centered monoclinic',
                       42:'face-centered orthorhombic',160:'rhombohedral/trigonal',
                       166:'rhombohedral/trigonal',225:'fcc'}

        # BMDL, KSTR, SHAPE, KGRN and KFCD class instances
        self.input_system = pyemto.System(folder=self.folder,EMTOdir=self.EMTOdir)
        return

    def calc_ws_radius(self,struct):
        bohr2angst = 0.52917721
        vol_unit = struct.volume/struct.num_sites
        sws = (3*vol_unit/4.0/np.pi)**(1.0/3)/bohr2angst
        return sws

    def make_basis_array(self,struct):
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

    def init_structure(self,prims=None,basis=None,atoms=None,latpath=None,
                       coords_are_cartesian=False,latname=None,kappaw=None):
        if prims is None:
            sys.exit('EMTO.init_structure(): \'prims\' has to be given!')
        if basis is None:
            sys.exit('EMTO.init_structure(): \'basis\' has to be given!')
        if atoms is None:
            sys.exit('EMTO.init_structure(): \'atoms\' has to be given!')
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
        self.stm = StructureMatcher(ltol=0.001,stol=0.001,angle_tol=0.01)
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
        if self.sg2ibz[self.finder_space] == 7:
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
            self.output_boa = 0.0
            self.output_coa = self.output_prima[2]
            self.output_alpha = 0.0
            self.output_beta = 0.0
            self.output_gamma = 0.0
            #print('emto_coa = ',emto_coa)
        elif self.sg2ibz[self.finder_space] == 11:
            self.primaa = prim_struct.lattice.matrix[0,:]
            self.primbb = prim_struct.lattice.matrix[1,:]
            self.primcc = prim_struct.lattice.matrix[2,:]
            self.output_basis = self.make_basis_array(self.prim_struct)
            lat_a = 2*self.primbb[0]
            self.output_prima = self.primbb/lat_a
            self.output_primb = self.primcc/lat_a
            self.output_primc = self.primaa/lat_a
            # Apply transformation on the basis atoms
            self.output_basis = self.output_basis/lat_a
            self.output_boa = 2*self.output_primc[1]
            self.output_coa = 2*self.output_primc[2]
            self.output_alpha = 0.0
            self.output_beta = 0.0
            self.output_gamma = 0.0
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
            self.KGRN_atoms = np.asarray(self.output_sites)
            index_array = np.ones(len(self.KGRN_atoms),dtype='int32')
            index_len = np.sum(index_array)
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
