# -*- coding: utf-8 -*-
"""

pyEmto v 0.9

Created on Wed Dec  3 14:25:06 2014

@author: Matti Ropo
@author: Henrik Levämäki
"""

from __future__ import print_function
import time
import os
import sys
import re
import numpy as np
import pyemto.common.common as common

class System:
    """The main class which provides the basis for the pyEMTO scripts.

    Somewhere in the beginning of a pyEMTO script a new instance of
    the system class should be created. All subsequent communication
    with the newly created system should be through the class methods,
    which are described below.

    :param folder: Main folder where the input and output files will
                   be stored. Use of absolute paths is recommended
                   (Default value = current working directory)
    :type folder: str
    :param EMTOdir: Path to the folder of the EMTO installation.
                    This entry can and should be modified by the user
                    inside the System.__init__ function
                    (Default value = /home/user/EMTO5.8)
    :type EMTOdir: str
    :param xc: Choice for the xc-functional can be set here.
               (Default value = PBE)
    :type xc: str
    :returns: None
    :rtype: None
    """

    def __init__(self, folder=None, EMTOdir=None, xc=None):
        # Import necessary packages
        from pyemto.latticeinputs.latticeinputs import Latticeinputs
        from pyemto.emtoinputs.emtoinputs import Emtoinputs

        # Check input arguments
        if folder is None:
            self.folder = os.getcwd()  # Use current folder
        else:
            self.folder = folder
        if EMTOdir is None:
            self.EMTOdir = "/home/hpleva/EMTO5.8"
        else:
            self.EMTOdir = EMTOdir

        # Initialize default parameters
        self.ca_range_default = np.linspace(1.50, 1.70, 7)
        self.elastic_constants_points = 6
        self.elastic_constants_deltas = np.linspace(0.0, 0.05,
                                        self.elastic_constants_points)
        self.RyBohr3_to_GPa = 14710.5065722
        self.kappaw_default = [0.0, -20.0]
        self.hcpo_relax_points = 5
        self.hcpm_relax_points = 5

        if xc is None:
            self.xc = 'PBE'
        else:
            self.xc = xc

        # Create working folders
        common.check_folders(folder + '/kgrn', folder + '/kgrn/tmp',
                           folder + '/kfcd', folder + '/fit')

        # BMDL, KSTR, SHAPE, KGRN and KFCD class instances
        self.lattice = Latticeinputs()
        self.emto = Emtoinputs()

        return

    def bulk(self, jobname=None, lat=None, atoms=None, concs=None, splts=None, sws=None,
             latname=None, latpath=None, emtopath=None, ibz=None, bmod=None, xc=None, ca=None,
             **kwargs):
        """Initializes the basic parameters for bulk systems.

        Basic information concerning the system,
        such as the types of atoms and the crystal structure should be given to this function and
        it should be called right after the class instance has been created.

        :param jobname: Name of the system (Default value = None)
        :type jobname:
        :param lat: The type of lattice structure (Default value = None)
        :type lat:
        :param atoms: List of atoms in the system (Default value = None)
        :type atoms:
        :param concs: List of concentrations of the elements in the
                      'atoms' list. This information is only used in CPA
                      calculations (Default value = None)
        :type concs:
        :param splts: List of initial magnetic moments of the elements in the
                      'atoms' list (Default value = None)
        :type splts:
        :param sws: The Wigner-Seitz radius of the system (Default value = None)
        :type sws: float
        :param latname: The 'jobname' of the BMDL, KSTR and SHAPE output files. These
                        structure output files have to be located in the 'latpath'
                        directory and they have to be named jobname.extention
                        (Default value = None)
        :type latname:
        :param latpath: The absolute path to the folder where the 'bmdl', 'kstr' and 'shape'
                        folders are located, which in turn contain the output files of the
                        structure calculation (Default value = None)
        :type latpath:
        :param emtopath: The absolute path to the folder where the EMTO installation is
                         located (Default value = None)
        :type emtopath:
        :param ibz: The code number indicating the Bravais lattice that the crystal
                    structure of the system has. For a list of possible values, please consult the
                    EMTO manual (Default value = None)
        :type ibz:
        :param bmod: The bulk modulus can be inputed here and if it is given,
                     it will be used by the elastic modulus routines (Default value = None)
        :type bmod:
        :param xc: The choice of the xc-functional. If None, PBE will be used as default
                   (Default value = None)
        :type xc:
        :param ca: The c/a ratio of hcp structures can be inputed here and if it is given,
                   it will be used by the elastic modulus routines (Default value = None)
        :type ca:
        :param **kwargs: Arbitrary other KGRN and KFCD input parameters can be given here
                         as keyword arguments. They will be passed down to the
                         self.emto.set_values() function
        :type **kwargs: str,int,float,list(str),list(int),list(float)
        :returns: None
        :rtype: None
        """

        if lat is None:
            sys.exit('System.bulk(): \'lat\' has to be given!')
        else:
            self.lat = lat
        if latname is None:
            self.latname = self.lat
        else:
            self.latname = latname
        if latpath is None:
            self.latpath = "./"
        else:
            self.latpath = latpath
        if emtopath is None:
            self.emtopath = self.folder
        else:
            self.emtopath = emtopath
        if atoms is None:
            sys.exit('System.bulk(): \'atoms\' has to be given!')
        else:
            self.atoms = atoms
        if concs is None:
            # Assume equal concentrations for each element
            self.concs = np.zeros(len(atoms))
            self.concs[:] = 1.0 / float(len(atoms))
        else:
            self.concs = concs
        if splts is None:
            self.splts = np.zeros(len(atoms))
        else:
            self.splts = np.asarray(splts)
        if sws is None:
            self.sws = 0.0
            #sys.exit('System.bulk(): \'sws\' has to be given!')
        else:
            self.sws = sws
        if jobname is None:
            self.jobname, self.fulljobname = self.create_jobname()
        else:
            self.jobname = jobname
            self.fulljobname = self.create_jobname(jobname)
        if ibz is None:
            self.ibz = common.lat_to_ibz(self.lat)
        else:
            self.ibz = ibz

        # Knowledge of the c/a lattice parameter for hcp systems
        if ca is not None:
            self.ca = ca
        else:
            self.ca = None

        # Knowledge of the xc-functional we want to use
        if xc is None:
            self.xc = 'PBE'
        else:
            self.xc = xc

        # Knowledge of the value of the bulk modulus, which
        # is mainly needed in the elastic constant functions
        self.bmod = bmod

        # hcp requires that we "double" the atoms array and
        # create a non-trivial iqs-array because hcp has a
        # non-trivial two-atom basis.
        if self.lat == 'hcp':
            self.atoms = np.array([self.atoms, self.atoms]).flatten()
            if concs is None:
                self.concs = np.zeros(len(self.atoms))
                self.concs[:] = 2.0 / float(len(self.atoms))
            else:
                self.concs = np.array([self.concs, self.concs]).flatten()
            self.iqs = np.zeros(len(self.atoms), dtype='int32')
            self.iqs[:len(self.iqs) / 2] = 1
            self.iqs[len(self.iqs) / 2:] = 2
            if splts is None:
                self.splts = np.zeros(len(self.atoms))
            else:
                self.splts = np.array([self.splts, self.splts]).flatten()
            self.itas = np.arange(1, len(self.atoms) / 2 + 1, dtype='int32')
            self.itas = np.array([self.itas, self.itas]).flatten()

            self.emto.set_values(jobname=self.fulljobname, sws=self.sws, atoms=self.atoms,
                                 iqs=self.iqs, itas=self.itas, concs=self.concs, splts=self.splts,
                                 ibz=self.ibz, latname=self.latname, latpath=self.latpath,
                                 emtopath=self.emtopath, EMTOdir=self.EMTOdir, **kwargs)

        # Special settings for the B2 structure. CPA currently not supported!!!
        elif self.lat == 'B2':
            self.iqs   = np.array([1,2],dtype='int32')
            self.concs = np.array([1.0,1.0])
            self.its   = np.array([1,2],dtype='int32')
            self.itas  = np.array([1,1],dtype='int32')
            
            self.emto.set_values(jobname=self.fulljobname, sws=self.sws, atoms=self.atoms,
                                 iqs=self.iqs, its=self.its, itas=self.itas, concs=self.concs,
                                 splts=self.splts, ibz=self.ibz, latname=self.latname, latpath=self.latpath,
                                 emtopath=self.emtopath, EMTOdir=self.EMTOdir, **kwargs)

        else:
            self.emto.set_values(jobname=self.fulljobname, sws=self.sws, atoms=self.atoms,
                                 concs=self.concs, splts=self.splts, ibz=self.ibz,
                                 latname=self.latname, latpath=self.latpath, emtopath=self.emtopath,
                                 EMTOdir=self.EMTOdir, **kwargs)
        return


    def write_inputs(self,folder=None, batch=True):
        """Write kgrn and kfcd inputs files and possible batch file """
        if folder == None:
            folder = self.folder
        self.emto.kgrn.write_input_file(folder=folder)
        self.emto.kfcd.write_input_file(folder=folder)
        if batch:
            self.emto.batch.write_input_file(folder=folder)

    
    def lattice_constants_analyze(self, sws=None, ca=None,prn=True,debug=False,method='morse',return_error=False):
        """Analyzes the output files generated using the
        lattice_constants_batch_generate function.

        The results are printed on screen.

        :param sws: List of WS-radii (Default value = None)
        :type sws: list(float)
        :param ca: List hpc c/a ratios  (Default value = None)
        :type ca: list(float)
        :param prn: True if results should be printed on screen, False if not (Default value = True)
        :type prn: boolean
        :returns: Equilibrium WS-radius, c/a (only hcp), bulk modulus, energy,
                  R (only hcp) and cs (only hcp)
        :rtype: float, float (only hcp), float, float, float (only hcp), float (only hcp)
        """

        from pyemto.EOS.EOS import EOS

        eos = EOS(name=self.jobname, xc=self.xc, method=method, units='bohr')

        if prn:
            print('')
            print('*****lattice_constants_analyze*****')
            print('')

        if sws is None:
            sys.exit('System.lattice_constants_analyze(): An array of' +
                     ' WS-radii \'sws\' has to be given!')
        else:
            self.lc_analyze_sws_range = np.asarray(sws)
        if ca is None:
            self.lc_analyze_ca_range = self.ca_range_default
        else:
            self.lc_analyze_ca_range = np.asarray(ca)

        if self.lat == 'bcc' or self.lat == 'fcc' or self.lat == 'trig' or self.lat == 'stric':
            energies = []
            swses = []

            for j in range(len(self.lc_analyze_sws_range)):
                self.sws = self.lc_analyze_sws_range[j]
                job = self.create_jobname(self.jobname)
                en = self.get_energy(job, folder=self.folder, func=self.xc)
                if isinstance(en, type(None)):
                    print('System.lattice_constants_analyze(): Warning:' +
                          ' No output energy found for {0}'.format(job))
                else:
                    energies.append(en)
                    swses.append(self.lc_analyze_sws_range[j])

            if prn:
                self.print_sws_ens(
                    'lattice_constants_analyze(cubic)', swses, energies)

            sws0, e0, B0, grun, R_squared = eos.fit(swses, energies)

            # These functions create files on disk about the data to be fitted
            # as well as the results of the fit.
            # eos.prepareData()
            #sws0,e0,B,grun = eos.fit2file()

            if prn:
                print('lattice_constants_analyze(cubic):')
                print('sws0 = {0:13.6f}'.format(sws0))
                print('B0   = {0:13.6f}'.format(B0))
                print('E0   = {0:13.6f}'.format(e0))
                print('')

            if return_error:
                return sws0, B0, e0, R_squared
            else:
                return sws0, B0, e0

        if self.lat == 'hcp':

            # Fit an n'th order polynomial to the energy vs. c/a data.
            ca_fit_order = 2

            #latnames = ['hcp_ca1', 'hcp_ca2', 'hcp_ca3', 'hcp_ca4', 'hcp_ca5', 'hcp_ca6',
            #            'hcp_ca7', 'hcp_ca8', 'hcp_ca9', 'hcp_ca10', 'hcp_ca11', 'hcp_ca12']

            caname = 'hcp_ca'
            # For the energies a 2D-array [i,j], where i = c/a axis and j = sws
            # axis
            energies = np.zeros(
                (len(self.lc_analyze_ca_range), len(self.lc_analyze_sws_range)))
            energies0 = []  # c/a optimized energy for a given WS-radius
            swses = []  # List of WS-radii
            cas0 = []  # Energetically optimized c/a's for a given WS-radius

            # First collect all the output energies into the 2D array.
            for i in range(len(self.lc_analyze_ca_range)):
                for j in range(len(self.lc_analyze_sws_range)):

                    self.sws = self.lc_analyze_sws_range[j]
                    job = self.create_jobname(self.jobname + "_" + caname +str(i+1))
                    en = self.get_energy(job, folder=self.folder, func=self.xc)
                    if isinstance(en, type(None)):
                        print('System.lattice_constants_analyze(): Warning:' +
                              ' No output energy found for {0}'.format(job))
                    else:
                        energies[i, j] = en

            if debug:
                print('Energy matrix (y axis = c/a axis, x axis = sws axis):'+"\n")

                formaatti = "                   "
                for i in range(len(self.lc_analyze_sws_range)):
                    formaatti = formaatti + "{0}{1}{2}{3}      ".format("{",i,":8.6f","}")
                print(formaatti.format(*self.lc_analyze_sws_range))

                formaatti = ""
                for i in range(len(self.lc_analyze_sws_range)+1):
                    formaatti = formaatti + "{0}{1}{2}{3} ".format("{",i,":13.6f","}")
                for i in range(len(energies[:,0])):
                    print(formaatti.format(self.lc_analyze_ca_range[i],*energies[i,:]))
                

            # Now we can start processing the 2D energy array.
            # There might be some calculations that didn't converge.
            # For those the energy will be zero in the 2D array so
            # we have to leave those points out.
            for j in range(len(self.lc_analyze_sws_range)):
                good_energies = []
                good_cas = []
                for i in range(len(self.lc_analyze_ca_range)):
                    if energies[i, j] < -1.0:
                        good_energies.append(energies[i, j])
                        good_cas.append(self.lc_analyze_ca_range[i])
                if len(good_energies) >= 3:
                    ca0, en0 = eos.ca_fit(good_cas, good_energies, ca_fit_order,
                                          debug=debug,title='sws{0}'.format(j+1))
                    cas0.append(ca0)
                    energies0.append(en0)
                    swses.append(self.lc_analyze_sws_range[j])

            if prn:
                self.print_sws_ens_hcp(
                    'lattice_constants_analyze(hcp)', swses, energies0, cas0)

            print('#'*80)
            print('# Ground state EOS fit:'+' '*56+'#')
            print('#'*80)
            sws0, e0, B0, grun, R_squared = eos.fit(swses, energies0)
            print('*'*80)
            print('*'*80)
            print('*'*80+'\n')

            # These functions create files on disk about the data to be fitted
            # as well as the results of the fit.
            # eos.prepareData()
            #sws0,e0,B0,grun = eos.fit2file()

            # Now that we have the ground state WS-radius sws0 we can use the cas0 array
            # to compute the corresponding ground state ca0
            ca0_vs_sws = np.polyfit(swses, cas0, 2)
            ca0_vs_sws = np.poly1d(ca0_vs_sws)
            c_over_a0 = ca0_vs_sws(sws0)

            # In order to calculate R = (c33 - c11 - c12 + c13) / cs
            # we use Eq. (6.56) p. 108 in Vitos' book.
            dca0_vs_dsws = np.polyder(ca0_vs_sws)
            dca0dsws0 = dca0_vs_dsws(sws0)
            R0 = -sws0 / 3.0 / c_over_a0 * dca0dsws0

            # In order to calculate cs = c11 + c12 + 2c33 - 4c13
            # we use Eq. (6.68) p. 109 in Vitos' book.
            # We need the total energies as a function of c/a
            # where sws is fixed to sws0. These can be obtained
            # by doing a Morse fit for each c/a and evaluating
            # the fitting function at sws0.

            energies_cs = []

            # There might be some calculations that didn't converge.
            # For those the energy will be zero in the 2D array so
            # we have to leave those points out.
            for i in range(len(self.lc_analyze_ca_range)):
                good_energies = []
                good_swses = []
                for j in range(len(self.lc_analyze_sws_range)):
                    if energies[i, j] < -1.0:
                        good_energies.append(energies[i, j])
                        good_swses.append(self.lc_analyze_sws_range[j])
                print('#'*80)
                print('# c/a{0} EOS fit:'.format(i+1)+' '*64+'#')
                print('#'*80)
                # _tmp variables are just dummies, we only want to
                # update the EOS parameters of the "eos" instance.
                sws_tmp, e_tmp, B_tmp, grun_tmp, R_squared_tmp = eos.fit(
                    good_swses, good_energies)
                print('*'*80)
                print('*'*80)
                print('*'*80+'\n')
                #e_cs = eos.fit_eval(sws0)
                #e_cs = eos.fit_eval(sws_tmp)
                e_cs = e_tmp
                energies_cs.append(e_cs)

            e_vs_ca_at_sws0 = np.polyfit(
                self.lc_analyze_ca_range, energies_cs, 2)
            e_vs_ca_at_sws0 = np.poly1d(e_vs_ca_at_sws0)

            d2e_vs_dca2_at_sws0 = np.polyder(e_vs_ca_at_sws0, 2)
            d2edca02sws0 = d2e_vs_dca2_at_sws0(c_over_a0)

            vol0 = 4.0 / 3.0 * np.pi * sws0**3
            cs0 = 9.0 / 2.0 * c_over_a0**2 / vol0 * \
                d2edca02sws0 * self.RyBohr3_to_GPa

            if prn:
                print('hcp_lattice_constants_analyze(hcp):')
                print('sws0 = {0:13.6f}'.format(sws0))
                print('c/a0 = {0:13.6f}'.format(c_over_a0))
                print('B0   = {0:13.6f}'.format(B0))
                print('E0   = {0:13.6f}'.format(e0))
                print('R    = {0:13.6f}'.format(R0))
                print('cs   = {0:13.6f}'.format(cs0))
                print('')

            if return_error:
                return sws0, c_over_a0, B0, e0, R0, cs0, R_squared
            else:
                return sws0, c_over_a0, B0, e0, R0, cs0

    def lattice_constants_batch_generate(self, sws=None, ca=None):
        """Generates input files and writes them to disk.

        Based on the input *sws* and *ca* lists jobnames are created and
        then corresponding input files are generated and written
        on disk. List of jobnames are returned.

        :param sws: List of WS-radii (Default value = None)
        :type sws: list(float)
        :param ca: List of hcp c/a ratios (Default value = None)
        :type ca: list(float)
        :returns: List of jobnames
        :rtype: list(str)
        """

        if sws is None:
            sys.exit('System.lattice_constants_batch(): An array of' +
                     ' WS-radii \'sws\' has to be given!')
        else:
            self.lc_batch_sws_range = np.asarray(sws)
        if ca is None:
            self.lc_batch_ca_range = self.ca_range_default
        else:
            self.lc_batch_ca_range = np.asarray(ca)

        jobnames = []

        if self.lat == 'bcc' or self.lat == 'fcc':
            for j in range(len(self.lc_batch_sws_range)):
                self.sws = self.lc_batch_sws_range[j]
                job = self.create_jobname(self.jobname)
                jobnames.append(job)
                self.emto.set_values(sws=self.sws, jobname=job)

                common.check_folders(
                    self.folder, self.folder + "/kgrn", self.folder + "/kgrn/tmp")
                common.check_folders(self.folder + "/kfcd")
                common.check_folders(self.folder + "/fit")

                self.emto.kgrn.write_input_file(folder=self.folder)
                self.emto.kfcd.write_input_file(folder=self.folder)
                self.emto.batch.write_input_file(folder=self.folder)

            return jobnames

        elif self.lat == 'trig' or self.lat == 'stric':
            for j in range(len(self.lc_batch_sws_range)):
                self.sws = self.lc_batch_sws_range[j]
                job = self.create_jobname(self.jobname)
                jobnames.append(job)
                self.emto.set_values(sws=self.sws, jobname=job)

                common.check_folders(
                    self.folder, self.folder + "/kgrn", self.folder + "/kgrn/tmp")
                common.check_folders(self.folder + "/kfcd")
                common.check_folders(self.folder + "/fit")

                self.emto.kgrn.write_input_file(folder=self.folder)
                self.emto.kfcd.write_input_file(folder=self.folder)
                self.emto.batch.write_input_file(folder=self.folder)

            return jobnames

        elif self.lat == 'hcp':
            #latnames = ['hcp_ca1', 'hcp_ca2', 'hcp_ca3', 'hcp_ca4', 'hcp_ca5', 'hcp_ca6',
            #            'hcp_ca7', 'hcp_ca8', 'hcp_ca9', 'hcp_ca10', 'hcp_ca11', 'hcp_ca12']

            for i in range(len(self.lc_batch_ca_range)):
                for j in range(len(self.lc_batch_sws_range)):

                    caname = "hcp_ca"+str(i+1)
                    self.sws = self.lc_batch_sws_range[j]
                    job = self.create_jobname(self.jobname + "_" + caname)
                    jobnames.append(job)
                    self.emto.set_values(sws=self.sws, jobname=job)
                    self.emto.set_values(latname=caname)

                    common.check_folders(
                        self.folder, self.folder + "/kgrn", self.folder + "/kgrn/tmp")
                    common.check_folders(self.folder + "/kfcd")
                    common.check_folders(self.folder + "/fit")

                    self.emto.kgrn.write_input_file(folder=self.folder)
                    self.emto.kfcd.write_input_file(folder=self.folder)
                    self.emto.batch.write_input_file(folder=self.folder)

            return jobnames

    def lattice_constants_batch_calculate(self, sws=None, ca=None):
        """Calculates the ground state WS-radius using the parallelism of
        the batch system by submitting one job for each entry in the *sws* list.

        This is a combination of the batch_generate and batch_analyze functions.
        At the end results are printed on screen.

        :param sws: List of WS-radii (Default value = None)
        :type sws: list(float)
        :param ca: hpc c/a ratio  (Default value = None)
        :type ca: float
        :returns: WS-radius, c/a (only hcp), bulk modulus, energy,
                  R (only hcp), cs (only hcp)
        :rtype: float, float (only hcp), float, float,
                float (only hcp), float (only hcp)
        """

        if sws is None:
            sys.exit('System.lattice_constants_calculate(): An array of' +
                     ' WS-radii \'sws\' has to be given!')
        else:
            self.lc_batch_sws_range = np.asarray(sws)
        if ca is None:
            self.lc_batch_ca_range = self.ca_range_default
        else:
            self.lc_batch_ca_range = np.asarray(ca)

        # Create input files
        jobnames = self.lattice_constants_batch_generate(
            sws=self.lc_batch_sws_range, ca=self.lc_batch_ca_range)

        # Submit calculations to the batch queue system
        jobIDs = self.submit_jobs(jobnames, folder=self.folder)

        # Wait until all the jobs have finished
        self.wait_for_jobs(jobIDs)

        # Now we can analyze the results
        if self.lat == 'bcc' or self.lat == 'fcc':
            sws0, B0, e0 = self.lattice_constants_analyze(
                sws=sws, ca=self.lc_batch_ca_range)
            return sws0, B0, e0
        elif self.lat == 'hcp':
            sws0, c_over_a0, B0, e0, R0, cs0 = self.lattice_constants_analyze(
                sws=sws, ca=self.lc_batch_ca_range)
            return sws0, c_over_a0, B0, e0, R0, cs0

    def lattice_constants_serial_calculate(self, sws=None, stype="simple", rerun=False, skip=False,
                                           delta=0.01, refine=True):
        """Calculates the equilibrium lattice constants on one CPU **WITHOUT** the batch system.

        Also the eq. bulk modulus, c/a ratio and energy are returned.

        :param sws: Initial guess for the eq. WS-radius (Default value = None)
        :type sws: float
        :param stype: Type of the energy minimisation algorithm (Default value = "simple")
        :type stype: str
        :param rerun: True if  (Default value = False)
        :type rerun: boolean
        :param skip: True if (Default value = False)
        :type skip: boolean
        :param delta: Step size (Default value = 0.01)
        :type delta: float
        :param refine: True if an optimized WS-radius vs. energy curve
                       should be computed, False if not (Default value = True)
        :type refine: boolean
        :returns: WS-radius, c/a, bulk modulus, energy
        :rtype: float, float, float, float
        """

        if sws is None:
            sys.exit('System.lattice_constants_serial_calculate():' +
                     ' Starting point \'sws\' has to be given!')
        else:
            self.lc_initial_sws = sws
        self.lc_stype = stype
        self.lc_rerun = rerun
        self.lc_skip = skip

        # Find initial volume and bulk module and refine those results.

        if self.lat == 'bcc' or self.lat == 'fcc':
            c_over_a = 0.0  # sc,bcc and fcc structures only need a
            print('Running self.finc_lc()')
            v, b, energy = self.find_lc(delta=delta, xc=self.xc)
            print('System.find_lc(): sws,Bmod,energy = ', v, b, energy)
            if refine:
                print('Running self.refine_lc()')
                v, b, energy = self.refine_lc(v, delta=delta, xc=self.xc)
            print('lattice_constants_serial_calculate: sws,c/a,Bmod,energy = ',
                  v, c_over_a, b, energy)
        elif self.lat == 'hcp':
            print('Running self.find_lc_hcp()')
            v, b, c_over_a, energy = self.find_lc_hcp(delta=delta, xc=self.xc)
            print('System.find_lc_hcp(): sws,c/a,Bmod,energy = ',
                  v, c_over_a, b, energy)
            if refine:
                print('Running self.refine_lc_hcp()')
                v, b, c_over_a, energy = self.refine_lc_hcp(
                    v, c_over_a, delta=delta, xc=self.xc)
            print('lattice_constants_serial_calculate: sws,c/a,Bmod,energy = ',
                  v, c_over_a, b, energy)
            print()

        return v, c_over_a, b, energy

    def elastic_constants_analyze(
            self, sws=None, bmod=None, ca=None, R=None, cs=None, relax=True, debug=False):
        """Analyzes the output files generated using the
        elastic_constants_batch_generate function.

        The results are printed on screen.

        :param sws: WS-radius (Default value = None)
        :type sws: float
        :param bmod: Bulk modulus (Default value = None)
        :type bmod: float
        :param ca: hpc c/a ratio  (Default value = None)
        :type ca: float
        :param R: Dimensionless quantity of hcp systems (Default value = None)
        :type R: float
        :param cs: Second order c/a-derivative of the energy (Default value = None)
        :type cs: float
        :returns: None
        :rtype: None
        """

        from pyemto.EOS.EOS import EOS

        eos = EOS(name=self.jobname, xc=self.xc, method='morse', units='bohr')

        # Mission critical parameters
        if sws is None:
            sys.exit(
                'System.elastic_constants_analyze(): \'sws\' has not been specified!')
        elif sws is not None:
            self.sws = sws

        if bmod is None:
            sys.exit(
                'System.elastic_constants_analyze(): \'sws\' has not been specified!')
        elif bmod is not None:
            self.bmod = bmod

        if ca is None and self.lat == 'hcp':
            sys.exit(
                'System.elastic_constants_analyze(): \'ca\' (c/a) has not been specified!')
        else:
            self.ec_analyze_ca = ca

        if R is None and self.lat == 'hcp':
            sys.exit(
                'System.elastic_constants_analyze(): \'R\' has not been specified!')
        else:
            self.ec_analyze_R = R

        if cs is None and self.lat == 'hcp':
            sys.exit(
                'System.elastic_constants_analyze(): \'cs\' has not been specified!')
        else:
            self.ec_analyze_cs = cs

        deltas = self.elastic_constants_deltas

        if self.lat == 'bcc' or self.lat == 'fcc' or self.lat == 'B2':
            # Orthorhombic distortion for c' first

            if self.lat == 'bcc':
                jobname_dist = [
                    '_bcco0', '_bcco1', '_bcco2', '_bcco3', '_bcco4', '_bcco5']
                latname_dist = [
                    'bcco0', 'bcco1', 'bcco2', 'bcco3', 'bcco4', 'bcco5']
            elif self.lat == 'fcc':
                jobname_dist = [
                    '_fcco0', '_fcco1', '_fcco2', '_fcco3', '_fcco4', '_fcco5']
                latname_dist = [
                    'fcco0', 'fcco1', 'fcco2', 'fcco3', 'fcco4', 'fcco5']
            elif self.lat == 'B2':
                jobname_dist = [
                    '_B2o0', '_B2o1', '_B2o2', '_B2o3', '_B2o4', '_B2o5']
                latname_dist = [
                    'B2o0', 'B2o1', 'B2o2', 'B2o3', 'B2o4', 'B2o5']

            en_cprime = []
            good_deltas_cprime = []

            for i in range(len(jobname_dist)):
                already = False
                job = self.create_jobname(self.jobname + jobname_dist[i])

                en = self.get_energy(job, folder=self.folder, func=self.xc)
                if isinstance(en, type(None)):
                    print('System.cubic_elastic_constants_analyze(): Warning:' +
                          ' No output energy found for {0}'.format(job))
                else:
                    en_cprime.append(en)
                    good_deltas_cprime.append(deltas[i])

            # Exit if we don't have enough points to do the fit
            if len(en_cprime) < 3:
                sys.exit(
                    'System.elastic_constants_analyze(): Not enough energy points to fit c\'!')

            # Convert into a numpy array
            en_cprime = np.asarray(en_cprime)
            good_deltas_cprime = np.asarray(good_deltas_cprime)

            # Next the monoclinic distortion for c44

            if self.lat == 'bcc':
                jobname_dist = [
                    '_bccm0', '_bccm1', '_bccm2', '_bccm3', '_bccm4', '_bccm5']
                latname_dist = [
                    'bccm0', 'bccm1', 'bccm2', 'bccm3', 'bccm4', 'bccm5']
            elif self.lat == 'fcc':
                jobname_dist = [
                    '_fccm0', '_fccm1', '_fccm2', '_fccm3', '_fccm4', '_fccm5']
                latname_dist = [
                    'fccm0', 'fccm1', 'fccm2', 'fccm3', 'fccm4', 'fccm5']
            elif self.lat == 'B2':
                jobname_dist = [
                    '_B2m0', '_B2m1', '_B2m2', '_B2m3', '_B2m4', '_B2m5']
                latname_dist = [
                    'B2m0', 'B2m1', 'B2m2', 'B2m3', 'B2m4', 'B2m5']

            en_c44 = []
            good_deltas_c44 = []

            for i in range(len(jobname_dist)):
                already = False
                job = self.create_jobname(self.jobname + jobname_dist[i])

                en = self.get_energy(job, folder=self.folder, func=self.xc)
                if isinstance(en, type(None)):
                    print('System.elastic_constants_analyze(): Warning:' +
                          ' No output energy found for {0}'.format(job))
                else:
                    en_c44.append(en)
                    good_deltas_c44.append(deltas[i])

            # Exit if we don't have enough points to do the fit
            if len(en_c44) < 3:
                sys.exit(
                    'System.elastic_constants_analyze(): Not enough energy points to fit c44!')

            # Convert into a numpy array
            en_c44 = np.asarray(en_c44)
            good_deltas_c44 = np.asarray(good_deltas_c44)

            # All calculations have been done, now it's time to fit the results
            popt_cprime, cprime_rsq = eos.distortion_fit(
                good_deltas_cprime, en_cprime,title='cprime')
            popt_c44, c44_rsq = eos.distortion_fit(
				good_deltas_c44, en_c44,title='c44')

            volume = 4.0 / 3.0 * np.pi * self.sws**3

            cprime = popt_cprime[0] / 2.0 / volume * self.RyBohr3_to_GPa
            c44 = popt_c44[0] / 2.0 / volume * self.RyBohr3_to_GPa

            c11 = self.bmod + 4.0 / 3.0 * cprime
            c12 = self.bmod - 2.0 / 3.0 * cprime

            # Polycrystalline elastic constants
            #
            # B = bulk modulus
            # G = shear modulus
            # E = Young modulus
            # v = Poisson ratio

            # Voigt average
            BV = (c11 + 2 * c12) / 3.0
            #BV = self.bmod
            GV = (c11 - c12 + 3 * c44) / 5.0
            EV = 9 * BV * GV / (3 * BV + GV)
            vV = (3 * BV - 2 * GV) / (6 * BV + 2 * GV)

            # Reuss average
            BR = BV
            #BR = self.bmod
            GR = 5 * (c11 - c12) * c44 / (4 * c44 + 3 * (c11 - c12))
            ER = 9 * BR * GR / (3 * BR + GR)
            vR = (3 * BR - 2 * GR) / (6 * BR + 2 * GR)

            # Hill average
            BH = (BV + BR) / 2.0
            #BH = self.bmod
            GH = (GV + GR) / 2.0
            EH = 9 * BH * GH / (3 * BH + GH)
            vH = (3 * BH - 2 * GH) / (6 * BH + 2 * GH)

            # Elastic anisotropy
            AVR = (GV - GR) / (GV + GR)

            print("")
            print('***cubic_elastic_constants***')
            print("")
            print(self.jobname)
            print("")
            print('sws(bohr)      = {0:7.3f}'.format(self.sws))
            print('B(GPa)         = {0:6.2f}'.format(self.bmod))
            print('c11(GPa)       = {0:6.2f}'.format(c11))
            print('c12(GPa)       = {0:6.2f}'.format(c12))
            print('c\'(GPa)        = {0:6.2f}'.format(cprime))
            print('c44(GPa)       = {0:6.2f}'.format(c44))
            print('R-squared(c\')  = {0:8.6f}'.format(cprime_rsq))
            print('R-squared(c44) = {0:8.6f}'.format(c44_rsq))
            print("")
            print('Voigt average:')
            print("")
            print('BV(GPa)  = {0:6.2f}'.format(BV))
            print('GV(GPa)  = {0:6.2f}'.format(GV))
            print('EV(GPa)  = {0:6.2f}'.format(EV))
            print('vV(GPa)  = {0:6.2f}'.format(vV))
            print("")
            print('Reuss average:')
            print("")
            print('BR(GPa)  = {0:6.2f}'.format(BR))
            print('GR(GPa)  = {0:6.2f}'.format(GR))
            print('ER(GPa)  = {0:6.2f}'.format(ER))
            print('vR(GPa)  = {0:6.2f}'.format(vR))
            print("")
            print('Hill average:')
            print("")
            print('BH(GPa)  = {0:6.2f}'.format(BH))
            print('GH(GPa)  = {0:6.2f}'.format(GH))
            print('EH(GPa)  = {0:6.2f}'.format(EH))
            print('vH(GPa)  = {0:6.2f}'.format(vH))
            print("")
            print('Elastic anisotropy:')
            print("")
            print('AVR(GPa)  = {0:6.2f}'.format(AVR))

            return

        elif self.lat == 'hcp':
            # Orthorhombic distortion for c66 first

            if relax:
                jobname_dist = []
                latname_dist = []
                for i in range(self.elastic_constants_points):
                    tmp_array1 = []
                    tmp_array2 = []
                    for j in range(self.hcpo_relax_points):
                        tmp_str1 = 'hcpo{0}_r{1}'.format(i,j)
                        tmp_str2 = '_hcpo{0}_r{1}'.format(i,j)
                        tmp_array1.append(tmp_str1)
                        tmp_array2.append(tmp_str2)
                        # We don't need to relax the first structure
                        # because it's undistorted
                        if i == 0 and j == 0:
                            break
                    latname_dist.append(tmp_array1)
                    jobname_dist.append(tmp_array2)

            else:
                jobname_dist = ['_hcpo0_ur',
                                '_hcpo1_ur',
                                '_hcpo2_ur',
                                '_hcpo3_ur',
                                '_hcpo4_ur',
                                '_hcpo5_ur']
                latname_dist = ['hcpo0_ur',
                                'hcpo1_ur',
                                'hcpo2_ur',
                                'hcpo3_ur',
                                'hcpo4_ur',
                                'hcpo5_ur']

            en_c66 = []
            good_deltas_c66 = []

            if relax:
                # First we have to determine how relaxation affects
                # the energy. The first structure corresponds to the
                # undistorted lattice, which means we don't need to
                # relax it:
                job = self.create_jobname(self.jobname + jobname_dist[0][0])
                        
                en_orig = self.get_energy(job, folder=self.folder, func=self.xc)
                if isinstance(en_orig, type(None)):
                    print('System.elastic_constants_analyze(): Warning:' +
                          ' No output energy found for {0}'.format(job))
                else:
                    en_c66.append(en_orig)
                    good_deltas_c66.append(deltas[0])
                
                # For all the rest of the structures we compute the
                # relaxed ground state energy by fitting a polynomial
                # to the atomic pos. vs. energy data:
                for i in range(1,self.elastic_constants_points):
                    ens_tmp = []
                    xind_tmp = []
                    for j in range(len(jobname_dist[i])):
                        already = False
                        job = self.create_jobname(self.jobname + jobname_dist[i][j])
                        
                        en = self.get_energy(job, folder=self.folder, func=self.xc)
                        if isinstance(en, type(None)):
                            print('System.elastic_constants_analyze(): Warning:' +
                                  ' No output energy found for {0}'.format(job))
                        else:
                            ens_tmp.append(en)
                            xind_tmp.append(j)
                    # Now we should find the relaxed energy by fitting
                    # a polynomial to the energy vs. atom coordinates data.
                    # Leave out the point if not enough calculations have
                    # converged.
                    if len(ens_tmp) >= 3:
                        relax_xmin, relax_emin = eos.relax_fit(xind_tmp,ens_tmp,2,debug=debug,
                                                               title='c66: hcpo{0}'.format(i))
                        en_c66.append(relax_emin)
                        # TEST: Use unrelaxed energies
                        #en_c66.append(ens_tmp[0])
                        good_deltas_c66.append(deltas[i])

            else:
                for i in range(self.elastic_constants_points):
                    already = False
                    job = self.create_jobname(self.jobname + jobname_dist[i])
                    
                    en = self.get_energy(job, folder=self.folder, func=self.xc)
                    if isinstance(en, type(None)):
                        print('System.elastic_constants_analyze(): Warning:' +
                                  ' No output energy found for {0}'.format(job))
                    else:
                        en_c66.append(en)
                        good_deltas_c66.append(deltas[i])

            # Exit if we don't have enough points to do the fit
            if len(en_c66) < 3:
                sys.exit('System.elastic_constants_analyze():' +
                         ' Not enough energy points to fit hcp c66!')

            # Convert into a numpy array
            en_c66 = np.asarray(en_c66)
            good_deltas_c66 = np.asarray(good_deltas_c66)

            # Next the monoclinic distortion for c44

            # For the monoclinic distortion relaxation effects
            # are negligible:
            relax = False
            if relax:
                jobname_dist = []
                latname_dist = []
                for i in range(self.elastic_constants_points):
                    tmp_array1 = []
                    tmp_array2 = []
                    for j in range(self.hcpo_relax_points):
                        tmp_str1 = 'hcpm{0}_r{1}'.format(i,j)
                        tmp_str2 = '_hcpm{0}_r{1}'.format(i,j)
                        tmp_array1.append(tmp_str1)
                        tmp_array2.append(tmp_str2)
                    latname_dist.append(tmp_array1)
                    jobname_dist.append(tmp_array2)

            else:
                jobname_dist = ['_hcpm0_ur',
                                '_hcpm1_ur',
                                '_hcpm2_ur',
                                '_hcpm3_ur',
                                '_hcpm4_ur',
                                '_hcpm5_ur']
                latname_dist = ['hcpm0_ur',
                                'hcpm1_ur',
                                'hcpm2_ur',
                                'hcpm3_ur',
                                'hcpm4_ur',
                                'hcpm5_ur']

            en_c44 = []
            good_deltas_c44 = []

            if relax:
                # First we have to determine how relaxation affects
                # the energy.
                
                for i in range(self.elastic_constants_points):
                    ens_tmp = []
                    xind_tmp = []
                    for j in range(self.hcpo_relax_points):
                        already = False
                        job = self.create_jobname(self.jobname + jobname_dist[i][j])
                        
                        en = self.get_energy(job, folder=self.folder, func=self.xc)
                        if isinstance(en, type(None)):
                            print('System.elastic_constants_analyze(): Warning:' +
                                  ' No output energy found for {0}'.format(job))
                        else:
                            ens_tmp.append(en)
                            xind_tmp.append(j)
                    # Now we should find the relaxed energy by fitting
                    # a polynomial to the energy vs. atom coordinates data.
                    # Leave out the point if not enough calculations have
                    # converged.
                    if len(ens_tmp) >= 3:
                        relax_xmin, relax_emin = eos.relax_fit(xind_tmp,ens_tmp,3,debug=debug,
                                                               title='c44: hcpm{0}'.format(i))
                        en_c44.append(relax_emix)
                        good_deltas_c44.append(deltas[i])

            else:
                for i in range(self.elastic_constants_points):
                    already = False
                    job = self.create_jobname(self.jobname + jobname_dist[i])
                    
                    en = self.get_energy(job, folder=self.folder, func=self.xc)
                    if isinstance(en, type(None)):
                        print('System.elastic_constants_analyze(): Warning:' +
                                  ' No output energy found for {0}'.format(job))
                    else:
                        en_c44.append(en)
                        good_deltas_c44.append(deltas[i])

            # Exit if we don't have enough points to do the fit
            if len(en_c44) < 3:
                sys.exit('System.elastic_constants_analyze():' +
                         ' Not enough energy points to fit hcp c44!')

            # Convert into a numpy array
            en_c44 = np.asarray(en_c44)
            good_deltas_c44 = np.asarray(good_deltas_c44)

            # All calculations have been done, now it's time to fit the results
            popt_c66, c66_rsq = eos.distortion_fit(good_deltas_c66, en_c66, title='c66')
            popt_c44, c44_rsq = eos.distortion_fit(good_deltas_c44, en_c44, title='c44')

            volume = 4.0 / 3.0 * np.pi * self.sws**3

            c66 = popt_c66[0] / 2.0 / volume * self.RyBohr3_to_GPa
            c44 = popt_c44[0] / 2.0 / volume * self.RyBohr3_to_GPa

            c11 = self.bmod + c66 + self.ec_analyze_cs * \
                (2 * self.ec_analyze_R - 1)**2 / 18.0
            c12 = self.bmod - c66 + self.ec_analyze_cs * \
                (2 * self.ec_analyze_R - 1)**2 / 18.0
            c13 = self.bmod + 1.0 / 9.0 * self.ec_analyze_cs * \
                (2 * self.ec_analyze_R**2 + self.ec_analyze_R - 1)
            c33 = self.bmod + 2.0 / 9.0 * \
                self.ec_analyze_cs * (self.ec_analyze_R + 1)**2
            c2 = c33 * (c11 + c12) - 2.0 * c13**2

            # Polycrystalline elastic constants
            #
            # B = bulk modulus
            # G = shear modulus
            # E = Young modulus
            # v = Poisson ratio

            # Voigt average
            BV = (2 * c11 + 2 * c12 + 4 * c13 + c33) / 9.0
            GV = (12 * c44 + 12 * c66 + self.ec_analyze_cs) / 30.0
            EV = 9 * BV * GV / (3 * BV + GV)
            vV = (3 * BV - 2 * GV) / (6 * BV + 2 * GV)

            # Reuss average
            BR = self.bmod
            GR = 5.0 / 2.0 * (c44 * c66 * c2) / \
                ((c44 + c66) * c2 + 3.0 * BV * c44 * c66)
            ER = 9 * BR * GR / (3 * BR + GR)
            vR = (3 * BR - 2 * GR) / (6 * BR + 2 * GR)

            # Hill average
            BH = (BV + BR) / 2.0
            #BH = self.bmod
            GH = (GV + GR) / 2.0
            EH = 9 * BH * GH / (3 * BH + GH)
            vH = (3 * BH - 2 * GH) / (6 * BH + 2 * GH)

            # Elastic anisotropy
            AVR = (GV - GR) / (GV + GR)

            print("")
            print('***hcp_elastic_constants***')
            print("")
            print(self.jobname)
            print("")
            print('sws(bohr)      = {0:7.3f}'.format(self.sws))
            print('B(GPa)         = {0:6.2f}'.format(self.bmod))
            print('c11(GPa)       = {0:6.2f}'.format(c11))
            print('c12(GPa)       = {0:6.2f}'.format(c12))
            print('c13(GPa)       = {0:6.2f}'.format(c13))
            print('c33(GPa)       = {0:6.2f}'.format(c33))
            print('c44(GPa)       = {0:6.2f}'.format(c44))
            print('c66(GPa)       = {0:6.2f}'.format(c66))
            print('R-squared(c44) = {0:8.6f}'.format(c44_rsq))
            print('R-squared(c66) = {0:8.6f}'.format(c66_rsq))
            print("")
            print('Voigt average:')
            print("")
            print('BV(GPa)  = {0:6.2f}'.format(BV))
            print('GV(GPa)  = {0:6.2f}'.format(GV))
            print('EV(GPa)  = {0:6.2f}'.format(EV))
            print('vV(GPa)  = {0:6.2f}'.format(vV))
            print("")
            print('Reuss average:')
            print("")
            print('BR(GPa)  = {0:6.2f}'.format(BR))
            print('GR(GPa)  = {0:6.2f}'.format(GR))
            print('ER(GPa)  = {0:6.2f}'.format(ER))
            print('vR(GPa)  = {0:6.2f}'.format(vR))
            print("")
            print('Hill average:')
            print("")
            print('BH(GPa)  = {0:6.2f}'.format(BH))
            print('GH(GPa)  = {0:6.2f}'.format(GH))
            print('EH(GPa)  = {0:6.2f}'.format(EH))
            print('vH(GPa)  = {0:6.2f}'.format(vH))
            print("")
            print('Elastic anisotropy:')
            print("")
            print('AVR(GPa)  = {0:6.2f}'.format(AVR))

            return

    def elastic_constants_batch_generate(self, sws=None, ca=None, relax=True):
        """Generates all the necessary input files based on the class data.

        :param sws: WS-radius (Default value = None)
        :type sws: float
        :param ca: hpc c/a ratio  (Default value = None)
        :type ca: float
        :returns: List of jobnames
        :rtype: list(str)
        """

        # Mission critical parameters
        if sws is None and self.sws is None:
            sys.exit(
                'System.elastic_constants_batch_generate(): \'sws\' has not been specified!')
        elif sws is not None:
            self.sws = sws

        if ca is None and self.lat == 'hcp':
            sys.exit('System.elastic_constants_batch_generate():' +
                     ' \'ca\' (c/a) has not been specified!')
        elif ca is not None:
            self.ca = ca

        jobnames = []
        if self.lat == 'B2':
            # Orthorhombic distortion input files for c' first

            jobname_dist = ['_B2o0','_B2o1','_B2o2','_B2o3','_B2o4','_B2o5']
            latname_dist = ['B2o0','B2o1','B2o2','B2o3','B2o4','B2o5']
            self.emto.set_values(ibz=8,nkx=17,nky=17,nkz=17) # High-quality
            #self.emto.set_values(ibz=10,nkx=15,nky=15,nkz=15) # Normal quality

            for i in range(len(jobname_dist)):
                job = self.create_jobname(self.jobname + jobname_dist[i])
                jobnames.append(job)
                self.emto.set_values(
                    sws=self.sws, jobname=job, latname=latname_dist[i])

                common.check_folders(
                    self.folder, self.folder + "/kgrn", self.folder + "/kgrn/tmp")
                common.check_folders(self.folder + "/kfcd")
                common.check_folders(self.folder + "/fit")
                self.emto.kgrn.write_input_file(folder=self.folder)
                self.emto.kfcd.write_input_file(folder=self.folder)
                self.emto.batch.write_input_file(folder=self.folder)

            # Next produce the input files of monoclinic distortion for c44

            jobname_dist = ['_B2m0','_B2m1','_B2m2','_B2m3','_B2m4','_B2m5']
            latname_dist = ['B2m0','B2m1','B2m2','B2m3','B2m4','B2m5']
            self.emto.set_values(ibz=9,nkx=17,nky=17,nkz=17) # High-quality
            #self.emto.set_values(ibz=11,nkx=15,nky=15,nkz=21) # Normal quality

            for i in range(len(jobname_dist)):
                job = self.create_jobname(self.jobname + jobname_dist[i])
                jobnames.append(job)
                self.emto.set_values(
                    sws=self.sws, jobname=job, latname=latname_dist[i])

                self.emto.kgrn.write_input_file(folder=self.folder)
                self.emto.kfcd.write_input_file(folder=self.folder)
                self.emto.batch.write_input_file(folder=self.folder)

            return jobnames

        if self.lat == 'bcc' or self.lat == 'fcc':
            # Orthorhombic distortion input files for c' first

            if self.lat == 'bcc':
                jobname_dist = ['_bcco0','_bcco1','_bcco2','_bcco3','_bcco4','_bcco5']
                latname_dist = ['bcco0','bcco1','bcco2','bcco3','bcco4','bcco5']
                self.emto.set_values(ibz=10,nkx=31,nky=31,nkz=31) # High-quality
                #self.emto.set_values(ibz=10,nkx=15,nky=15,nkz=15) # Normal quality
            elif self.lat == 'fcc':
                jobname_dist = ['_fcco0','_fcco1','_fcco2','_fcco3','_fcco4','_fcco5']
                latname_dist = ['fcco0','fcco1','fcco2','fcco3','fcco4','fcco5']
                #self.emto.set_values(ibz=11,nkx=41,nky=41,nkz=41) # Super-high-quality
                self.emto.set_values(ibz=11,nkx=31,nky=31,nkz=31) # High-quality
                #self.emto.set_values(ibz=11,nkx=17,nky=17,nkz=17) # Normal-quality

            for i in range(len(jobname_dist)):
                job = self.create_jobname(self.jobname + jobname_dist[i])
                jobnames.append(job)
                self.emto.set_values(
                    sws=self.sws, jobname=job, latname=latname_dist[i])

                common.check_folders(
                    self.folder, self.folder + "/kgrn", self.folder + "/kgrn/tmp")
                common.check_folders(self.folder + "/kfcd")
                common.check_folders(self.folder + "/fit")
                self.emto.kgrn.write_input_file(folder=self.folder)
                self.emto.kfcd.write_input_file(folder=self.folder)
                self.emto.batch.write_input_file(folder=self.folder)

            # Next produce the input files of monoclinic distortion for c44

            if self.lat == 'bcc':
                jobname_dist = ['_bccm0','_bccm1','_bccm2','_bccm3','_bccm4','_bccm5']
                latname_dist = ['bccm0','bccm1','bccm2','bccm3','bccm4','bccm5']
                self.emto.set_values(ibz=11,nkx=27,nky=27,nkz=37) # High-quality
                #self.emto.set_values(ibz=11,nkx=15,nky=15,nkz=21) # Normal quality
            elif self.lat == 'fcc':
                jobname_dist = ['_fccm0','_fccm1','_fccm2','_fccm3','_fccm4','_fccm5']
                latname_dist = ['fccm0','fccm1','fccm2','fccm3','fccm4','fccm5']
                #self.emto.set_values(ibz=10,nkx=37,nky=53,nkz=37) # Super-high-quality
                self.emto.set_values(ibz=10,nkx=27,nky=37,nkz=27) # High-quality
                #self.emto.set_values(ibz=10,nkx=15,nky=21,nkz=15) # Normal quality

            for i in range(len(jobname_dist)):
                job = self.create_jobname(self.jobname + jobname_dist[i])
                jobnames.append(job)
                self.emto.set_values(
                    sws=self.sws, jobname=job, latname=latname_dist[i])

                self.emto.kgrn.write_input_file(folder=self.folder)
                self.emto.kfcd.write_input_file(folder=self.folder)
                self.emto.batch.write_input_file(folder=self.folder)

            return jobnames

        elif self.lat == 'hcp':

            # Orthorhombic distortion input files for c66 first

            # With hcp the structure depends on the c/a ratio. Therefore we also have
            # to generate the corresponding structure files.

            if relax:
                jobname_dist = []
                latname_dist = []
                for i in range(self.elastic_constants_points):
                    tmp_array1 = []
                    tmp_array2 = []
                    for j in range(self.hcpo_relax_points):
                        tmp_str1 = 'hcpo{0}_r{1}'.format(i,j)
                        tmp_str2 = '_hcpo{0}_r{1}'.format(i,j)
                        tmp_array1.append(tmp_str1)
                        tmp_array2.append(tmp_str2)
                        # We don't need to relax the first structure
                        # because it's undistorted
                        if i == 0 and j == 0:
                            break
                    latname_dist.append(tmp_array1)
                    jobname_dist.append(tmp_array2)
            else:
                jobname_dist = ['_hcpo0_ur', '_hcpo1_ur',
                                '_hcpo2_ur', '_hcpo3_ur', '_hcpo4_ur', '_hcpo5_ur']
                latname_dist = ['hcpo0_ur', 'hcpo1_ur',
                                'hcpo2_ur', 'hcpo3_ur', 'hcpo4_ur', 'hcpo5_ur']

            # Check whether Two-center Taylor expansion is on/off
            if self.emto.kgrn.expan == 'M':
                kappaw = self.kappaw_default
                self.lattice.set_values(kappaw=kappaw)

            common.check_folders(self.folder)
            common.check_folders(self.folder + '/bmdl')
            common.check_folders(self.folder + '/kstr')
            common.check_folders(self.folder + '/shape')

            if relax:
                for i in range(self.elastic_constants_points):
                    for j in range(len(jobname_dist[i])):
                        if i == 0:
                            do_i_relax = False
                        else:
                            do_i_relax = True

                        self.lattice.distortion(lat='hcp', dist='ortho', ca=self.ca, index=i,
                                                deltas=self.elastic_constants_deltas,
                                                relax=do_i_relax,relax_index=j)
                        
                        self.lattice.set_values(jobname=latname_dist[i][j],latpath=self.folder)
                        self.lattice.bmdl.write_input_file(folder=self.folder)
                        self.lattice.kstr.write_input_file(folder=self.folder)
                        self.lattice.shape.write_input_file(folder=self.folder)
                        self.lattice.batch.write_input_file(folder=self.folder)

            else:
                for i in range(self.elastic_constants_points):
                    self.lattice.distortion(lat='hcp', dist='ortho', ca=self.ca, index=i,
                                            deltas=self.elastic_constants_deltas,
                                            relax=False)
                    
                    self.lattice.set_values(jobname=latname_dist[i],latpath=self.folder)
                    self.lattice.bmdl.write_input_file(folder=self.folder)
                    self.lattice.kstr.write_input_file(folder=self.folder)
                    self.lattice.shape.write_input_file(folder=self.folder)
                    self.lattice.batch.write_input_file(folder=self.folder)

            #self.emto.set_values(ibz=9, nkx=25, nky=25, nkz=21)
            self.emto.set_values(ibz=9, nkx=31, nky=31, nkz=25)

            if relax:
                for i in range(self.elastic_constants_points):
                    for j in range(len(jobname_dist[i])):

                        job = self.create_jobname(self.jobname + jobname_dist[i][j])
                        jobnames.append(job)
                        self.emto.set_values(sws=self.sws, jobname=job, latname=latname_dist[i][j],
                                             latpath=self.folder)
                        
                        common.check_folders(self.folder, self.folder + "/kgrn")
                        common.check_folders(self.folder + "/kgrn/tmp")
                        common.check_folders(self.folder + "/kfcd")
                        common.check_folders(self.folder + "/fit")
                        self.emto.kgrn.write_input_file(folder=self.folder)
                        self.emto.kfcd.write_input_file(folder=self.folder)
                        self.emto.batch.write_input_file(folder=self.folder)
                        
            else:
                for i in range(self.elastic_constants_points):
                    job = self.create_jobname(self.jobname + jobname_dist[i])
                    jobnames.append(job)
                    self.emto.set_values(sws=self.sws, jobname=job, latname=latname_dist[i],
                                         latpath=self.folder)
                    
                    common.check_folders(self.folder, self.folder + "/kgrn")
                    common.check_folders(self.folder + "/kgrn/tmp")
                    common.check_folders(self.folder + "/kfcd")
                    common.check_folders(self.folder + "/fit")
                    self.emto.kgrn.write_input_file(folder=self.folder)
                    self.emto.kfcd.write_input_file(folder=self.folder)
                    self.emto.batch.write_input_file(folder=self.folder)


            # Monoclinic distortion input files for c44 next
            # For the monoclinic distortion relaxation effects are
            # negligibly small; therefore
            relax = False

            # With hcp the structure depends on the c/a ratio. Therefore we also have
            # to generate the corresponding structure files.

            if relax:
                jobname_dist = []
                latname_dist = []
                for i in range(self.elastic_constants_points):
                    tmp_array1 = []
                    tmp_array2 = []
                    for j in range(self.hcpm_relax_points):
                        tmp_str1 = 'hcpm{0}_r{1}'.format(i,j)
                        tmp_str2 = '_hcpm{0}_r{1}'.format(i,j)
                        tmp_array1.append(tmp_str1)
                        tmp_array2.append(tmp_str2)
                        # We don't need to relax the first structure
                        # because it's undistorted
                        if i == 0 and j == 0:
                            break
                    latname_dist.append(tmp_array1)
                    jobname_dist.append(tmp_array2)
            else:
                jobname_dist = ['_hcpm0_ur', '_hcpm1_ur',
                                '_hcpm2_ur', '_hcpm3_ur', '_hcpm4_ur', '_hcpm5_ur']
                latname_dist = ['hcpm0_ur', 'hcpm1_ur',
                                'hcpm2_ur', 'hcpm3_ur', 'hcpm4_ur', 'hcpm5_ur']

            # Check whether Two-center Taylor expansion is on/off
            if self.emto.kgrn.expan == 'M':
                kappaw = self.kappaw_default
                self.lattice.set_values(kappaw=kappaw)

            common.check_folders(self.folder)
            common.check_folders(self.folder + '/bmdl')
            common.check_folders(self.folder + '/kstr')
            common.check_folders(self.folder + '/shape')

            if relax:
                for i in range(self.elastic_constants_points):
                    for j in range(len(jobname_dist[i])):

                        if i == 0:
                            do_i_relax = False
                        else:
                            do_i_relax = True

                        self.lattice.distortion(lat='hcp', dist='mono', ca=self.ca, index=i,
                                                deltas=self.elastic_constants_deltas,
                                                relax=do_i_relax,relax_index=j)
                        
                        self.lattice.set_values(jobname=latname_dist[i][j],latpath=self.folder)
                        self.lattice.bmdl.write_input_file(folder=self.folder)
                        self.lattice.kstr.write_input_file(folder=self.folder)
                        self.lattice.shape.write_input_file(folder=self.folder)
                        self.lattice.batch.write_input_file(folder=self.folder)

            else:
                for i in range(self.elastic_constants_points):
                    self.lattice.distortion(lat='hcp', dist='mono', ca=self.ca, index=i,
                                            deltas=self.elastic_constants_deltas,
                                            relax=False)
                    
                    self.lattice.set_values(jobname=latname_dist[i],latpath=self.folder)
                    self.lattice.bmdl.write_input_file(folder=self.folder)
                    self.lattice.kstr.write_input_file(folder=self.folder)
                    self.lattice.shape.write_input_file(folder=self.folder)
                    self.lattice.batch.write_input_file(folder=self.folder)

            # Two-atom basis
            #self.emto.set_values(ibz=13, nkx=41, nky=41, nkz=41)

            if relax:
                for i in range(self.elastic_constants_points):
                    for j in range(len(jobname_dist[i])):

                        job = self.create_jobname(self.jobname + jobname_dist[i][j])
                        jobnames.append(job)
                        self.emto.set_values(sws=self.sws, jobname=job, latname=latname_dist[i][j],
                                             latpath=self.folder)
                        
                        common.check_folders(self.folder, self.folder + "/kgrn")
                        common.check_folders(self.folder + "/kgrn/tmp")
                        common.check_folders(self.folder + "/kfcd")
                        common.check_folders(self.folder + "/fit")
                        self.emto.kgrn.write_input_file(folder=self.folder)
                        self.emto.kfcd.write_input_file(folder=self.folder)
                        self.emto.batch.write_input_file(folder=self.folder)

            else:
                # Four-atom basis
                #self.emto.set_values(ibz=12, nkx=18, nky=18, nkz=12)
                self.emto.set_values(ibz=12, nkx=22, nky=22, nkz=14)
                ###################################################################
                # Atconf related arrays need to be modified because we now have   #
                # a four atom basis.                                              #
                ###################################################################

                self.atoms = np.array([self.atoms, self.atoms]).flatten()
                self.concs = np.array([self.concs, self.concs]).flatten()

                self.iqs = np.zeros(len(self.atoms), dtype='int32')
                len_div = len(self.iqs) // 4
                for i in range(4):
                    self.iqs[i * len_div:(i + 1) * len_div] = i + 1

                self.splts = np.array([self.splts, self.splts]).flatten()
                self.itas = np.array([self.itas, self.itas]).flatten()
                
                self.emto.set_values(atoms=self.atoms, iqs=self.iqs, itas=self.itas,
                                     concs=self.concs, splts=self.splts)
                ####################################################################
                
                for i in range(self.elastic_constants_points):
                    job = self.create_jobname(self.jobname + jobname_dist[i])
                    jobnames.append(job)
                    self.emto.set_values(sws=self.sws, jobname=job, latname=latname_dist[i],
                                         latpath=self.folder)
                    
                    common.check_folders(self.folder, self.folder + "/kgrn")
                    common.check_folders(self.folder + "/kgrn/tmp")
                    common.check_folders(self.folder + "/kfcd")
                    common.check_folders(self.folder + "/fit")
                    self.emto.kgrn.write_input_file(folder=self.folder)
                    self.emto.kfcd.write_input_file(folder=self.folder)
                    self.emto.batch.write_input_file(folder=self.folder)

            # These following lines are for the four-atom basis with
            # simple monoclinic lattice.
            """
            for i in range(len(self.elastic_constants_deltas)):
                self.lattice.distortion(lat='hcp', dist='mono', ca=self.ca, index=i,
                                        deltas=self.elastic_constants_deltas)

                self.lattice.set_values(
                    jobname=latname_dist[i], latpath=self.folder)
                self.lattice.bmdl.write_input_file(folder=self.folder)
                self.lattice.kstr.write_input_file(folder=self.folder)
                self.lattice.shape.write_input_file(folder=self.folder)
                self.lattice.batch.write_input_file(folder=self.folder)

            self.emto.set_values(ibz=12, nkx=30, nky=20, nkz=20)

            ################################################################
            # Atconf related arrays need to be modified because we now have
            # a four atom basis.
            ################################################################

            self.atoms = np.array([self.atoms, self.atoms]).flatten()
            self.concs = np.array([self.concs, self.concs]).flatten()

            self.iqs = np.zeros(len(self.atoms), dtype='int32')
            len_div = len(self.iqs) // 4
            for i in range(4):
                self.iqs[i * len_div:(i + 1) * len_div] = i + 1

            self.splts = np.array([self.splts, self.splts]).flatten()
            self.itas = np.array([self.itas, self.itas]).flatten()

            self.emto.set_values(atoms=self.atoms, iqs=self.iqs, itas=self.itas,
                                 concs=self.concs, splts=self.splts)

            for i in range(len(jobname_dist)):
                job = self.create_jobname(self.jobname + jobname_dist[i])
                jobnames.append(job)
                self.emto.set_values(sws=self.sws, jobname=job, latname=latname_dist[i],
                                     latpath=self.folder)

                common.check_folders(
                    self.folder, self.folder + "/kgrn", self.folder + "/kgrn/tmp")
                common.check_folders(self.folder + "/kfcd")
                common.check_folders(self.folder + "/fit")
                self.emto.kgrn.write_input_file(folder=self.folder)
                self.emto.kfcd.write_input_file(folder=self.folder)
                self.emto.batch.write_input_file(folder=self.folder)
            """

            return jobnames

    def elastic_constants_batch_calculate(
            self, sws=None, bmod=None, ca=None, R=None, cs=None):
        """Calculates the elastic constants of a system using the parallelism
        of the batch system.

        This is a combination of the batch_generate and batch_analyze functions.

        :param sws: WS-radius (Default value = None)
        :type sws: float
        :param bmod: Bulk modulus (Default value = None)
        :type bmod: float
        :param ca: hpc c/a ratio  (Default value = None)
        :type ca: float
        :param R: The dimensionless quantity of hcp systems (Default value = None)
        :type R: float
        :param cs: Second order c/a-derivative of the energy (Default value = None)
        :type cs: float
        :returns: None
        :rtype: None
        """

        import time

        # Mission critical parameters
        if sws is None and self.sws is None:
            sys.exit(
                'System.elastic_constants_batch_calculate(): \'sws\' has not been specified!')
        elif sws is not None:
            self.sws = sws

        if bmod is None and self.bmod is None:
            sys.exit(
                'System.elastic_constants_batch_calculate(): \'bmod\' has not been specified!')
        elif bmod is not None:
            self.bmod = bmod

        if ca is None and self.lat == 'hcp':
            sys.exit('System.elastic_constants_batch_calculate():' +
                     ' \'ca\' (c/a) has not been specified!')
        else:
            self.ec_batch_calculate_ca = ca

        if R is None and self.lat == 'hcp':
            sys.exit(
                'System.elastic_constants_batch_calculate(): \'R\' has not been specified!')
        else:
            self.ec_batch_calculate_R = R

        if cs is None and self.lat == 'hcp':
            sys.exit(
                'System.elastic_constants_batch_calculate(): \'cs\' has not been specified!')
        else:
            self.ec_batch_calculate_cs = cs

        # Generate input files
        if self.lat == 'bcc' or self.lat == 'fcc':
            jobnames = self.elastic_constants_batch_generate(sws=self.sws)

            # Submit calculation to the batch system
            jobIDs = self.submit_jobs(jobnames, folder=self.folder)

            # Wait until all the jobs have finished
            self.wait_for_jobs(jobIDs)

            # Now we can analyze the results
            self.elastic_constants_analyze(sws=self.sws, bmod=self.bmod)
            return

        elif self.lat == 'hcp':
            jobnames = self.elastic_constants_batch_generate(
                sws=self.sws, ca=self.ec_batch_calculate_ca)
            jobnames_lat = ['hcpo0_ca0',
                            'hcpo1_ca0',
                            'hcpo2_ca0',
                            'hcpo3_ca0',
                            'hcpo4_ca0',
                            'hcpo5_ca0',
                            'hcpm0_ca0',
                            'hcpm1_ca0',
                            'hcpm2_ca0',
                            'hcpm3_ca0',
                            'hcpm4_ca0',
                            'hcpm5_ca0']

            # First submit and run the lattice calculations
            jobIDs_lat = self.submit_jobs(jobnames_lat, folder=self.folder)

            # Wait until all the jobs have finished
            self.wait_for_jobs(jobIDs_lat)

            # Structure calculations have finished, now submit
            # KGRN and KFCD calculation.
            jobIDs = self.submit_jobs(jobnames, folder=self.folder)

            # Wait until all the jobs have finished
            self.wait_for_jobs(jobIDs)

            # Now we can analyze the results
            self.elastic_constants_analyze(
                sws=self.sws, bmod=self.bmod, ca=self.ec_batch_calculate_ca,
                R=self.ec_batch_calculate_R, cs=self.ec_batch_calculate_cs)
            return

    def elastic_constants_serial_calculate(self, sws=None, bmod=None, ca=None, R=None, cs=None):
        """Calculates elastic constants on one CPU without using the batch system.

        At the end the results are printed on screen.

        :param sws: WS-radius (Default value = None)
        :type sws: float
        :param bmod: Bulk modulus (Default value = None)
        :type bmod: float
        :param ca: hpc c/a ratio  (Default value = None)
        :type ca: float
        :returns: None
        :rtype: None
        """

        from pyemto.EOS.EOS import EOS

        eos = EOS(name=self.jobname, xc=self.xc, method='morse', units='bohr')

        # Mission critical parameters
        if sws is None and self.sws is None:
            sys.exit(
                'System.elastic_constants_serial_calculate(): \'sws\' has not been specified!')
        elif sws is not None:
            self.sws = sws

        if bmod is None and self.bmod is None:
            sys.exit('System.elastic_constants_serial_calculate():' +
                     ' \'bmod\' has not been specified!')
        elif bmod is not None:
            self.bmod = bmod

        if ca is None and self.lat == 'hcp':
            sys.exit('System.elastic_constants_serial_calculate():' +
                     ' \'ca\' (c/a) has not been specified!')
        elif ca is not None:
            self.ca = ca

        if R is None and self.lat == 'hcp':
            sys.exit(
                'System.elastic_constants_analyze(): \'R\' has not been specified!')
        else:
            self.ec_analyze_R = R

        if cs is None and self.lat == 'hcp':
            sys.exit(
                'System.elastic_constants_analyze(): \'cs\' has not been specified!')
        else:
            self.ec_analyze_cs = cs

        deltas = self.elastic_constants_deltas

        if self.lat == 'bcc' or self.lat == 'fcc':
            # Orthorhombic distortion for c' first

            if self.lat == 'bcc':
                jobname_dist = [
                    '_bcco0', '_bcco1', '_bcco2', '_bcco3', '_bcco4', '_bcco5']
                latname_dist = [
                    'bcco0', 'bcco1', 'bcco2', 'bcco3', 'bcco4', 'bcco5']
                self.emto.set_values(ibz=10, nkx=27, nky=27, nkz=27)
            elif self.lat == 'fcc':
                jobname_dist = [
                    '_fcco0', '_fcco1', '_fcco2', '_fcco3', '_fcco4', '_fcco5']
                latname_dist = [
                    'fcco0', 'fcco1', 'fcco2', 'fcco3', 'fcco4', 'fcco5']
                self.emto.set_values(ibz=11, nkx=27, nky=27, nkz=27)

            en_cprime = []

            for i in range(len(jobname_dist)):
                already = False
                job = self.create_jobname(self.jobname + jobname_dist[i])
                self.emto.set_values(
                    sws=self.sws, jobname=job, latname=latname_dist[i])

                # check if calculations are already done
                already = self.check_conv(job, folder=self.folder)

                if all(already):
                    conv = (True, True)
                else:
                    if already[0] and not already[1]:
                        conv = self.runemto(
                            jobname=job, folder=self.folder, onlyKFCD=True)
                    else:
                        conv = self.runemto(jobname=job, folder=self.folder)

                # KGRN has crashed, find out why
                if conv[0] == False:
                    self.which_error(job, folder=self.folder)
                    quit()

                en = self.get_energy(job, folder=self.folder, func=self.xc)
                en_cprime.append(en)

            # Convert energies into a numpy array
            en_cprime = np.asarray(en_cprime)

            # Next calculate the monoclinic distortion for c44

            if self.lat == 'bcc':
                jobname_dist = [
                    '_bccm0', '_bccm1', '_bccm2', '_bccm3', '_bccm4', '_bccm5']
                latname_dist = [
                    'bccm0', 'bccm1', 'bccm2', 'bccm3', 'bccm4', 'bccm5']
                self.emto.set_values(ibz=11, nkx=27, nky=27, nkz=37)
            elif self.lat == 'fcc':
                jobname_dist = [
                    '_fccm0', '_fccm1', '_fccm2', '_fccm3', '_fccm4', '_fccm5']
                latname_dist = [
                    'fccm0', 'fccm1', 'fccm2', 'fccm3', 'fccm4', 'fccm5']
                self.emto.set_values(ibz=10, nkx=27, nky=37, nkz=27)

            en_c44 = []

            for i in range(len(jobname_dist)):
                already = False
                job = self.create_jobname(self.jobname + jobname_dist[i])
                self.emto.set_values(
                    sws=self.sws, jobname=job, latname=latname_dist[i])

                # check if calculations are already done
                already = self.check_conv(job, folder=self.folder)

                if all(already):
                    conv = (True, True)
                else:
                    if already[0] and not already[1]:
                        conv = self.runemto(
                            jobname=job, folder=self.folder, onlyKFCD=True)
                    else:
                        conv = self.runemto(jobname=job, folder=self.folder)

                # KGRN has crashed, find out why
                if conv[0] == False:
                    self.which_error(job, folder=self.folder)
                    quit()

                en = self.get_energy(job, folder=self.folder, func=self.xc)
                en_c44.append(en)

            # Convert energies into a numpy array
            en_c44 = np.asarray(en_c44)

            # All calculations have been done, now it's time to fit the results
            popt_cprime, cprime_rsq = eos.distortion_fit(deltas, en_cprime)
            popt_c44, c44_rsq = eos.distortion_fit(deltas, en_c44)

            volume = 4.0 / 3.0 * np.pi * self.sws**3

            cprime = popt_cprime[0] / 2.0 / volume * self.RyBohr3_to_GPa
            c44 = popt_c44[0] / 2.0 / volume * self.RyBohr3_to_GPa

            c11 = self.bmod + 4.0 / 3.0 * cprime
            c12 = self.bmod - 2.0 / 3.0 * cprime

            # Polycrystalline elastic constants
            #
            # B = bulk modulus
            # G = shear modulus
            # E = Young modulus
            # v = Poisson ratio

            # Voigt average
            BV = (c11 + 2 * c12) / 3.0
            #BV = self.bmod
            GV = (c11 - c12 + 3 * c44) / 5.0
            EV = 9 * BV * GV / (3 * BV + GV)
            vV = (3 * BV - 2 * GV) / (6 * BV + 2 * GV)

            # Reuss average
            BR = BV
            #BR = self.bmod
            GR = 5 * (c11 - c12) * c44 / (4 * c44 + 3 * (c11 - c12))
            ER = 9 * BR * GR / (3 * BR + GR)
            vR = (3 * BR - 2 * GR) / (6 * BR + 2 * GR)

            # Hill average
            BH = (BV + BR) / 2.0
            #BH = self.bmod
            GH = (GV + GR) / 2.0
            EH = 9 * BH * GH / (3 * BH + GH)
            vH = (3 * BH - 2 * GH) / (6 * BH + 2 * GH)

            # Elastic anisotropy
            AVR = (GV - GR) / (GV + GR)

            print("")
            print('***cubic_elastic_constants***')
            print("")
            print(self.jobname)
            print("")
            print('c11(GPa) = {0:6.2f}'.format(c11))
            print('c12(GPa) = {0:6.2f}'.format(c12))
            print(
                'c44(GPa) = {0:6.2f}, R-squared = {1:8.6f}'.format(c44, c44_rsq))
            print(
                'c\' (GPa) = {0:6.2f}, R-squared = {1:8.6f}'.format(cprime, cprime_rsq))
            print('B  (GPa) = {0:6.2f}'.format(self.bmod))
            print("")
            print('Voigt average:')
            print("")
            print('BV(GPa)  = {0:6.2f}'.format(BV))
            print('GV(GPa)  = {0:6.2f}'.format(GV))
            print('EV(GPa)  = {0:6.2f}'.format(EV))
            print('vV(GPa)  = {0:6.2f}'.format(vV))
            print("")
            print('Reuss average:')
            print("")
            print('BR(GPa)  = {0:6.2f}'.format(BR))
            print('GR(GPa)  = {0:6.2f}'.format(GR))
            print('ER(GPa)  = {0:6.2f}'.format(ER))
            print('vR(GPa)  = {0:6.2f}'.format(vR))
            print("")
            print('Hill average:')
            print("")
            print('BH(GPa)  = {0:6.2f}'.format(BH))
            print('GH(GPa)  = {0:6.2f}'.format(GH))
            print('EH(GPa)  = {0:6.2f}'.format(EH))
            print('vH(GPa)  = {0:6.2f}'.format(vH))
            print("")
            print('Elastic anisotropy:')
            print("")
            print('AVR(GPa)  = {0:6.2f}'.format(AVR))

            return

        elif self.lat == 'hcp':

            # Check whether Two-center Taylor expansion is on/off
            if self.emto.kgrn.expan == 'M':
                kappaw = self.kappaw_default
                self.lattice.set_values(kappaw=kappaw)

            # Check whether we need to create some folders
            # for the structure output files.
            common.check_folders(self.folder)
            common.check_folders(self.folder + '/bmdl')
            common.check_folders(self.folder + '/kstr')
            common.check_folders(self.folder + '/shape')

            # Orthorhombic distortion for c66 first
            jobname_dist = ['_hcpo0_ca0', '_hcpo1_ca0',
                            '_hcpo2_ca0', '_hcpo3_ca0', '_hcpo4_ca0', '_hcpo5_ca0']
            latname_dist = ['hcpo0_ca0', 'hcpo1_ca0',
                            'hcpo2_ca0', 'hcpo3_ca0', 'hcpo4_ca0', 'hcpo5_ca0']
            self.emto.set_values(ibz=9, nkx=31, nky=19, nkz=19)

            en_c66 = []

            for i in range(len(jobname_dist)):
                # With hcp the structure depends on the c/a ratio. Therefore we also have
                # to generate the corresponding structure files.

                self.lattice.distortion(lat='hcp', dist='ortho', ca=self.ca, index=i,
                                        deltas=self.elastic_constants_deltas)

                self.lattice.set_values(
                    jobname=latname_dist[i], latpath=self.folder)
                self.runlattice(jobname=latname_dist[i], folder=self.folder)

                already = False
                job = self.create_jobname(self.jobname + jobname_dist[i])
                self.emto.set_values(
                    sws=self.sws, jobname=job, latname=latname_dist[i])

                # check if calculations are already done
                already = self.check_conv(job, folder=self.folder)

                if all(already):
                    conv = (True, True)
                else:
                    if already[0] and not already[1]:
                        conv = self.runemto(
                            jobname=job, folder=self.folder, onlyKFCD=True)
                    else:
                        conv = self.runemto(jobname=job, folder=self.folder)

                # KGRN has crashed, find out why
                if conv[0] == False:
                    self.which_error(job, folder=self.folder)
                    quit()

                en = self.get_energy(job, folder=self.folder, func=self.xc)
                en_c66.append(en)

            # Convert energies into a numpy array
            en_c66 = np.asarray(en_cprime)

            # Monoclinic distortion for c44 next
            jobname_dist = ['_hcpm0_ca0',
                            '_hcpm1_ca0',
                            '_hcpm2_ca0',
                            '_hcpm3_ca0',
                            '_hcpm4_ca0',
                            '_hcpm5_ca0']
            latname_dist = ['hcpm0_ca0',
                            'hcpm1_ca0',
                            'hcpm2_ca0',
                            'hcpm3_ca0',
                            'hcpm4_ca0',
                            'hcpm5_ca0']
            self.emto.set_values(ibz=12, nkx=30, nky=20, nkz=20)

            en_c44 = []

            for i in range(len(jobname_dist)):
                # With hcp the structure depends on the c/a ratio. Therefore we also have
                # to generate the corresponding structure files.

                self.lattice.distortion(lat='hcp', dist='mono', ca=self.ca, index=i,
                                        deltas=self.elastic_constants_deltas)

                self.lattice.set_values(
                    jobname=latname_dist[i], latpath=self.folder)
                self.runlattice(jobname=latname_dist[i], folder=self.folder)

                ###############################################################
                # Atconf related arrays need to be modified because we now have
                # a four atom basis.
                ###############################################################

                self.atoms = np.array([self.atoms, self.atoms]).flatten()
                self.concs = np.array([self.concs, self.concs]).flatten()

                self.iqs = np.zeros(len(self.atoms), dtype='int32')
                len_div = len(self.iqs) // 4
                for i in range(4):
                    self.iqs[i * len_div:(i + 1) * len_div] = i + 1

                self.splts = np.array([self.splts, self.splts]).flatten()
                self.itas = np.array([self.itas, self.itas]).flatten()

                self.emto.set_values(atoms=self.atoms, iqs=self.iqs, itas=self.itas,
                                     concs=self.concs, splts=self.splts)

                already = False
                job = self.create_jobname(self.jobname + jobname_dist[i])
                self.emto.set_values(
                    sws=self.sws, jobname=job, latname=latname_dist[i])

                # check if calculations are already done
                already = self.check_conv(job, folder=self.folder)

                if all(already):
                    conv = (True, True)
                else:
                    if already[0] and not already[1]:
                        conv = self.runemto(
                            jobname=job, folder=self.folder, onlyKFCD=True)
                    else:
                        conv = self.runemto(jobname=job, folder=self.folder)

                # KGRN has crashed, find out why
                if conv[0] == False:
                    self.which_error(job, folder=self.folder)
                    quit()

                en = self.get_energy(job, folder=self.folder, func=self.xc)
                en_c44.append(en)

            # Convert energies into a numpy array
            en_c44 = np.asarray(en_c44)

            # All calculations have been done, now it's time to fit the results
            popt_c66, c66_rsq = eos.distortion_fit(deltas, en_c66)
            popt_c44, c44_rsq = eos.distortion_fit(deltas, en_c44)

            volume = 4.0 / 3.0 * np.pi * self.sws**3

            c66 = popt_c66[0] / 2.0 / volume * self.RyBohr3_to_GPa
            c44 = popt_c44[0] / 2.0 / volume * self.RyBohr3_to_GPa

            c11 = self.bmod + c66 + self.ec_analyze_cs * \
                (2 * self.ec_analyze_R - 1)**2 / 18.0
            c12 = self.bmod - c66 + self.ec_analyze_cs * \
                (2 * self.ec_analyze_R - 1)**2 / 18.0
            c13 = self.bmod + 1.0 / 9.0 * self.ec_analyze_cs * (
                2 * self.ec_analyze_R**2 + self.ec_analyze_R - 1)
            c33 = self.bmod + 2.0 / 9.0 * \
                self.ec_analyze_cs * (self.ec_analyze_R + 1)**2
            c2 = c33 * (c11 + c12) - 2.0 * c13**2

            # Polycrystalline elastic constants
            #
            # B = bulk modulus
            # G = shear modulus
            # E = Young modulus
            # v = Poisson ratio

            # Voigt average
            BV = (2 * c11 + 2 * c12 + 4 * c13 + c33) / 9.0
            GV = (12 * c44 + 12 * c66 + self.ec_analyze_cs) / 30.0
            EV = 9 * BV * GV / (3 * BV + GV)
            vV = (3 * BV - 2 * GV) / (6 * BV + 2 * GV)

            # Reuss average
            BR = self.bmod
            GR = 5.0 / 2.0 * (c44 * c66 * c2) / \
                ((c44 + c66) * c2 + 3.0 * BV * c44 * c66)
            ER = 9 * BR * GR / (3 * BR + GR)
            vR = (3 * BR - 2 * GR) / (6 * BR + 2 * GR)

            # Hill average
            BH = (BV + BR) / 2.0
            #BH = self.bmod
            GH = (GV + GR) / 2.0
            EH = 9 * BH * GH / (3 * BH + GH)
            vH = (3 * BH - 2 * GH) / (6 * BH + 2 * GH)

            # Elastic anisotropy
            AVR = (GV - GR) / (GV + GR)

            print("")
            print('***hcp_elastic_constants***')
            print("")
            print(self.jobname)
            print("")
            print('c11(GPa) = {0:6.2f}'.format(c11))
            print('c12(GPa) = {0:6.2f}'.format(c12))
            print('c13(GPa) = {0:6.2f}'.format(c13))
            print('c33(GPa) = {0:6.2f}'.format(c33))
            print(
                'c44(GPa) = {0:6.2f}, R-squared = {1:8.6f}'.format(c44, c44_rsq))
            print(
                'c66(GPa) = {0:6.2f}, R-squared = {1:8.6f}'.format(c66, c66_rsq))
            print('B  (GPa) = {0:6.2f}'.format(self.bmod))
            print("")
            print('Voigt average:')
            print("")
            print('BV(GPa)  = {0:6.2f}'.format(BV))
            print('GV(GPa)  = {0:6.2f}'.format(GV))
            print('EV(GPa)  = {0:6.2f}'.format(EV))
            print('vV(GPa)  = {0:6.2f}'.format(vV))
            print("")
            print('Reuss average:')
            print("")
            print('BR(GPa)  = {0:6.2f}'.format(BR))
            print('GR(GPa)  = {0:6.2f}'.format(GR))
            print('ER(GPa)  = {0:6.2f}'.format(ER))
            print('vR(GPa)  = {0:6.2f}'.format(vR))
            print("")
            print('Hill average:')
            print("")
            print('BH(GPa)  = {0:6.2f}'.format(BH))
            print('GH(GPa)  = {0:6.2f}'.format(GH))
            print('EH(GPa)  = {0:6.2f}'.format(EH))
            print('vH(GPa)  = {0:6.2f}'.format(vH))
            print("")
            print('Elastic anisotropy:')
            print("")
            print('AVR(GPa)  = {0:6.2f}'.format(AVR))

            return

    ##########################################################################
    #                                                                        #
    # Internal routines start here                                           #
    #                                                                        #
    ##########################################################################

    def find_lc(self, delta=0.01, prn=True, xc='PBE'):
        """Computes initial estimates for the ground state quantities for cubic systems.

        :param delta: Step size for the volume vs. energy array (Default value = 0.01)
        :type delta: float
        :param prn: True if results should be printed, False if not (Default value = True)
        :type prn: boolean
        :param xc: Choice of the xc-functional (Default value = 'PBE')
        :type xc: str
        :returns: WS-radius, bulk modulus, c/a and energy
        :rtype: float, float, float, float
        """

        from pyemto.EOS.EOS import EOS

        eos = EOS(name=self.jobname, xc=xc, method='morse', units='bohr')

        if prn:
            print('')
            print('*****find_lc*****')
            print('')
        #enough = False
        energies = []
        swses = []

        # Compute first two points
        # Start at the initial sws
        already = False
        self.sws = self.lc_initial_sws
        job = self.create_jobname(self.jobname)
        self.emto.set_values(sws=self.sws, jobname=job)

        # check if calculations are already done
        already = self.check_conv(job, folder=self.folder)

        if self.lc_skip or (all(already) and not self.lc_rerun):
            conv = (True, True)
        else:
            if already[0] and not already[1]:
                conv = self.runemto(
                    jobname=job, folder=self.folder, onlyKFCD=True)
            else:
                conv = self.runemto(jobname=job, folder=self.folder)

        # KGRN has crashed, find out why
        if conv[0] == False:
            self.which_error(job, folder=self.folder)
            quit()

        en = self.get_energy(job, folder=self.folder, func=xc)

        energies.append(en)
        swses.append(self.sws)

        # Compute second point at sws+delta
        already = False
        self.sws = self.lc_initial_sws + delta
        job = self.create_jobname(self.jobname)
        self.emto.set_values(sws=self.sws, jobname=job)
        # check if calculations are already done
        already = self.check_conv(job, folder=self.folder)

        #print('already = ',already)

        if self.lc_skip or (all(already) and not self.lc_rerun):
            conv = (True, True)
        else:
            if already[0] and not already[1]:
                conv = self.runemto(
                    jobname=job, folder=self.folder, onlyKFCD=True)
            else:
                conv = self.runemto(jobname=job, folder=self.folder)

        #print('conv = ',conv)

        # KGRN has crashed, find out why
        if conv[0] == False:
            self.which_error(job, folder=self.folder)
            quit()

        en = self.get_energy(job, folder=self.folder, func=xc)
        energies.append(en)
        swses.append(self.sws)

        # The initial 2 points are now ready and
        # we use them to predict where to calculate next point

        next_sws, minfound = self.predict_next_sws(swses, energies, delta)

        # 2 points is not enought to calculate equation of state
        enough = False

        # Loop until we have enough points to calculate initial sws
        iteration = 0
        while not enough:
            iteration += 1

            # if 25 is not enough one should check initial sws
            if iteration > 25:
                print(
                    "SWS loop did not converge in {0} iterations!".format(iteration))
                quit()

            # Calculate next point
            self.sws = next_sws
            job = self.create_jobname(self.jobname)
            self.emto.set_values(sws=self.sws, jobname=job)

            # First check if calculations are already done
            already = False
            already = self.check_conv(job, folder=self.folder)
            # Use existing calculations if available.
            if self.lc_skip or (all(already) and not self.lc_rerun):
                conv = (True, True)
            else:
                if already[0] and not already[1]:
                    conv = self.runemto(
                        jobname=job, folder=self.folder, onlyKFCD=True)
                else:
                    conv = self.runemto(jobname=job, folder=self.folder)

            # KGRN has crashed, find out why
            if conv[0] == False:
                self.which_error(job, folder=self.folder)
                quit()

            en = self.get_energy(job, folder=self.folder, func=xc)
            energies.append(en)
            swses.append(self.sws)

            # Check do we have minumun and predict next sws
            next_sws, minfound = self.predict_next_sws(swses, energies, delta)

            # Check if we have mimimum and enough points
            # 9 points should be enough for EOS
            if minfound and len(swses) > 9:
                enough = True

        #print('debug find_lc: ',swses)

        if prn:
            self.print_sws_ens('find_lc', swses, energies)

        sws0, e0, B, grun = eos.fit(swses, energies)

        # These functions create files on disk about the data to be fitted
        # as well as the results of the fit.
        # eos.prepareData()
        #sws0,e0,B,grun = eos.fit2file()

        return sws0, B, e0

    def refine_lc(self, sws, delta=0.01, prn=True, xc='PBE'):
        """Calculates a more accurate equilibrium volume for cubic systems.

        :param sws: WS-radius
        :type sws: float
        :param delta: Step size for the volume vs. energy array (Default value = 0.01)
        :type delta: float
        :param prn: True if results should be printed, False if not (Default value = True)
        :type prn: boolean
        :param xc: Choice of the xc-functional (Default value = 'PBE')
        :type xc: str
        :returns: Ws-radius, bulk modulus, c/a and energy
        :rtype: float, float, float, float
        """

        from pyemto.EOS.EOS import EOS

        eos = EOS(name=self.jobname, xc=xc, method='morse', units='bohr')

        if prn:
            print('')
            print('*****refine_lc*****')
            print('')

        # make sws ranges around given sws
        points = []
        for i in range(-3, 7):
            points.append(i * delta)

        energies = []
        swses = []
        for p in points:
            already = False
            next_sws = sws + p
            # Make inputs and job name
            self.sws = next_sws
            job = self.create_jobname(self.jobname)
            self.emto.set_values(sws=self.sws, jobname=job)
            # check if calculations are already done
            already = self.check_conv(job, folder=self.folder)
            if self.lc_skip or (all(already) and not self.lc_rerun):
                conv = (True, True)
            else:
                if already[0] and not already[1]:
                    conv = self.runemto(
                        jobname=job, folder=self.folder, onlyKFCD=True)
                else:
                    conv = self.runemto(jobname=job, folder=self.folder)

            # KGRN has crashed, find out why
            if conv[0] == False:
                self.which_error(job, folder=self.folder)
                quit()

            en = self.get_energy(job, folder=self.folder, func=xc)
            energies.append(en)
            swses.append(self.sws)

        if prn:
            self.print_sws_ens('refine_lc', swses, energies)

        sws0, e0, B, grun = eos.fit(swses, energies)

        # These functions create files on disk about the data to be fitted
        # as well as the results of the fit.
        # eos.prepareData()
        #sws0,e0,B,grun = eos.fit2file()

        return sws0, B, e0

    def find_lc_hcp(self, delta=0.01, prn=True, xc='PBE'):
        """Computes initial estimates for the ground state quantities for hcp systems.

        :param delta: Step size for the volume vs. energy array (Default value = 0.01)
        :type delta: float
        :param prn: True if results should be printed, False if not (Default value = True)
        :type prn: boolean
        :param xc: Choice of the xc-functional (Default value = 'PBE')
        :type xc: str
        :returns: WS-radius, bulk modulus, c/a and energy
        :rtype: float, float, float, float
        """

        from pyemto.EOS.EOS import EOS

        eos_hcp = EOS(name=self.jobname, xc=xc, method='morse', units='bohr')

        if prn:
            print('')
            print('*****find_lc_hcp*****')
            print('')
        #enough = False

        # Before anything else happens we generate the
        # structure output files (BMDL, KSTR and SHAPE) for a reasonable c over
        # a mesh
        ca_mesh = self.ca_range_default
        ca_mesh_len = len(ca_mesh)
        # sws-optimization is done using this particular c/a value,
        ca_sws_ind = 2
        # all the other c/a's will use these same volumes.
        # Fit an n'th order polynomial to the energy vs. c/a data.
        ca_fit_order = 2
        ca_latpaths = []
        ca_prefixes = []

        # A 2D-array [i,j], where i = c/a axis and j = sws axis
        energies = np.zeros((ca_mesh_len, 50))
        energies0 = []  # c/a optimized energy for a given WS-radius
        swses = []     # List of WS-radii
        cas0 = []      # Energetically optimized c/a's for a given WS-radius

        for i in range(0, ca_mesh_len):
            ca_prefixes.append("/ca_{0:8.6f}".format(ca_mesh[i]))
            ca_latpaths.append(self.latpath + ca_prefixes[i])

        # Check whether structure files already exists
        # and if not run the calculations.
        for i in range(0, ca_mesh_len):
            self.lattice.set_values(ca=ca_mesh[i])
            # self.lattice.bmdl.write_input_file(folder=ca_latpaths[i])
            str_exists = self.check_str(self.latname, folder=ca_latpaths[i])
            if str_exists == False:
                self.runlattice(jobname=self.latname, folder=ca_latpaths[i])

        # Compute first two points
        # Start at the initial sws

        already = False
        self.sws = self.lc_initial_sws
        job = self.create_jobname(self.jobname)
        self.emto.set_values(sws=self.sws, jobname=job)

        # check if calculations are already done
        for i in range(ca_mesh_len):
            # We have to remember to update the lattice path every
            # time we change c/a
            self.emto.set_values(latpath=ca_latpaths[i])
            hcp_subfolder = self.folder + "/{0}".format(ca_prefixes[i])
            already = self.check_conv(job, folder=hcp_subfolder)

            if self.lc_skip or (all(already) and not self.lc_rerun):
                conv = (True, True)
            else:
                if already[0] and not already[1]:
                    conv = self.runemto(
                        jobname=job, folder=hcp_subfolder, onlyKFCD=True)
                else:
                    conv = self.runemto(jobname=job, folder=hcp_subfolder)

            # KGRN has crashed, find out why
            if conv[0] == False:
                self.which_error(job, folder=hcp_subfolder)
                quit()

            en = self.get_energy(job, folder=hcp_subfolder, func=xc)
            energies[i, 0] = en

        swses.append(self.sws)

        # Calculate energetically optimized c/a and the corresponding energy
        # print(energies[:,:15])
        ca0, en0 = eos_hcp.ca_fit(ca_mesh, energies[:, 0], ca_fit_order)
        cas0.append(ca0)
        energies0.append(en0)

        # Compute second point at initial sws + delta
        already = False
        self.sws = self.lc_initial_sws + delta
        job = self.create_jobname(self.jobname)
        self.emto.set_values(sws=self.sws, jobname=job)

        # check if calculations are already done
        for i in range(ca_mesh_len):
            # We have to remember to update the lattice path every
            # time we change c/a
            self.emto.set_values(latpath=ca_latpaths[i])
            hcp_subfolder = self.folder + "/{0}".format(ca_prefixes[i])
            already = self.check_conv(job, folder=hcp_subfolder)

            if self.lc_skip or (all(already) and not self.lc_rerun):
                conv = (True, True)
            else:
                if already[0] and not already[1]:
                    conv = self.runemto(
                        jobname=job, folder=hcp_subfolder, onlyKFCD=True)
                else:
                    conv = self.runemto(jobname=job, folder=hcp_subfolder)

            # KGRN has crashed, find out why
            if conv[0] == False:
                self.which_error(job, folder=hcp_subfolder)
                quit()

            en = self.get_energy(job, folder=hcp_subfolder, func=xc)
            energies[i, 1] = en

        swses.append(self.sws)

        # Calculate energetically optimized c/a and the corresponding energy
        # print(energies[:,:15])
        ca0, en0 = eos_hcp.ca_fit(ca_mesh, energies[:, 1], ca_fit_order)
        cas0.append(ca0)
        energies0.append(en0)

        # The initial 2 points are now ready and
        # we use them to predict where to calculate next point

        next_sws, minfound = self.predict_next_sws(swses, energies0, delta)

        # 2 points is not enought to calculate equation of state
        enough = False

        # Loop until we have enough points to calculate initial sws
        iteration = 0
        while not enough:
            iteration += 1

            # if 25 is not enough one should check initial sws
            if iteration > 25:
                print(
                    "SWS loop did not converge in {0} iterations!".format(iteration))
                quit()

            # Calculate next point
            already = False
            self.sws = next_sws
            job = self.create_jobname(self.jobname)
            self.emto.set_values(sws=self.sws, jobname=job)

            # First check if calculations are already done
            for i in range(ca_mesh_len):
                # We have to remember to update the lattice path every
                # time we change c/a
                self.emto.set_values(latpath=ca_latpaths[i])
                hcp_subfolder = self.folder + "/{0}".format(ca_prefixes[i])
                already = self.check_conv(job, folder=hcp_subfolder)

                if self.lc_skip or (all(already) and not self.lc_rerun):
                    conv = (True, True)
                else:
                    if already[0] and not already[1]:
                        conv = self.runemto(
                            jobname=job, folder=hcp_subfolder, onlyKFCD=True)
                    else:
                        conv = self.runemto(jobname=job, folder=hcp_subfolder)

                # KGRN has crashed, find out why
                if conv[0] == False:
                    self.which_error(job, folder=hcp_subfolder)
                    quit()

                en = self.get_energy(job, folder=hcp_subfolder, func=xc)
                energies[i, iteration + 1] = en

            swses.append(self.sws)

            # Calculate energetically optimized c/a and the corresponding energy
            # print(energies[:,:15])
            ca0, en0 = eos_hcp.ca_fit(
                ca_mesh, energies[:, iteration + 1], ca_fit_order)
            cas0.append(ca0)
            energies0.append(en0)

            # Check do we have minumun and predict next sws
            next_sws, minfound = self.predict_next_sws(swses, energies0, delta)

            # Check if we have mimimum and enough points
            # 9 points should be enough for EOS
            if minfound and len(swses) > 9:
                enough = True

        #print('debug find_lc_hcp: ',swses)

        if prn:
            self.print_sws_ens_hcp('find_lc_hcp', swses, energies0, cas0)

        sws0, e0, B0, grun = eos_hcp.fit(swses, energies0)

        # These functions create files on disk about the data to be fitted
        # as well as the results of the fit.
        # eos.prepareData()
        #sws0,e0,B0,grun = eos.fit2file()

        # Now that we have the ground state WS-radius sws0 we can use the cas0 array
        # to compute the corresponding ground state ca0
        ca0_vs_sws = np.polyfit(swses, cas0, 2)
        ca0_vs_sws = np.poly1d(ca0_vs_sws)
        c_over_a0 = ca0_vs_sws(sws0)

        return sws0, c_over_a0, B0, e0

    def refine_lc_hcp(self, sws, ca0, delta=0.01, prn=True, xc='PBE'):
        """Calculates a more accurate equilibrium volume for hcp systems.

        :param sws: WS-radius
        :type sws: float
        :param ca0: Previously computed eq. c/a ratio
        :type ca0: float
        :param delta: Step size for the volume vs. energy array (Default value = 0.01)
        :type delta: float
        :param prn: True if results should be printed, False if not (Default value = True)
        :type prn: boolean
        :param xc: Choice of the xc-functional (Default value = 'PBE')
        :type xc: str
        :returns: Ws-radius, bulk modulus, c/a and energy
        :rtype: float, float, float, float
        """

        from pyemto.EOS.EOS import EOS

        eos_hcp = EOS(name=self.jobname, xc=xc, method='morse', units='bohr')

        if prn:
            print('')
            print('*****refine_lc_hcp*****')
            print('')

        # First compute the structure with the optimized
        # c/a (found by the find_lc_hcp() function).
        self.lattice.set_values(ca=ca0)

        ca_prefix = "/ca_{0:8.6f}".format(ca0)
        ca_latpath = self.latpath + ca_prefix

        str_exists = self.check_str(self.latname, folder=ca_latpath)
        if str_exists == False:
            self.runlattice(jobname=self.latname, folder=ca_latpath)

        # make sws ranges around given sws
        points = []
        for i in range(-3, 7):
            points.append(i * delta)

        energies = []
        swses = []

        self.emto.set_values(latpath=ca_latpath)
        hcp_subfolder = self.folder + "/{0}".format(ca_prefix)

        for p in points:
            already = False
            next_sws = sws + p
            # Make inputs and job name
            self.sws = next_sws
            job = self.create_jobname(self.jobname)
            self.emto.set_values(sws=self.sws, jobname=job)
            # check if calculations are already done
            already = self.check_conv(job, folder=hcp_subfolder)
            if self.lc_skip or (all(already) and not self.lc_rerun):
                conv = (True, True)
            else:
                if already[0] and not already[1]:
                    conv = self.runemto(
                        jobname=job, folder=hcp_subfolder, onlyKFCD=True)
                else:
                    conv = self.runemto(jobname=job, folder=hcp_subfolder)

            # KGRN has crashed, find out why
            if conv[0] == False:
                self.which_error(job, folder=hcp_subfolder)
                quit()

            en = self.get_energy(job, folder=hcp_subfolder, func=xc)
            energies.append(en)
            swses.append(self.sws)

        if prn:
            self.print_sws_ens('refine_lc_hcp', swses, energies)

        sws0, e0, B, grun = eos_hcp.fit(swses, energies)

        c_over_a = ca0

        # These functions create files on disk about the data to be fitted
        # as well as the results of the fit.
        # eos.prepareData()
        #sws0,e0,B,grun = eos.fit2file()

        return sws0, B, c_over_a, e0

    def predict_next_sws(self, swses, en, maxdelta=0.05):
        """Predict next WS-radius based on a simple gradient descent algorithm.

        :param swses: List of current WS-radii
        :type swses: list(float)
        :param en: List of current energies
        :type en: list(float)
        :param maxdelta: Maximum step size (Default value = 0.05)
        :type maxdelta: float
        :returns: Next WS-radii and True if energy minimum has been found,
                  False if not yet found
        :rtype: float, boolean
        """

        # Check if we do have minimum here and predict directon for next point
        m = 0  # The first one is minimum at start
        minsws = 10000.0
        maxsws = 0.0
        for i in range(len(swses)):
            # Check if we have minimum volume
            if swses[i] < minsws:
                minsws = swses[i]
            elif swses[i] > maxsws:
                maxsws = swses[i]
            if en[i] < en[m]:
                m = i

        wehaveminumum = False

        # Possible cases
        if swses[m] == min(swses):  # No minimum decrease sws
            # One should make delta depend on energy slope at some point.
            delta = -maxdelta
            newsws = swses[m] + delta
        elif swses[m] == max(swses):  # No minimum increase sws
            delta = maxdelta
            newsws = swses[m] + delta
        else:
            wehaveminumum = True
            # Minimum is inside the sws range. Decide where to add new point
            larger = [s for s in swses if swses[m] < s]
            smaller = [s for s in swses if swses[m] > s]

            if len(larger) > len(smaller):
                # Decrease volume
                sws = minsws
                delta = -maxdelta
                newsws = sws + delta
            else:
                # Increase volume
                sws = maxsws
                delta = maxdelta
                newsws = sws + delta

        return newsws, wehaveminumum

    def print_sws_ens(self, string, swses, energies):
        """Prints the WS-radii and calculated energies of cubic systems

        :param string: Header for the printout
        :type string: str
        :param swses: List of WS-radii
        :type swses: list(float)
        :param energies: List of energies
        :type energies: list(float)
        :returns: None
        :rtype: None
        """

        str_len = len(string)
        print('' + '\n')
        print('************')
        print(str_len * '*')
        print(string)
        print(str_len * '*')
        print('************')
        print('  SWS        Energy')
        for i in range(len(swses)):
            print('{0:8.6f}  {1:12.6f}'.format(swses[i], energies[i]))
        print('' + '\n')
        return

    def print_sws_ens_hcp(self, string, swses, energies, cas):
        """Prints the WS-radii, ca/a ratios and calculated energies of hcp systems

        :param string: Header for the printout
        :type string: str
        :param swses: List of WS-radii
        :type swses: list(float)
        :param energies: List of energies
        :type energies: list(float)
        :param cas: List of c/a ratios
        :type cas: list(float)
        :returns: None
        :rtype: None
        """

        str_len = len(string)
        print('' + '\n')
        print(str_len * '*')
        print(string)
        print(str_len * '*')
        print('  SWS        Energy0        c/a0')
        for i in range(len(swses)):
            print('{0:8.6f}  {1:12.6f}  {2:8.6f}'.format(
                swses[i], energies[i], cas[i]))
        print('' + '\n')
        return

    def check_conv(self, jobname, folder="./"):
        """Checks the convergence of given KGRN and KFCD calculations by reading
        their output files.

        :param jobname: Name of the output files
        :type jobname: str
        :param folder: Name of the folder where the output files are located (Default value = "./")
        :type folder: str
        :returns: Convergence of KGRN (True/False), convergence of KFCD  (True/False)
        :rtype: boolean, boolean
        """

        folderKGRN = folder + '/kgrn/'
        # Check if we got convergence in KGRN
        prntfile = jobname + ".prn"

        convergedKGRN = False
        fn = os.path.join(folderKGRN, prntfile)
        if not os.path.isfile(fn):
            pass
        else:
            pfile = open(fn, "r")
            for line in pfile:
                if "Converged in" in line:
                    convergedKGRN = True
                    break
            pfile.close()

        folderKFCD = folder + '/kfcd/'
        # Check if we got convergence in KFCD
        prntfile = jobname + ".prn"

        convergedKFCD = False
        fn = os.path.join(folderKFCD, prntfile)
        if not os.path.isfile(fn):
            pass
        else:
            pfile = open(fn, "r")
            for line in pfile:
                if "Finished at:" in line:
                    convergedKFCD = True
                    break
            pfile.close()

        return convergedKGRN, convergedKFCD

    def runlattice(self, jobname=None, folder="./", EMTODIR=None):
        """Run BMDL, KSTR and SHAPE calculation **WITHOUT** using the batch system.

        :param jobname: Name of the input files (Default value = None)
        :type jobname: str
        :param folder: Name of the folder where the input files are located (Default value = "./")
        :type folder: str
        :param EMTODIR: Path to the EMTO installation folder (Default value = None)
        :type EMTODIR: str
        :returns: None
        :rtype: None
        """

        if EMTODIR is None:
            EMTODIR = self.EMTOdir

        if jobname is None:
            sys.exit("System.runlattice: \'jobname\' has to be given!")
        # Make sure folders exist
        common.check_folders(folder, folder + "/bmdl/")
        # Create input file
        self.lattice.bmdl.write_input_file(folder=folder)
        # Run BMDL
        job = os.path.join(folder, jobname)
        command = "cd {0}; ".format(
            folder) + EMTODIR + "/bmdl/source/bmdl < " + job + ".bmdl"
        print("running:          " + command)
        starttime = time.time()
        os.system(command)
        endtime = time.time()
        timetaken = endtime - starttime
        timehours = int(timetaken // 3600)
        timeminutes = int((timetaken - timehours * 3600) // 60)
        timeseconds = (timetaken - timehours * 3600) - timeminutes * 60
        print("Finished running: " + command + '. Time: {0:3}h {1:2}m {2:5.2f}s '
              .format(timehours, timeminutes, timeseconds) + '\n')

        # Make sure folders exist
        common.check_folders(folder, folder + "/kstr/")
        # Create input file
        self.lattice.kstr.write_input_file(folder=folder)
        # Run KSTR
        job = os.path.join(folder, jobname)
        command = "cd {0}; ".format(
            folder) + EMTODIR + "/kstr/source/kstr < " + job + ".kstr"
        print("running:          " + command)
        starttime = time.time()
        os.system(command)
        endtime = time.time()
        timetaken = endtime - starttime
        timehours = int(timetaken // 3600)
        timeminutes = int((timetaken - timehours * 3600) // 60)
        timeseconds = (timetaken - timehours * 3600) - timeminutes * 60
        print("Finished running: " + command + '. Time: {0:3}h {1:2}m {2:5.2f}s '
              .format(timehours, timeminutes, timeseconds) + '\n')

        if self.lattice.kstr.twocenter == True:
            command = "cd {0}; ".format(
                folder) + EMTODIR + "/kstr/source/kstr < " + job + "2.kstr"
            print("running:          " + command)
            starttime = time.time()
            os.system(command)
            endtime = time.time()
            timetaken = endtime - starttime
            timehours = int(timetaken // 3600)
            timeminutes = int((timetaken - timehours * 3600) // 60)
            timeseconds = (timetaken - timehours * 3600) - timeminutes * 60
            print("Finished running: " + command + '. Time: {0:3}h {1:2}m {2:5.2f}s '
                  .format(timehours, timeminutes, timeseconds) + '\n')

        # Make sure folders exist
        common.check_folders(folder, folder + "/shape/")
        # Create input file
        self.lattice.shape.write_input_file(folder=folder)
        # Run SHAPE
        job = os.path.join(folder, jobname)
        command = "cd {0}; ".format(
            folder) + EMTODIR + "/shape/source/shape < " + job + ".shape"
        print("running:          " + command)
        starttime = time.time()
        os.system(command)
        endtime = time.time()
        timetaken = endtime - starttime
        timehours = int(timetaken // 3600)
        timeminutes = int((timetaken - timehours * 3600) // 60)
        timeseconds = (timetaken - timehours * 3600) - timeminutes * 60
        print("Finished running: " + command + '. Time: {0:3}h {1:2}m {2:5.2f}s '
              .format(timehours, timeminutes, timeseconds) + '\n')

        return

    def runemto(self, jobname, folder="./", EMTODIR=None, onlyKFCD=False):
        """Run KGRN and KFCD **WITHOUT** using the batch system and check convergence.

        
        :param jobname: Name of the input files
        :type jobname: str
        :param folder: Name of the folder where the input files are located (Default value = "./")
        :type folder: str
        :param EMTODIR: Path to the EMTO installation folder (Default value = None)
        :type EMTODIR: str
        :param onlyKFCD: True if only KFCD needs to calculated,
                         False if also KGRN (Default value = False)
        :type onlyKFCD: boolean
        :returns: True if the calculations converged, False if not
        :rtype: boolean
        """

        if jobname is None:
            sys.exit("System.runemto: \'jobname\' has to be given!")

        if EMTODIR is None:
            EMTODIR = self.EMTOdir

        if onlyKFCD:
            # Make sure folders exist
            common.check_folders(folder, folder + "/kfcd")
            # Run KFCD
            self.emto.kfcd.write_input_file(folder=folder)
            job = os.path.join(folder, jobname)
            command = "cd {0}; ".format(
                folder) + EMTODIR + "/kfcd/source/kfcd_cpa < " + job + ".kfcd"
            print("running:          " + command)
            starttime = time.time()
            os.system(command)
            endtime = time.time()
            timetaken = endtime - starttime
            timehours = int(timetaken // 3600)
            timeminutes = int((timetaken - timehours * 3600) // 60)
            timeseconds = (timetaken - timehours * 3600) - timeminutes * 60
            print("Finished running: " + command + '. Time: {0:3}h {1:2}m {2:5.2f}s '
                  .format(timehours, timeminutes, timeseconds) + '\n')

            # Check if we got convergence
            converged = self.check_conv(jobname, folder=folder)

        else:
            # Make sure folders exist
            common.check_folders(folder, folder + "/kgrn", folder + "/kgrn/tmp")
            # Run KGRN
            self.emto.kgrn.write_input_file(folder=folder)
            job = os.path.join(folder, jobname)
            command = "cd {0}; ".format(
                folder) + EMTODIR + "/kgrn/source/kgrn_cpa < " + job + ".kgrn"
            print("running:          " + command)
            starttime = time.time()
            os.system(command)
            endtime = time.time()
            timetaken = endtime - starttime
            timehours = int(timetaken // 3600)
            timeminutes = int((timetaken - timehours * 3600) // 60)
            timeseconds = (timetaken - timehours * 3600) - timeminutes * 60
            print("Finished running: " + command + '. Time: {0:3}h {1:2}m {2:5.2f}s '
                  .format(timehours, timeminutes, timeseconds))

            # Check if we got convergence
            converged = self.check_conv(jobname, folder=folder)

            if converged[0]:
                # Make sure folders exist
                common.check_folders(folder, folder + "/kfcd")
                # Run KFCD
                self.emto.kfcd.write_input_file(folder=folder)
                job = os.path.join(folder, jobname)
                command = "cd {0}; ".format(folder) + EMTODIR +\
                          "/kfcd/source/kfcd_cpa < " + job + ".kfcd"
                print("running:          " + command)
                starttime = time.time()
                os.system(command)
                endtime = time.time()
                timetaken = endtime - starttime
                timehours = int(timetaken // 3600)
                timeminutes = int((timetaken - timehours * 3600) // 60)
                timeseconds = (timetaken - timehours * 3600) - timeminutes * 60
                print("Finished running: " + command + '. Time: {0:3}h {1:2}m {2:5.2f}s '
                      .format(timehours, timeminutes, timeseconds) + '\n')
                converged = self.check_conv(jobname, folder=folder)

        return converged

    def which_error(self, jobname, folder="./"):
        """Tries to determine the reason why a given KGRN calculation did not converge.

        The reason for the crash will be printed on screen.
        
        :param jobname: Name of the KGRN output file
        :type jobname: str
        :param folder: Name of the folder where the output file is located (Default value = "./")
        :type folder: str
        :returns: None
        :rtype: None
        """

        folder = folder + '/kgrn/'

        errdict = {'EFXMOM': 'EFXMOM:** Fermi level not found',
                   'NOTCONV': 'Not converged'}

        prntfile = jobname + ".prn"

        fn = os.path.join(folder, prntfile)
        if not os.path.isfile(fn):
            print('Function which_error: file {0} does not exist!'.format(fn))
            quit()
        pfile = open(fn, "r")
        done = False
        # knownError = failure caused by a known reason in the errdict.
        # Otherwise cause is a more exotic KGRN error or the run crashed.
        knownError = False
        for line in pfile:
            for key in errdict:
                if errdict[key] in line:
                    errorkey = key
                    errormesg = line
                    done = True
                    knownError = True
                    break
            if done:
                break
        pfile.close()
        if knownError == True:
            print('Problem in KGRN for {0}: {1}'.format(fn, errormesg))
        else:
            print('Problem in KGRN for {0}: Program has crashed or stopped due to a' +
                  ' less common KGRN error. Check the output file.'.format(fn))
        return

    def get_energy(self, jobname=None, func="PBE", folder=None):
        """Extracts total energy from the KFCD output file.

        Different total energies given by different xc-functionals can
        be selected using the *func* input parameter. Default value is 
        'PBE'.

        :param jobname: Name of the KFCD output file
        :type jobname: str
        :param func: Name of the xc-functional (Default value = "PBE")
        :type func: str
        :param folder: Name of the folder where the output file is located (Default value = "./")
        :type folder: str
        :returns: Total energy if it is found, otherwise return None
        :rtype: float or None
        """

        if folder == None:
            folder = self.folder
        if jobname == None:
            jobname = self.fulljobname
        
        fn = os.path.join(folder, "kfcd/")
        fn = os.path.join(fn, jobname + ".prn")
        try:
            kfcdfile = open(fn, "r")
            energy = "TOT-" + func
            energyFound = False
            for line in kfcdfile:
                if energy in line:
                    energyFound = True
                    break
            # Total energy of the unit cell
            # return float(line.split()[1])
            # Energy per WS-radius: Better to use this to get Bmod correctly
            if energyFound:
                return float(line.split()[3])
            else:
                return
        except IOError:
            print('System.get_energy(): {0} does not exist!'.format(fn))
            return

    def get_moments(self, jobname=None, func="PBE", folder=None):
        """Extracts magnetic moments from the KFCD output file.

        :param jobname: Name of the KFCD output file
        :type jobname: str
        :param func: Name of the xc-functional (Default value = "PBE")
        :type func: str
        :param folder: Name of the folder where the output file is located (Default value = "./")
        :type folder: str
        :returns: Total energy if it is found, otherwise return None
        :rtype: float or None
        """

        if folder == None:
            folder = self.folder
        if jobname == None:
            jobname = self.fulljobname
        
        fn = os.path.join(folder, "kfcd/")
        fn = os.path.join(fn, jobname + ".prn")

        readyTag = "KFCD: Finished at:"
        momTag = "Magnetic moment for IQ ="

        try:
            kfcdfile = open(fn, "r")
            lines = kfcdfile.readlines()
            kfcdfile.close()
            moms = []
            ready = False
            for line in lines:
                if readyTag in line:
                    ready = True
            if ready:
                for line in lines:
                    #if enTag in line:
                    #    linesplit = line.split()
                    #    sws = float(linesplit[6])
                    if momTag in line:
                        linesplit = line.split()
                        moms.append(float(linesplit[7]))
                return moms
            else:
                return None
        except IOError:
            print('System.get_moments(): {0} does not exist!'.format(fn))
            return None

    def create_jobname(self, jobname=None):
        """Creates jobnames based on system information.

        Creates a jobname based on the optional input prefix *jobname*.
        The format of the full jobname in this case will be:
        jobname_3.000000, where 3.000000 is the sws.

        If *jobname* prefix is not given, full jobname is generated automatically
        based on system data self.sws, self.atoms and self.concs.
        The format of the jobname is: au0.50pd0.50_3.000000, where 0.50 are the
        atomic concentrations and 3.000000 is the sws.

        :param jobname: Optional prefix for the full jobname (Default value = None)
        :type jobname: str
        :returns: Newly created full jobname
        :rtype: str
        """

        if jobname is None:
            # Prefix not specified => create a default jobname based on
            # the types of atoms and their concentrations.
            jobname = ''
            for i in range(len(self.atoms)):
                if jobname == "":
                    pass
                else:
                    jobname = jobname + "_"
                jobname = jobname + \
                    self.atoms[i].lower() + "%4.2f" % (self.concs[i])
            fulljobname = jobname + "_%8.6f" % (self.sws)
            return jobname, fulljobname

        else:
            fulljobname = jobname + "_%8.6f" % (self.sws)

        return fulljobname

    def check_str(self, jobname, folder="./"):
        """Checks if a KSTR calculation file exists and has converged.

        Only KSTR file is checked because BMDL and SHAPE are fast to
        rerun in any case.

        :param jobname: Filename of the structure output file
        :type jobname: str
        :param folder: Folder where the output file is located (Default value = "./")
        :type folder: str
        :returns: True if KSTR has finished and False if not
        :rtype: boolean
        """

        folderKSTR = folder + '/kstr/'
        # Check if we got convergence in KSTR
        prntfile = jobname + ".prn"

        convergedKSTR = False
        fn = os.path.join(folderKSTR, prntfile)
        if not os.path.isfile(fn):
            pass
        else:
            pfile = open(fn, "r")
            for line in pfile:
                if "KSTR:     Finished at :" in line:
                    convergedKSTR = True
                    break
            pfile.close()

        return convergedKSTR

    def wait_for_jobs(self, jobsdict, restart_partition='general', sleeptime=60, restart_z=None,
                      restart_stragglers_after=0.75, kill_if_all_ssusp=False):
        """Loops checking status until no jobs are waiting or running / all are finished.
        
        wait/run states:

        ======= =========== ==================================================================================================
        **Key** **Meaning** **Description**
        ------- ----------- --------------------------------------------------------------------------------------------------
        CF      CONFIGURING Job has been allocated resources, but are waiting for them to become ready for use (e.g. booting).
        CG      COMPLETING  Job is in the process of completing. Some processes on some nodes may still be active.
        PD      PENDING     Job is awaiting resource allocation.
        R       RUNNING     Job currently has an allocation.
        RS      RESIZING    Job is about to change size.
        S       SUSPENDED   Job has an allocation, but execution has been suspended.
        ======= =========== ==================================================================================================

        done states:

        ======= =========== ==============================================================================================================
        **Key** **Meaning** **Description**
        ------- ----------- --------------------------------------------------------------------------------------------------------------
        CA      CANCELLED   Job was explicitly cancelled by the user or system administrator.  The job may or may not have been initiated.
        CD      COMPLETED   Job has terminated all processes on all nodes.
        F       FAILED      Job terminated with non-zero exit code or other failure condition.
        NF      NODE_FAIL   Job terminated due to failure of one or more allocated nodes.
        PR      PREEMPTED   Job terminated due to preemption.
        TO      TIMEOUT     Job terminated upon reaching its time limit.
        ======= =========== ==============================================================================================================

        :param jobsdict:
        :type jobsdict:
        :param restart_partition:  (Default value = 'general')
        :type restart_partition:
        :param sleeptime:  (Default value = 60)
        :type sleeptime:
        :param restart_z:  (Default value = None)
        :type restart_z:
        :param restart_stragglers_after:  (Default value = 0.75)
        :type restart_stragglers_after:
        :param kill_if_all_ssusp:  (Default value = False)
        :type kill_if_all_ssusp:
        :returns: None
        :rtype: None
        """

        run_status = ['CONFIGURING', 'COMPLETING',
                      'PENDING', 'RUNNING', 'RESIZING', 'SUSPENDED']
        done_status = ['CANCELLED', 'COMPLETED',
                       'FAILED', 'NODE_FAIL', 'PREEMPTED', 'TIMEOUT']

        import time
        import datetime
        jobs_amount = len(jobsdict)
        print()
        print('wait_for_jobs: Submitted {0} jobs'.format(jobs_amount))
        status = self.get_status_counts(jobsdict)
        t = time.time()
        maxllen = 0

        print('wait_for_jobs: Will be requesting job statuses' +
              ' every {0} seconds'.format(sleeptime) + "\n")

        while any([k in run_status for k in status.keys()]):
            time.sleep(sleeptime)
            status = self.get_status_counts(jobsdict)
            pctdone = sum([status.get(rs, 0)
                           for rs in done_status]) / float(sum(status.values()))

            # CHECK SUSPENDED; RESTART STRAGGLERS, ETC

            outl = '%s %s (%3d%% completion)' % (str(
                datetime.timedelta(seconds=int(time.time() - t))), status.__repr__(), pctdone * 100)
            # if len(outl) < maxllen:
            #    pad = maxllen - len(outl)
            #    outl += ' '*pad
            # else:
            #    maxllen = len(outl)

            print(outl)
            # sys.stderr.write(outl)
            # sys.stderr.flush()

        print('completed {0} batch jobs in {1}'.format(
            jobs_amount, str(datetime.timedelta(seconds=int(time.time() - t)))))
        return

    def get_status_counts(self, jobids=None):
        """Returns the counts of all jobs by status category.

        :param jobids:  (Default value = None)
        :type jobids:
        :returns: 
        :rtype:
        """
        from collections import defaultdict

        jobs_status = self.get_jobs_status(jobids)

        status_counts = defaultdict(int)
        for jd in jobs_status.values():
            status_counts[jd['State'].split()[0]] += 1

        return dict(status_counts)

    def get_jobs_status(self, jobids=None, toplevel=True):
        """Returns status of the jobs indicated
        (jobsdict or list of job ids) or all jobs if no jobids supplied.
        Set toplevel=False for job step data.

        :param jobids: List of job IDs (Default value = None)
        :type jobids: dict or list
        :param toplevel:  (Default value = True)
        :type toplevel: boolean
        :returns: Job statuses
        :rtype: list(str)
        """

        import subprocess

        if jobids is None:
            sacct_return = subprocess.Popen(
                'sacct -p -l', shell=True, stdout=subprocess.PIPE).stdout.readlines()
        else:
            if isinstance(jobids, dict):
                qjobs = jobids.keys()
            else:
                qjobs = jobids
            sacct_return = subprocess.Popen(
                'sacct -j %s -p -l' % (
                    ','.join(qjobs),), shell=True, stdout=subprocess.PIPE).stdout.readlines()

        jobs_status = {}
        for el in sacct_return[1:]:
            d = dict(
                zip(sacct_return[0].strip().split('|'), el.strip().split('|')))
            if not '.' in d['JobID'] or not toplevel:
                jobs_status[d['JobID']] = d

        return jobs_status

    def submit_jobs(self, jobnames, folder=None):
        """Takes a list of jobnames and submits the corresponding
        batch scripts to the batch system.

        :param jobnames: List of jobnames
        :type jobnames: list(float)
        :param folder: Folder where the batch scripts are located (Default value = None)
        :type folder: str
        :returns: List of job ID numbers
        :rtype: list(str)
        """

        import time
        from pyemto.utilities.utils import run_emto

        sleeptime = 10

        if folder is None:
            folder = self.folder

        job_ids = []
        for i in range(len(jobnames)):
            job_ids.append(run_emto(jobnames[i], folder=self.folder))

        # Flatten job_ids list and convert the integers into strings
        job_ids = [item for sublist in job_ids for item in sublist]
        for i in range(len(job_ids)):
            job_ids[i] = str(job_ids[i])

        # Give SLURM some time to register the jobs.
        # If continued straight away the self.wait_for_jobs
        # script will likely think all the jobs have finished
        # since it cannot find them yet.
        time.sleep(sleeptime)

        return job_ids
