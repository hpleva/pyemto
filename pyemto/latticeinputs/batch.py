# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 15:09:24 2014

@author: Matti Ropo
@author: Henrik Levämäki

"""

from __future__ import print_function
import sys
import pyemto.common.common as common


class Batch:
    """Creates a batch script for running BMDL, KSTR and SHAPE calculations

    This class is used to to create batch scripts for a supercomputer environment (EMTO 5.8).
    !!! Currently only SLURM is supported. !!!

    :param jobname_lat:  (Default value = None)
    :type jobname_lat:
    :param lat:  (Default value = None)
    :type lat:
    :param runtime:  (Default value = None)
    :type runtime:
    :param latpath:  (Default value = None)
    :type latpath:
    :param EMTOdir:  (Default value = None)
    :type EMTOdir:
    :param runBMDL:  (Default value = None)
    :type runBMDL:
    :param runKSTR:  (Default value = None)
    :type runKSTR:
    :param runKSTR2:  (Default value = None)
    :type runKSTR2:
    :param runSHAPE:  (Default value = None)
    :type runSHAPE:
    :param kappaw:  (Default value = None)
    :type kappaw:
    :param kappalen:  (Default value = None)
    :type kappalen:
    :returns: None
    :rtype: None
    """

    def __init__(self, jobname_lat=None, lat=None, runtime=None, latpath=None,
                 EMTOdir=None, runBMDL=None, runKSTR=None, runKSTR2=None,
                 runSHAPE=None, kappaw=None, kappalen=None,
                 slurm_options=None, account=None):

        # Batch script related parameters
        self.jobname_lat = jobname_lat
        self.lat = lat
        self.latpath = latpath
        self.runtime = runtime
        self.EMTOdir = EMTOdir
        self.runBMDL = runBMDL
        self.runKSTR = runKSTR
        self.runKSTR2 = runKSTR2
        self.runSHAPE = runSHAPE
        self.kappaw = kappaw
        self.kappalen = kappalen
        self.account = account
        self.slurm_options = slurm_options
        self.use_module = False
        #print('BMDL self.slurm_options = ',self.slurm_options)

    def output(self):
        """(self) -> (str)

            Output first part of the kgrn input file in formated string

        :returns:
        :rtype:
        """

        # Clean up path names

        line = "#!/bin/bash" + "\n" + "\n"
        line += "#SBATCH -J " + self.jobname_lat + "\n"
        line += "#SBATCH -t " + self.runtime + "\n"
        line += "#SBATCH -o " + \
            common.cleanup_path(
                self.latpath + "/" + self.jobname_lat) + ".output" + "\n"
        line += "#SBATCH -e " + \
            common.cleanup_path(
                self.latpath + "/" + self.jobname_lat) + ".error" + "\n"
        if self.account is not None:
            line += "#SBATCH -A {0}".format(self.account) + "\n"

        self.use_module = False
        if self.slurm_options is not None:
            for tmp in self.slurm_options:
                if 'module load emto' in tmp:
                    self.use_module = True
                    break
            for so in self.slurm_options:
                # Do not use more than one CPU for the structure calculations
                if "#SBATCH -n " in so:
                    pass
                else:
                    line += so + "\n"
        line += "\n"

        #elapsed_time = "/usr/bin/time "
        elapsed_time = ""

        if not self.use_module:
            # BMDL_path = self.EMTOdir + "/bmdl/source/bmdl"
            # KSTR_path = self.EMTOdir + "/kstr/source/kstr"
            # SHAPE_path = self.EMTOdir + "/shape/source/shape"
            BMDL_path = self.EMTOdir + "/bmdl/bmdl"
            KSTR_path = self.EMTOdir + "/kstr/kstr"
            SHAPE_path = self.EMTOdir + "/shape/shape"
        else:
            BMDL_path = "bmdl"
            KSTR_path = "kstr"
            SHAPE_path = "shape"
        if self.runBMDL:
            line += elapsed_time + common.cleanup_path(BMDL_path + " < ") +\
                common.cleanup_path(self.latpath + "/" + self.jobname_lat) + ".bmdl > " +\
                common.cleanup_path(
                    self.latpath + "/" + self.jobname_lat) + "_bmdl.output" + "\n"
        if self.runKSTR:
            line += elapsed_time + common.cleanup_path(KSTR_path + " < ") +\
                common.cleanup_path(self.latpath + "/" + self.jobname_lat) + ".kstr > " +\
                common.cleanup_path(
                    self.latpath + "/" + self.jobname_lat) + "_kstr.output" + "\n"
        if self.runKSTR2:
            line += elapsed_time + common.cleanup_path(KSTR_path + " < ") +\
                common.cleanup_path(self.latpath + "/" + self.jobname_lat) + 'M' + ".kstr > " +\
                common.cleanup_path(
                    self.latpath + "/" + self.jobname_lat) + 'M' + "_kstr.output" + "\n"
        if self.runSHAPE:
            line += elapsed_time + common.cleanup_path(SHAPE_path + " < ") +\
                common.cleanup_path(self.latpath + "/" + self.jobname_lat) + ".shape > " +\
                common.cleanup_path(
                    self.latpath + "/" + self.jobname_lat) + "_shape.output" + "\n"

        return line

    def write_input_file(self, folder=None):
        """(self,str) ->(None)

            Save batch input data to file named filename

        :param folder:  (Default value = None)
        :type folder:
        :returns:
        :rtype:
        """

        # Check data integrity before anything is written on disk or run
        self.check_input_file()

        if folder is None:
            #sys.exit('Batch_lattice.write_input_file: \'folder\' has to be given!')
            folder = "./"
        else:
            common.check_folders(folder)

        fl = open(folder + '/{0}.sh'.format(self.jobname_lat), "w")
        fl.write(self.output())
        fl.close()

    def set_values(self, key, value):
        """

        :param key:
        :type key:
        :param value:
        :type value:
        :returns:
        :rtype:
        """

        if hasattr(self, key):
            setattr(self, key, value)
        else:
            print('WARNING: Batch_lattice() class has no attribute \'{0}\''.format(key))
        return

    def check_input_file(self):
        """Perform various checks on the class data to
            make sure that all necessary data exists
            before we attempt to write the input file to disk

        :returns:
        :rtype:
        """

        # Mission critical parameters
        if self.jobname_lat is None:
            if self.lat is not None:
                self.jobname_lat = self.lat
            else:
                sys.exit(
                    'Batch_lattice.check_input_file: \'jobname_lat\' or' +\
                    ' \'lat\' (jobname_lat = lat) has to be given!')

        if self.latpath is None:
            self.latpath = "./"
        if self.runtime is None:
            self.runtime = "01:00:00"
        if self.EMTOdir is None:
            self.EMTOdir = "$HOME/EMTO5.8/"
        if self.runBMDL is None:
            self.runBMDL = True
        if self.runKSTR is None:
            self.runKSTR = True
        if self.runSHAPE is None:
            self.runSHAPE = True
        if self.kappaw is None:
            self.kappaw = [0.0]

        self.kappalen = len(self.kappaw)

        if self.kappalen == 2:
            self.runKSTR2 = True
        return
