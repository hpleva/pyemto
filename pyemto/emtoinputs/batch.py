# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 15:09:24 2014

@author: Matti Ropo
@author: Henrik Levämäki

"""

from __future__ import print_function
import sys
import os
import datetime
import re
import pyemto.common.common as common

class Batch:
    """Creates a batch script for running KGRN and KFCD calculations on a
    supercomputer environment (EMTO 5.8).

    !!! Currently only SLURM is supported. !!!

    :param jobname: Name for the KGRN and KFCD jobs. This will become
                    the first part of the input and output file names.
    :type jobname: str
    :param runtime: Maximum running time for the individual batch jobs.
                    The format of this entry should be 'xx:yy:zz', where
                    xx is hours, yy is minutes and zz is seconds.
    :type runtime: str
    :param EMTOdir: Path to the EMTO installation (Default value = '$HOME/EMTO5.8')
    :type EMTOdir: str
    :param emtopath: Path to the folder where the KGRN and KFCD input files are
                     located
    :type emtopath: str
    :param runKGRN: True if KGRN should be run, False if KGRN should not be run
    :type runKGRN: boolean
    :param runKFCD: True if KFCD should be run, False if KFCD should not be run
    :type runKFCD: boolean
    :returns: None
    :rtype: None
    """

    def __init__(self, jobname=None, runtime=None, EMTOdir=None,
                 emtopath=None, runKGRN=None, runKFCD=None,
                 account=None, KGRN_file_type=None, KFCD_file_type=None,
                 slurm_options=None, parallel=None):

        # Batch script related parameters
        self.jobname = jobname
        self.runtime = runtime
        self.emtopath = emtopath
        self.EMTOdir = EMTOdir
        self.runKGRN = runKGRN
        self.runKFCD = runKFCD
        self.account = account
        if KGRN_file_type is not None:
            self.KGRN_file_type = KGRN_file_type
        else:
            self.KGRN_file_type = 'kgrn'
        if KFCD_file_type is not None:
            self.KFCD_file_type = KFCD_file_type
        else:
            self.KFCD_file_type = 'kfcd'
        self.slurm_options = slurm_options
        self.parallel = parallel
        self.use_module = False
        return

    def output(self):
        """
        Output first part of the kgrn input file in formated string

        :returns: Batch job script file in the form of a long string
        :rtype: str
        """

        line = "#!/bin/bash" + "\n" + "\n"
        line += "#SBATCH -J " + self.jobname + "\n"
        line += "#SBATCH -t " + self.runtime + "\n"
        line += "#SBATCH -o " + \
            common.cleanup_path(
                self.emtopath + "/" + self.jobname) + ".output" + "\n"
        line += "#SBATCH -e " + \
            common.cleanup_path(
                self.emtopath + "/" + self.jobname) + ".error" + "\n"
        if self.account is not None:
            line += "#SBATCH -A {0}".format(self.account) + "\n"
        
        self.use_module = False
        if self.slurm_options is not None:
            for tmp in self.slurm_options:
                if 'module load emto' in tmp:
                    self.use_module = True
                    break
            for so in self.slurm_options:
                line += so + "\n"
        line += "\n"

        elapsed_time = "/usr/bin/time "
        if self.parallel is True:
            kgrn_exe = 'kgrn_omp'
            kfcd_exe = 'kfcd_cpa'
        else:
            kgrn_exe = 'kgrn_cpa'
            kfcd_exe = 'kfcd_cpa'

        if not self.use_module:
            KGRN_path = self.EMTOdir + "/kgrn/source/"
            KFCD_path = self.EMTOdir + "/kfcd/source/"
        else:
            KGRN_path = ""
            KFCD_path = ""

        if self.runKGRN:
            line += elapsed_time + common.cleanup_path(KGRN_path + kgrn_exe + " < ") +\
                    common.cleanup_path(self.emtopath + "/" + self.jobname) + ".{0} > ".format(self.KGRN_file_type) +\
                    common.cleanup_path(self.emtopath + "/" + self.jobname) + "_kgrn.output" + "\n"
        if self.runKFCD:
            line += elapsed_time + common.cleanup_path(KFCD_path + kfcd_exe + " < ") +\
                common.cleanup_path(self.emtopath + "/" + self.jobname) + ".{0} > ".format(self.KFCD_file_type) +\
                common.cleanup_path(self.emtopath + "/" + self.jobname) + "_kfcd.output" + "\n"

        return line

    def write_input_file(self, folder=None):
        """(self,str) ->(None)

            Save BMDL input data to file named filename

        :param folder:  (Default value = None)
        :type folder:
        :returns:
        :rtype:
        """

        # Check data integrity before anything is written on disk or run
        self.check_input_file()

        if folder is None:
            #sys.exit('Batch_emto.write_input_file: \'folder\' has to be given!')
            folder = "./"
        else:
            common.check_folders(folder)

        fl = open(folder + '/{0}.sh'.format(self.jobname), "w")
        fl.write(self.output())
        fl.close()
        return

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
            print('WARNING: Batch_emto() class has no attribute \'{0}\''.format(key))
        return

    def check_input_file(self):
        """Perform various checks on the class data to
            make sure that all necessary data exists
            before we attempt to write the input file to disk
        """

        # Mission critical parameters
        if self.jobname is None:
            sys.exit('Batch_emto.check_input_file: \'jobname\' has to be given!')

        if self.runtime is None:
            self.runtime = "02:00:00"
        if self.emtopath is None:
            self.emtopath = "./"
        if self.EMTOdir is None:
            self.EMTOdir = "$HOME/EMTO5.8/"
        if self.runKGRN is None:
            self.runKGRN = True
        if self.runKFCD is None:
            self.runKFCD = True
        return
