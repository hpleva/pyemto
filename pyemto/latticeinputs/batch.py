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
    """Creates a batch script for running BMDL, KSTR and SHAPE calculations

    This class is used to to create batch scripts for a supercomputer environment (EMTO 5.8).
    !!! Currently only SLURM is supported. !!!

    :param jobname:  (Default value = None)
    :type jobname:
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

    def __init__(self, jobname=None, lat=None, runtime=None, latpath=None,
                 EMTOdir=None, runBMDL=None, runKSTR=None, runKSTR2=None,
                 runSHAPE=None, kappaw=None, kappalen=None):

        # Batch script related parameters
        self.jobname = jobname
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

    def output(self):
        """(self) -> (str)

            Output first part of the kgrn input file in formated string

        :returns:
        :rtype:
        """

        # Clean up path names

        line = "#!/bin/bash" + "\n" + "\n"
        line += "#SBATCH -J " + self.jobname + "\n"
        line += "#SBATCH -t " + self.runtime + "\n"
        line += "#SBATCH -o " + \
            common.cleanup_path(
                self.latpath + "/" + self.jobname) + ".output" + "\n"
        line += "#SBATCH -e " + \
            common.cleanup_path(
                self.latpath + "/" + self.jobname) + ".error" + "\n"
        # line += "#SBATCH -x pl1,pl11"+"\n"
        line += "\n"

        elapsed_time = "/usr/bin/time "
        if self.runBMDL:
            line += elapsed_time + common.cleanup_path(self.EMTOdir + "/bmdl/source/bmdl < ") +\
                common.cleanup_path(self.latpath + "/" + self.jobname) + ".bmdl > " +\
                common.cleanup_path(
                    self.latpath + "/" + self.jobname) + "_bmdl.output" + "\n"
        if self.runKSTR:
            line += elapsed_time + common.cleanup_path(self.EMTOdir + "/kstr/source/kstr < ") +\
                common.cleanup_path(self.latpath + "/" + self.jobname) + ".kstr > " +\
                common.cleanup_path(
                    self.latpath + "/" + self.jobname) + "_kstr.output" + "\n"
        if self.runKSTR2:
            line += elapsed_time + common.cleanup_path(self.EMTOdir + "/kstr/source/kstr < ") +\
                common.cleanup_path(self.latpath + "/" + self.jobname) + '2' + ".kstr > " +\
                common.cleanup_path(
                    self.latpath + "/" + self.jobname) + '2' + "_kstr.output" + "\n"
        if self.runSHAPE:
            line += elapsed_time + common.cleanup_path(self.EMTOdir + "/shape/source/shape < ") +\
                common.cleanup_path(self.latpath + "/" + self.jobname) + ".shape > " +\
                common.cleanup_path(
                    self.latpath + "/" + self.jobname) + "_shape.output" + "\n"

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

        fl = open(folder + '/{0}.cmd'.format(self.jobname), "w")
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
        if self.jobname is None:
            if self.lat is not None:
                self.jobname = self.lat
            else:
                sys.exit(
                    'Batch_lattice.check_input_file: \'jobname\' or' +\
                    ' \'lat\' (jobname = lat) has to be given!')

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
