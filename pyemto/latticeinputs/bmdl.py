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
import pyemto.common.common as common
import numpy as np

class Bmdl:
    """Contains information about BMDL input files for EMTO 5.8 program.

    :param jobname:  (Default value = None)
    :type jobname:
    :param lat:  (Default value = None)
    :type lat:
    :param latparams:  (Default value = None)
    :type latparams:
    :param latvectors:  (Default value = None)
    :type latvectors:
    :param basis:  (Default value = None)
    :type basis:
    :param msgl:  (Default value = None)
    :type msgl:
    :param nprn:  (Default value = None)
    :type nprn:
    :param bmdl_nl:  (Default value = None)
    :type bmdl_nl:
    :param lamda:  (Default value = None)
    :type lamda:
    :param amax:  (Default value = None)
    :type amax:
    :param bmax:  (Default value = None)
    :type bmax:
    :param nqr2:  (Default value = None)
    :type nqr2:
    :param ca: hcp c/a ratio (Default value = None)
    :type ca: float
    :returns: None
    :rtype: None
    """

    def __init__(self, jobname=None, lat=None, latparams=None, latvectors=None,
                 basis=None, msgl=None, nprn=None, bmdl_nl=None, lamda=None,
                 amax=None, bmax=None, nqr2=None, ca=None):

        self.lat = lat
        self.jobname = jobname
        self.latparams = latparams
        self.latvectors = latvectors
        self.basis = np.asarray(basis)
        self.msgl = msgl
        self.nprn = nprn
        self.bmdl_nl = bmdl_nl
        self.lamda = lamda
        self.amax = amax
        self.bmax = bmax
        self.nqr2 = nqr2
        self.ca = ca  # hcp's c/a

    def output(self):
        """Outputs BMDL input file as a formatted string

        Creates a long string containing the batch
        job script which will later be written on disk.

        :returns: batch script input file
        :rtype: str
        """

        mdlfile = "bmdl/"
        prnfile = "bmdl/"

        now = datetime.datetime.now()
        line = "BMDL      HP......=N                              "\
            + str(now.day) + "." + str(now.month) + "." + str(now.year) + "\n"
        JOBNAMline = "JOBNAM...=" + self.jobname
        MSGLline = "MSGL.=  " + str(self.msgl)
        NPRNline = "NPRN.=  " + str(self.nprn)
        line = line + \
            "{0:21s}{1:9s} {2:9s}".format(
                JOBNAMline, MSGLline, NPRNline) + "\n"
        line = line + "DIR001=" + mdlfile + "\n"
        line = line + "DIR006=" + prnfile + "\n"
        line = line + "Madelung potential for {0}".format(self.jobname) + "\n"
        line = line + "NL.....={0:2d}".format(self.bmdl_nl) + "\n"
        line = line + "LAMDA....=      {0:4.2f} AMAX....=      {1:4.2f} BMAX....=      {2:4.2f}"\
            .format(self.lamda, self.amax, self.bmax) + "\n"
        line = line + "NQ....={0:3d} LAT...={1:2d} IPRIM.= {2} NQR2..= {3}"\
            .format(self.nq, common.lat_to_ibz(self.lat), self.iprim, self.nqr2) + "\n"
        line = line + "A........={0:10.8f} B.......={1:10.8f} C.......={2:10.8f}"\
            .format(self.latparams[0], self.latparams[1], self.latparams[2]) + "\n"
        if self.iprim == 1:
            line = line + "ALPHA....={0:10.6f} BETA....={1:10.6f} GAMMA...={2:10.6f}"\
                .format(self.latvectors[0], self.latvectors[1], self.latvectors[2]) + "\n"
        elif self.iprim == 0:
            line = line + "BSX......={0:10.7f} BSY.....={1:10.7f} BSZ.....={2:10.7f}"\
                .format(self.latvectors[0][0], self.latvectors[0][1], self.latvectors[0][2]) + "\n"
            line = line + "BSX......={0:10.7f} BSY.....={1:10.7f} BSZ.....={2:10.7f}"\
                .format(self.latvectors[1][0], self.latvectors[1][1], self.latvectors[1][2]) + "\n"
            line = line + "BSX......={0:10.7f} BSY.....={1:10.7f} BSZ.....={2:10.7f}"\
                .format(self.latvectors[2][0], self.latvectors[2][1], self.latvectors[2][2]) + "\n"
        for i in range(self.nq):
            line = line + "QX({3:03})..={0:10.7f} QY({3:03}).={1:10.7f} QZ({3:03}).={2:10.7f}"\
                .format(self.basis[i,0], self.basis[i,1], self.basis[i,2], i + 1) + "\n"

        return line

    def write_input_file(self, folder=None):
        """Save BMDL input data to file named by self.jobname

        :param folder:  (Default value = None)
        :type folder:
        :returns:
        :rtype:
        """

        # Check data integrity before anything is written on disk or run
        self.check_input_file()

        if folder is None:
            #sys.exit('Bmdl.write_input_file: \'folder\' has to be given!')
            folder = "./"
        else:
            common.check_folders(folder)

        fl = open(folder + '/{0}.bmdl'.format(self.jobname), "w")
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
            if key == 'ca' and self.lat == 'hcp':
                # c/a has changed => update lattice parameters
                self.latparams = [1.0, 1.0, self.ca]
                self.latvectors = [
                    [1.0, 0.0, 0.0], [-0.5, 0.8660254, 0.0], [0.0, 0.0, self.ca]]
                self.basis = np.asarray([
                    [0.0, 0.0, 0.0], [0.0, 0.57735027, self.ca / 2.0]])

        else:
            print('WARNING: Bmdl() class has no attribute \'{0}\''.format(key))
        return

    def check_input_file(self):
        """Perform various checks on the class data.

        Makes sure that all necessary data exists
        before we attempt to write the input file to disk.

        :returns: None
        :rtype: None
        """

        if self.ca is None and self.lat == 'hcp':
            #sys.exit('Bmdl.check_input_file: for hcp \'ca\' = c/a has to be given!')
            self.ca = 1.632993

        if self.lat is None:
            sys.exit('Bmdl.check_input_file: \'lat\' has to be given!')
        elif self.jobname is None and self.lat is not None:
            self.jobname = self.lat

        if self.latparams is None:
            if self.lat == 'hcp':
                self.latparams = [1.0, 1.0, self.ca]
            else:
                self.latparams = [1.0, 1.0, 1.0]

        if self.latvectors is None:
            if self.lat == 'sc':
                self.latvectors = [
                    [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
            elif self.lat == 'bcc':
                self.latvectors = [
                    [0.5, 0.5, -0.5], [0.5, -0.5, 0.5], [-0.5, 0.5, 0.5]]
            elif self.lat == 'fcc':
                self.latvectors = [
                    [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]
            elif self.lat == 'hcp':
                self.latvectors = [
                    [1.0, 0.0, 0.0], [-0.5, 0.8660254, 0.0], [0.0, 0.0, self.ca]]

        if self.basis is None:
            if self.lat == 'hcp':
                self.basis = np.asarray([
                    [0.0, 0.0, 0.0], [0.0, 0.57735027, self.ca / 2.0]])
            else:
                self.basis = np.array([[0.0, 0.0, 0.0]])

        #if isinstance(self.basis[0], list):
        #    self.nq = len(self.basis)
        if type(self.basis[0]) == type(np.array([0.0,0.0,0.0])):
            self.nq = len(self.basis)
        else:
            self.nq = 1
            self.basis = np.asarray([self.basis])
        if isinstance(self.latvectors[0], list):
            if len(self.latvectors) == 1:
                self.iprim = 1
            else:
                self.iprim = 0
        else:
            self.iprim = 1

        if self.msgl is None:
            self.msgl = 1
        if self.nprn is None:
            self.nprn = 0
        if self.bmdl_nl is None:
            self.bmdl_nl = 7
        if self.lamda is None:
            self.lamda = 2.5
        if self.amax is None:
            self.amax = 4.5
        if self.bmax is None:
            self.bmax = 4.5
        if self.nqr2 is None:
            self.nqr2 = 0
        return

