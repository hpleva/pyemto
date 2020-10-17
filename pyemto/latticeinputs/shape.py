# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 15:10:00 2014

@author: Matti Ropo
@author: Henrik Levämäki

"""

from __future__ import print_function
import sys
import datetime
import numpy as np
import pyemto.common.common as common


class Shape:

    """Contains information about SHAPE input files for EMTO 5.8 program.

    :param jobname_lat:  (Default value = None)
    :type jobname_lat:
    :param lat:  (Default value = None)
    :type lat:
    :param lmax:  (Default value = None)
    :type lmax:
    :param nsr:  (Default value = None)
    :type nsr:
    :param nfi:  (Default value = None)
    :type nfi:
    :param ivef:  (Default value = None)
    :type ivef:
    :param msgl:  (Default value = None)
    :type msgl:
    :param nprn:  (Default value = None)
    :type nprn:
    :returns:
    :rtype:
    """

    def __init__(self, jobname_lat=None, lat=None, lmax=None, nsr=None, nfi=None,
                 ivef=None, msgl=None, nprn=None, basis=None, asr=None):

        self.jobname_lat = jobname_lat
        self.lat = lat
        self.lmax = lmax
        self.nsr = nsr
        self.nfi = nfi
        self.ivef = ivef
        self.msgl = msgl
        self.nprn = nprn
        if basis is None:
            self.basis = None
        else:
            self.basis = np.asarray(basis)
        self.asr = asr

    def output(self):
        """Output SHAPE input file in formatted string.

        :returns: SHAPE input file as a string.
        :rtype: str
        """

        slope = 'kstr/' + self.jobname_lat + ".tfh"
        shapef = "shape/"
        prn = "shape/"
        now = datetime.datetime.now()

        line = "SHAPE     HP......=N                            "\
            + str(now.day) + "." + str(now.month) + "." + str(now.year) + "\n"
        JOBNAMline = "JOBNAM...=" + self.jobname_lat
        MSGLline = "MSGL.=  " + str(self.msgl)
        line = line + "{0:21s}{1:9s}".format(JOBNAMline, MSGLline) + "\n"
        line = line + "FOR001=" + slope + "\n"
        line = line + "DIR002=" + shapef + "\n"
        line = line + "DIR006=" + prn + "\n"
        line = line + "Lmax..={0:3d} NSR..={1:3d} NFI..={2:3d}"\
            .format(self.lmax, self.nsr, self.nfi) + "\n"
        line = line + "NPRN..={0:3d} IVEF.={1:3d}"\
            .format(self.nprn, self.ivef) + "\n"
        line = line + "****** Relative atomic sphere radii ASR(1:NQ) ******\n"
        for i in range(self.nq):
            line += "ASR({0}).= {1:6.4f}\n".format(i+1, self.asr[i])
        return line

    def write_input_file(self, folder=None):
        """Save SHAPE input data to file named filename

        :param folder: directory to write  (Default value = None)
        :type folder:
        :returns:
        :rtype:
        """

        # Check data integrity before anything is written on disk or run
        self.check_input_file()

        if folder is None:
            #sys.exit('Shape.create_input_file: \'path\' has to be given!')
            folder = "./"
        else:
            common.check_folders(folder)

        fl = open(folder + '/{0}.shape'.format(self.jobname_lat), "w")
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
            print(
                'WARNING: Shape() class has no attribute \'{0}\''.format(key))
        return

    def check_input_file(self):
        """Perform various checks on the class data to
            make sure that all necessary data exists
            before we attempt to write the input file to disk

        :returns:
        :rtype:
        """

        if self.jobname_lat is None and self.lat is None:
            sys.exit(
                'Shape.check_input_file: \'jobname_lat\' OR \'lat\' has to be given!')
        elif self.jobname_lat is None and self.lat is not None:
            self.jobname_lat = self.lat

        if self.lmax is None:
            self.lmax = 30
        if self.nsr is None:
            self.nsr = 129
        if self.nfi is None:
            self.nfi = 11
        if self.ivef is None:
            self.ivef = 3
        if self.msgl is None:
            self.msgl = 1
        if self.nprn is None:
            self.nprn = 0

        if self.basis is None:
            if self.lat == 'hcp':
                self.basis = np.asarray([
                    [0.0, 0.0, 0.0], [0.0, 0.57735027, self.ca / 2.0]])
            else:
                self.basis = np.asarray([0.0, 0.0, 0.0])

        # Make sure basis is a numpy array of shape np.array([[xxx], [xxx]])
        if isinstance(self.basis, list):
            self.basis = np.asarray(self.basis)
        elif isinstance(self.basis, np.ndarray):
            pass
        if len(self.basis.shape) == 1:
            self.basis = np.asarray([self.basis])
        self.nq = self.basis.shape[0]

        if self.asr is None:
            self.asr = np.ones(self.nq, dtype=float)

        return
