# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 15:00:14 2014

@author: Matti Ropo
@author: Henrik Levämäki

"""

from __future__ import print_function
import sys
import os
import datetime
import numpy as np
import pyemto.common.common as common

class Kstr:
    """Contains information about KSTR input files for EMTO 5.8 program

    :param jobname:  (Default value = None)
    :type jobname:
    :param lat:  (Default value = None)
    :type lat:
    :param latparams:  (Default value = None)
    :type latparams:
    :param ca:  (Default value = None)
    :type ca:
    :param latvectors:  (Default value = None)
    :type latvectors:
    :param basis:  (Default value = None)
    :type basis:
    :param kappaw:  (Default value = None)
    :type kappaw:
    :param dmax:  (Default value = None)
    :type dmax:
    :param msgl:  (Default value = None)
    :type msgl:
    :param nprn:  (Default value = None)
    :type nprn:
    :param lamda:  (Default value = None)
    :type lamda:
    :param amax:  (Default value = None)
    :type amax:
    :param bmax:  (Default value = None)
    :type bmax:
    :param nqr2:  (Default value = None)
    :type nqr2:
    :param mode:  (Default value = None)
    :type mode:
    :param store:  (Default value = None)
    :type store:
    :param high:  (Default value = None)
    :type high:
    :param kstr_nl:  (Default value = None)
    :type kstr_nl:
    :param nlh:  (Default value = None)
    :type nlh:
    :param nlw:  (Default value = None)
    :type nlw:
    :param nder:  (Default value = None)
    :type nder:
    :param itrans:  (Default value = None)
    :type itrans:
    :param rwats:  (Default value = None)
    :type rwats:
    :param nghbp:  (Default value = None)
    :type nghbp:
    :param awIQ:  (Default value = None)
    :type awIQ:
    :returns: None
    :rtype: None
    """

    def __init__(self, jobname=None, lat=None, latparams=None, ca=None,
                 latvectors=None, basis=None, kappaw=None, dmax=None,
                 msgl=None, nprn=None, lamda=None, amax=None,
                 bmax=None, nqr2=None, mode=None, store=None, high=None,
                 kstr_nl=None, nlh=None, nlw=None, nder=None, itrans=None,
                 rwats=None, nghbp=None, awIQ=None):

        # Argument checking and default values

        self.lat = lat
        self.jobname = jobname
        self.latparams = latparams
        self.latvectors = latvectors
        self.basis = basis
        self.kappaw = kappaw
        self.dmax = dmax
        self.msgl = msgl
        self.nprn = nprn
        self.lamda = lamda
        self.amax = amax
        self.bmax = bmax
        self.nqr2 = nqr2
        self.mode = mode
        self.store = store
        self.high = high
        self.kstr_nl = kstr_nl
        self.nlh = nlh
        self.nlw = nlw
        self.nder = nder
        self.itrans = itrans
        self.rwats = rwats
        self.nghbp = nghbp
        self.ca = ca  # hcp's c/a ratio
        self.awIQ = awIQ

    def output(self, index):
        """Outputs KSTR input file as a formatted string.

        :param index:
        :type index:
        :returns: KSTR input file string
        :rtype: str
        """

        slope = 'kstr/'
        prn = 'kstr/'
        boolean = {True: "Y", False: "N"}

        if index == 1:
            kappaw = self.kappaw[0]
            jobname = self.jobname
        elif index == 2:
            kappaw = self.kappaw[1]
            jobname = self.jobname2

        now = datetime.datetime.now()
        line = "KSTR      HP......=N                              "\
            + str(now.day) + "." + str(now.month) + "." + str(now.year) + "\n"
        JOBNAMline = "JOBNAM...=" + jobname
        MSGLline = "MSGL.=  " + str(self.msgl)
        line = line + "{0:21s}{1:9s}".format(JOBNAMline, MSGLline) +\
            " MODE...={0} STORE..={1} HIGH...={2}"\
            .format(self.mode, self.store, self.high) + "\n"
        line = line + "DIR001=" + slope + "\n"
        line = line + "DIR006=" + prn + "\n"
        line = line + "Slope matrices, {0:10}, (kappa*w)^2= {1:5.1f}"\
            .format(jobname, kappaw) + "\n"
        line = line + "NL.....={0:2d} NLH...={1:2d} NLW...={2:2d} "\
            .format(self.kstr_nl, self.nlh, self.nlw) +\
            "NDER..={0:2d} ITRANS={1:2d} NPRN..={2:2d}"\
            .format(self.nder, self.itrans, self.nprn) + "\n"
        line = line + "(K*W)^2..={0:10.6f} DMAX....={1:10.4f} RWATS...={2:10.2f}"\
            .format(kappaw, self.dmax, self.rwats) + "\n"
        line = line + "NQ3...={0:3d} LAT...={1:2d} IPRIM.={2:2d} NGHBP.={3:2d} NQR2..={4:2d}"\
            .format(self.nq, common.lat_to_ibz(self.lat), self.iprim, self.nghbp, self.nqr2) + "\n"
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
            line = line + "QX.......={0:10.7f} QY......={1:10.7f} QZ......={2:10.7f}"\
                .format(self.basis[i][0], self.basis[i][1], self.basis[i][2]) + "\n"
        for i in range(self.nq):
            line = line + "a/w(IQ)..= {0:4.2f} {1:4.2f} {2:4.2f} {3:4.2f}"\
                .format(*self.awIQ[i, :]) + "\n"
        line = line + "LAMDA....={0:10.7f} AMAX....={1:10.7f} BMAX....={2:10.7f}"\
            .format(self.lamda, self.amax, self.bmax) + "\n"

        return line

    def write_input_file(self, folder=None):
        """Save KSTR input data to file named filename.

        :param folder:  (Default value = None)
        :type folder:
        :returns: None
        :rtype: None
        """

        # Check data integrity before anything is written on disk or run
        self.check_input_file()

        if folder is None:
            #sys.exit('Kstr.create_input_file: \'path\' has to be given!')
            folder = "./"
        else:
            common.check_folders(folder)

        fl = open(folder + '/{0}.kstr'.format(self.jobname), "w")
        fl.write(self.output(1))
        fl.close()
        if self.twocenter:
            fl = open(folder + '/{0}2.kstr'.format(self.jobname), "w")
            fl.write(self.output(2))
            fl.close()

        # Re-initialize some of the input parameters
        # to avoid situations, where an input parameter
        # doesn't get updated correctly, because it already
        # existed.
        self.finalize()
        return

    def set_values(self, key, value):
        """Set input parameter values.

        :param key:
        :type key:
        :param value:
        :type value:
        :returns: None
        :rtype: None
        """

        if hasattr(self, key):
            setattr(self, key, value)
            if key == 'ca' and self.lat == 'hcp':
                # c/a has changed => update lattice parameters
                self.latparams = [1.0, 1.0, self.ca]
                self.latvectors = [
                    [1.0, 0.0, 0.0], [-0.5, 0.8660254, 0.0], [0.0, 0.0, self.ca]]
                self.basis = [
                    [0.0, 0.0, 0.0], [0.0, 0.57735027, self.ca / 2.0]]
        else:
            print('WARNING: Kstr() class has no attribute \'{0}\''.format(key))
        return

    def check_input_file(self):
        """Perform various checks on the class data.

        Checking function to make sure that all necessary data exists
        before we attempt to write the input file to disk.

        :returns: None
        :rtype: None
        """

        if self.ca is None and self.lat == 'hcp':
            self.ca = 1.632993
            #sys.exit('Kstr.check_input_file: for hcp \'ca\' = c/a has to be given!')

        if self.lat is None:
            sys.exit('Kstr.check_input_file: \'lat\' has to be given!')
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
                self.basis = [
                    [0.0, 0.0, 0.0], [0.0, 0.57735027, self.ca / 2.0]]
            else:
                self.basis = [0.0, 0.0, 0.0]

        if isinstance(self.basis[0], list):
            self.nq = len(self.basis)
        else:
            self.nq = 1
            self.basis = [self.basis]
        if isinstance(self.latvectors[0], list):
            if len(self.latvectors) == 1:
                self.iprim = 1
            else:
                self.iprim = 0
        else:
            self.iprim = 1

        if self.kappaw is None:
            self.kappaw = [0.0]

        self.kappalen = len(self.kappaw)
        if self.kappalen == 2:
            self.twocenter = True
        else:
            self.twocenter = False
        if self.twocenter:
            self.jobname2 = self.jobname + '2'

        if self.dmax is None:
            if common.lat_to_ibz(self.lat) == 2:
                self.dmax = 1.7
            elif common.lat_to_ibz(self.lat) == 3:
                self.dmax = 2.2
            elif common.lat_to_ibz(self.lat) == 4:
                self.dmax = 2.4

        if self.msgl is None:
            self.msgl = 1
        if self.nprn is None:
            self.nprn = 0
        if self.lamda is None:
            self.lamda = 2.5
        if self.amax is None:
            self.amax = 4.5
        if self.bmax is None:
            self.bmax = 4.5
        if self.nqr2 is None:
            self.nqr2 = 0
        if self.mode is None:
            self.mode = 'B'
        if self.store is None:
            self.store = 'Y'
        if self.high is None:
            self.high = 'Y'
        if self.kstr_nl is None:
            self.kstr_nl = 4
        if self.nlh is None:
            self.nlh = 11
        if self.nlw is None:
            self.nlw = 9
        if self.nder is None:
            self.nder = 6
        if self.itrans is None:
            self.itrans = 3
        if self.rwats is None:
            self.rwats = 0.1
        if self.nghbp is None:
            self.nghbp = 13
        if self.awIQ is None:
            self.awIQ = np.ones((self.nq, 4)) * 0.7
        return

    def finalize(self):
        """Re-initializes input parameters."""

        self.awIQ = None
        return

