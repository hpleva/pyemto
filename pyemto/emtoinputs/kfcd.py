# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 15:05:43 2014

@author: Matti Ropo
@author: Henrik Levämäki

"""

from __future__ import print_function
import sys
import os
import datetime
import re
import pyemto.common.common as common

class Kfcd:

    """Handles the information and writing of kfcd file.

    :param jobname:  (Default value = None)
    :type jobname:
    :param latname:  (Default value = None)
    :type latname:
    :param latpath:  (Default value = None)
    :type latpath:
    :param msgl:  (Default value = None)
    :type msgl:
    :param nprn:  (Default value = None)
    :type nprn:
    :param lmaxs:  (Default value = None)
    :type lmaxs:
    :param nth:  (Default value = None)
    :type nth:
    :param kfcd_nfi:  (Default value = None)
    :type kfcd_nfi:
    :param fpot:  (Default value = None)
    :type fpot:
    :param ovcor:  (Default value = None)
    :type ovcor:
    :param ubg:  (Default value = None)
    :type ubg:
    :param DIR001:  (Default value = None)
    :type DIR001:
    :param DIR002:  (Default value = None)
    :type DIR002:
    :param DIR003:  (Default value = None)
    :type DIR003:
    :param DIR004:  (Default value = None)
    :type DIR004:
    :param DIR006:  (Default value = None)
    :type DIR006:
    :param sws:  (Default value = None)
    :type sws:
    :returns: None
    :rtype: None
    """

    def __init__(self, jobname=None, latname=None, latpath=None, msgl=None, nprn=None,
                 lmaxs=None, nth=None, kfcd_nfi=None, fpot=None, ovcor=None,
                 ubg=None, DIR001=None, DIR002=None, DIR003=None, DIR004=None,
                 DIR006=None, sws=None, CQNA=None, KFCD_file_type=None):

        self.jobname = jobname
        self.latname = latname
        self.latpath = latpath
        self.msgl = msgl
        self.nprn = nprn
        self.lmaxs = lmaxs
        self.nth = nth
        self.kfcd_nfi = kfcd_nfi
        self.fpot = fpot
        self.ovcor = ovcor
        self.ubg = ubg
        self.DIR001 = DIR001
        self.DIR002 = DIR002
        self.DIR003 = DIR003
        self.DIR004 = DIR004
        self.DIR006 = DIR006
        self.CQNA = CQNA
        self.KFCD_file_type = KFCD_file_type
        if KFCD_file_type is None:
            self.KFCD_file_type = 'kfcd'
        
    def output(self):
        """Outputs KFCD input file as a formatted string.

        Outputs EMTO5.8 KFCD input file

        :returns: Formatted string
        :rtype: str
        """

        now = datetime.datetime.now()
        #line = "KFCD      MSGL..=  %1i                        " % (self.msgl) +\        
        line = "KFCD      MSGL..=  {0:1}".format(self.msgl) + " CQNA=  {0}      ".format(self.CQNA) +\
               str(now.day) + "." + str(now.month) + "." + str(now.year) + "\n"
        line = line + "JOBNAM...=" + self.jobname + "\n"
        line = line + "STRNAM...=" + self.latname + "\n"
        line = line + "DIR001=" + self.DIR001 + "\n"
        line = line + "DIR002=" + self.DIR002 + "\n"
        line = line + "DIR003=" + self.DIR003 + "\n"
        line = line + "DIR004=" + self.DIR004 + "\n"
        line = line + "DIR006=" + self.DIR006 + "\n"
        line = line + "Lmaxs.=%3i NTH..=%3i NFI..=%3i FPOT..= %1s"\
            % (self.lmaxs, self.nth, self.kfcd_nfi, self.fpot) + "\n"
        line = line + "OVCOR.=  %1s UBG..=  %1s NPRN.=  %1s"\
            % (self.ovcor, self.ubg, self.nprn) + "\n"

        return line

    def write_input_file(self, folder=None):
        """Writes input file to disk.

        Save KFCD input data to file named filename.

        :param folder: Folder where the data will be written (Default value = None)
        :type folder: str
        :returns: None
        :rtype: None
        """

        # Check data integrity before anything is written on disk or run
        self.check_input_file()

        if folder is None:
            #sys.exit('Kfcd.create_input_file: \'folder\' has to be given!')
            folder = './'
        else:
            common.check_folders(folder)

        fl = open(folder + '/{0}.{1}'.format(self.jobname,self.KFCD_file_type), "w")
        fl.write(self.output())
        fl.close()
        return

    def set_values(self, key, value):
        """Changes values of the class variables.

        :param key: name of the variable
        :type key: str
        :param value: value of the variable
        :type value: str, int or float
        :returns: None
        :rtype: None
        """
        if not hasattr(self, key):
            print('WARNING: Kfcd() class has no attribute \'{0}\''.format(key))
            return

        # Lattice name or path has changed => we have to update the FOR and DIR
        # information
        elif key == 'latpath':
            setattr(self, key, value)
            self.DIR001 = self.latpath + '/kstr/'
            self.DIR001 = common.cleanup_path(self.DIR001)
            self.DIR003 = self.latpath + '/shape/'
            self.DIR003 = common.cleanup_path(self.DIR003)
            self.DIR004 = self.latpath + '/bmdl/'
            self.DIR004 = common.cleanup_path(self.DIR004)

        else:
            setattr(self, key, value)
        return

    def check_input_file(self):
        """Perform various checks on the class data.

        Makes sure that all necessary data exists
        before we attempt to write the input file to disk

        :returns: None
        :rtype: None
        """

        if self.jobname is None:
            sys.exit('Kfcd: \'jobname\' has to be given!')
        if self.latname is None:
            sys.exit('Kfcd: \'latname\' has to be given!')
        if self.latpath is None:
            self.latpath = './'
        if self.msgl is None:
            self.msgl = 0
        if self.nprn is None:
            self.nprn = 0
        if self.lmaxs is None:
            self.lmaxs = 30
        if self.nth is None:
            self.nth = 41
        if self.kfcd_nfi is None:
            self.kfcd_nfi = 81
        if self.fpot is None:
            self.fpot = 'N'
        if self.ovcor is None:
            self.ovcor = 'Y'
        if self.ubg is None:
            self.ubg = 'N'
        if self.DIR001 is None:
            self.DIR001 = self.latpath + '/kstr/'
            self.DIR001 = common.cleanup_path(self.DIR001)
        if self.DIR002 is None:
            self.DIR002 = 'kgrn/'
            self.DIR002 = common.cleanup_path(self.DIR002)
        if self.DIR003 is None:
            self.DIR003 = self.latpath + '/shape/'
            self.DIR003 = common.cleanup_path(self.DIR003)
        if self.DIR004 is None:
            self.DIR004 = self.latpath + '/bmdl/'
            self.DIR004 = common.cleanup_path(self.DIR004)
        if self.DIR006 is None:
            self.DIR006 = 'kfcd/'
            self.DIR006 = common.cleanup_path(self.DIR006)
        if self.CQNA is None:
            self.CQNA = 'N'
        return

