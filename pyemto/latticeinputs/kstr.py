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
                 rwats=None, nghbp=None, awIQ=None, numvec_target=None):

        # Argument checking and default values

        self.lat = lat
        self.jobname = jobname
        self.latparams = latparams
        self.latvectors = latvectors
        self.basis = np.asarray(basis)
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
        self.numvec_target = numvec_target

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
                .format(self.basis[i,0], self.basis[i,1], self.basis[i,2]) + "\n"
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
            fl = open(folder + '/{0}.kstr'.format(self.jobname2), "w")
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
                self.basis = np.asarray([
                    [0.0, 0.0, 0.0], [0.0, 0.57735027, self.ca / 2.0]])
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
                self.basis = np.asarray([
                    [0.0, 0.0, 0.0], [0.0, 0.57735027, self.ca / 2.0]])
            else:
                self.basis = np.asarray([0.0, 0.0, 0.0])

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

        if self.kappaw is None:
            self.kappaw = [0.0]

        self.kappalen = len(self.kappaw)

        if self.kappalen == 2:
            self.twocenter = True
        elif self.kappalen == 1:
            self.twocenter = False    
        else:
            print("Kappa does not have correct number of values, please correct! : %s" % (str(self.kappaw)))
            exit()
        if self.twocenter:
            self.jobname2 = self.jobname + 'M'

        # Optimize dmax
        if self.nghbp is None:
            self.nghbp = 13
        if self.numvec_target is None:
            self.numvec_target = 89
        # Prevent self.basis from being modified by the omtimize_dmax function
        # by using np.copy() function
        self.dmax,numvec_tmp = self.optimize_dmax(self.latvectors,np.copy(self.basis))
        #if self.dmax is None:
        #    if common.lat_to_ibz(self.lat) == 2:
        #        self.dmax = 1.7
        #    elif common.lat_to_ibz(self.lat) == 3:
        #        self.dmax = 2.2
        #    elif common.lat_to_ibz(self.lat) == 4:
        #        self.dmax = 2.4

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
        if self.awIQ is None:
            self.awIQ = np.ones((self.nq, 4)) * 0.7
        return


    def generate_prims(self,angles):
        """Generates the primitive lattice vectors from the angle and Bravais-lattice
        type information.
        """
        deg_to_rad = np.pi/180.0
        alpha = angles[0] * deg_to_rad
        beta  = angles[1] * deg_to_rad
        gamma = angles[2] * deg_to_rad
        ibz = common.lat_to_ibz(self.lat)
        boa = self.latparams[1]/self.latparams[0]
        coa = self.latparams[2]/self.latparams[0]
        #
        if ibz == 1:
            prims = np.array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]])
        elif ibz == 2:
            prims = np.array([[0.5,0.5,0.0],[0.0,0.5,0.5],[0.5,0.0,0.5]])
        elif ibz == 3:
            prims = np.array([[0.5,0.5,-0.5],[-0.5,0.5,0.5],[0.5,-0.5,0.5]])
        elif ibz == 4:
            prims = np.array([[1.0,0.0,0.0],[-0.5,0.5*np.sqrt(3.0),0.0],[0.0,0.0,coa]])
        elif ibz == 5:
            prims = np.array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,coa]])
        elif ibz == 6:
            prims = np.array([[-0.5,0.5,0.5*coa],[0.5,-0.5,0.5*coa],[0.5,0.5,-0.5*coa]])
        elif ibz == 7:
            prims = np.array([[0.0,1.0,coa],[-0.5*np.sqrt(3.0),-0.5,coa],[0.5*np.sqrt(3.0),-0.5,coa]])
        elif ibz == 8:
            prims = np.array([[1.0,0.0,0.0],[0.0,boa,0.0],[0.0,0.0,coa]])
        elif ibz == 9:
            prims = np.array([[0.5,-0.5*boa,0.0],[0.5,0.5*boa,0.0],[0.0,0.0,coa]])
        elif ibz == 10:
            prims = np.array([[0.5,-0.5*boa,0.5*coa],[0.5,0.5*boa,-0.5*coa],[-0.5,0.5*boa,0.5*coa]])
        elif ibz == 11:
            prims = np.array([[0.5,0.0,0.5*coa],[0.5,0.5*boa,0.0],[0.0,0.5*boa,0.5*coa]])
        elif ibz == 12:
            prims = np.array([[1.0,0.0,0.0],[boa*np.cos(gamma),boa*np.sin(gamma),0.0],
                              [0.0,0.0,coa]])
        elif ibz == 13:
            prims = np.array([[0.0,-boa,0.0],[0.5*np.sin(gamma),-0.5*np.cos(gamma),-0.5*coa],
                              [0.5*np.sin(gamma),-0.5*np.cos(gamma),0.5*coa]])
        elif ibz == 14:
            tmp1 = (np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma)
            tmp2 = np.sqrt((1.0-np.cos(gamma)**2-np.cos(alpha)**2-np.cos(beta)**2
                            +2.0*np.cos(alpha)*np.cos(beta)*np.cos(gamma)))/np.sin(gamma)
            prims = np.array([[1.0,0.0,0.0],[boa*np.cos(gamma),boa*np.sin(gamma),0.0],
                              [coa*np.cos(beta),coa*tmp1,coa*tmp2]])
        #print('generated prims:')
        #print(prims)
        return prims

    def compute_num_of_vecs(self,prims,basis,nghbp,dmax):
        """Computes the number of vectors for a given dmax value."""
        ncrq = 0
        nq = self.nq
        #
        qix = basis[0,0]
        qiy = basis[0,1]
        qiz = basis[0,2]
        #
        for jq in range(nq):
            qmqpx=basis[jq,0]-qix
            qmqpy=basis[jq,1]-qiy
            qmqpz=basis[jq,2]-qiz
            for l in range(-nghbp,nghbp+1):
                for m in range(-nghbp,nghbp+1):
                    for n in range(-nghbp,nghbp+1):
                        sx=qmqpx+l*prims[0,0]+m*prims[1,0]+n*prims[2,0]
                        sy=qmqpy+l*prims[0,1]+m*prims[1,1]+n*prims[2,1]
                        sz=qmqpz+l*prims[0,2]+m*prims[1,2]+n*prims[2,2]
                        dx=np.sqrt(sx*sx+sy*sy+sz*sz)
                        if dx <= dmax:
                            ncrq += 1
                            #print('DX,NCRQ = ',dx,ncrq)
        return ncrq

    def optimize_dmax(self,prims,basis):
        """Calculates the best possible dmax value, which gives the closest number
        of vectors given some target value (which is typically 80-90)."""
        # Check if primitive vectors have been
        # explicitly given (iprim = 0).
        # If not (iprim = 1), we have to generate
        # from the alpha, beta, and gamma + ibz
        # information
        #print('basis:')
        #print(basis)
        if self.iprim == 0:
            prims = np.asarray(prims)
            basis = np.asarray(basis)
            # Take the first site as the origin
            basis[:,0] -= basis[0,0]
            basis[:,1] -= basis[0,1]
            basis[:,2] -= basis[0,2]
        elif self.iprim == 1:
            prims = self.generate_prims(prims)
            basis = np.asarray(basis)
            # Take the first site as the origin
            basis[:,0] -= basis[0,0]
            basis[:,1] -= basis[0,1]
            basis[:,2] -= basis[0,2]
        #print(prims)
        #print(basis)
        #
        nghbp = self.nghbp
        dmax_min = 0.5
        dmax_max = 8.0
        ncrq_target = self.numvec_target
        tol_target = 1.0e-4
        # Initial dmax guess
        dmax_mid_old = 1000.0
        dmax_mid = np.round((dmax_max + dmax_min)/2,4)
        f_closest = 1000
        #
        print('Latticeinputs.Kstr.optimize_dmax(): Optimizing dmax (target = {0:3d}) for {1}...'
              .format(ncrq_target,self.jobname))
        #
        while np.abs(dmax_mid_old - dmax_mid) > tol_target:
            # Compute the number of vectors that
            # corresponds to the trial choice of dmax
            # ncrq = Number of vectors for the first site
            # in the basis, which has been shifted to lie
            # in the origin.
            f_mid = self.compute_num_of_vecs(prims,basis,nghbp,dmax_mid)
            #
            if f_mid == ncrq_target:
                # Target number of vectors has been reached exactly
                dmax_closest = dmax_mid
                f_closest = f_mid
                break
            elif f_mid - ncrq_target < 0:
                # Number of vectors too small
                if np.abs(f_mid-ncrq_target) < np.abs(f_closest-ncrq_target):
                    dmax_closest = dmax_mid
                    f_closest = f_mid
                dmax_min = dmax_mid
            else:
                # Number of vectors too large.
                if np.abs(f_mid-ncrq_target) < np.abs(f_closest-ncrq_target):
                    dmax_closest = dmax_mid
                    f_closest = f_mid
                elif np.abs(f_mid-ncrq_target) == np.abs(f_closest-ncrq_target):
                    if f_mid > f_closest:
                        dmax_closest = dmax_mid
                        f_closest = f_mid
                dmax_max = dmax_mid
            #
            dmax_mid_old = dmax_mid
            dmax_mid = np.round((dmax_max + dmax_min)/2,4)
            #print('dmax, Number of vectors: {0:6.4f}, {1}'.format(dmax_mid_old,f_mid))
        print('dmax, Number of vectors: {0:6.4f}, {1}'.format(dmax_closest,f_closest))
        return dmax_closest,f_closest

    def finalize(self):
        """Re-initializes input parameters."""

        self.awIQ = None
        #self.dmax = None
        #self.kappaw = None
        #self.jobname = None
        #self.latparams = None
        #self.latvectors = None
        #self.basis = None
        return

