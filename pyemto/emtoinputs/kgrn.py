# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 14:48:25 2014

@author: Matti Ropo
@author: Henrik Levämäki

"""

import sys
import os
import datetime
import numpy as np
import re
import pyemto.common.common as common

class Kgrn:

    """A class which contains all the KGRN input file related information.

    :param jobname:  (Default value = None)
    :type jobname:
    :param latname:  (Default value = None)
    :type latname:
    :param latpath:  (Default value = None)
    :type latpath:
    :param ibz:  (Default value = None)
    :type ibz:
    :param atoms:  (Default value = None)
    :type atoms:
    :param concs:  (Default value = None)
    :type concs:
    :param iqs:  (Default value = None)
    :type iqs:
    :param its:  (Default value = None)
    :type its:
    :param itas:  (Default value = None)
    :type itas:
    :param qtrs:  (Default value = None)
    :type qtrs:
    :param splts:  (Default value = None)
    :type splts:
    :param fixs:  (Default value = None)
    :type fixs:
    :param sm_ss:  (Default value = None)
    :type sm_ss:
    :param s_wss:  (Default value = None)
    :type s_wss:
    :param ws_wsts:  (Default value = None)
    :type ws_wsts:
    :param atconf:  (Default value = None)
    :type atconf:
    :param sws:  (Default value = None)
    :type sws:
    :param strt:  (Default value = None)
    :type strt:
    :param msgl:  (Default value = None)
    :type msgl:
    :param expan:  (Default value = None)
    :type expan:
    :param fcd:  (Default value = None)
    :type fcd:
    :param func:  (Default value = None)
    :type func:
    :param niter:  (Default value = None)
    :type niter:
    :param nlin:  (Default value = None)
    :type nlin:
    :param nprn:  (Default value = None)
    :type nprn:
    :param ncpa:  (Default value = None)
    :type ncpa:
    :param mode:  (Default value = None)
    :type mode:
    :param frc:  (Default value = None)
    :type frc:
    :param dos:  (Default value = None)
    :type dos:
    :param ops:  (Default value = None)
    :type ops:
    :param afm:  (Default value = None)
    :type afm:
    :param crt:  (Default value = None)
    :type crt:
    :param lmaxh:  (Default value = None)
    :type lmaxh:
    :param lmaxt:  (Default value = None)
    :type lmaxt:
    :param kgrn_nfi:  (Default value = None)
    :type kgrn_nfi:
    :param fixg:  (Default value = None)
    :type fixg:
    :param shf:  (Default value = None)
    :type shf:
    :param sofc:  (Default value = None)
    :type sofc:
    :param kmsh:  (Default value = None)
    :type kmsh:
    :param nkx:  (Default value = None)
    :type nkx:
    :param nky:  (Default value = None)
    :type nky:
    :param nkz:  (Default value = None)
    :type nkz:
    :param fbz:  (Default value = None)
    :type fbz:
    :param kmsh2:  (Default value = None)
    :type kmsh2:
    :param ibz2:  (Default value = None)
    :type ibz2:
    :param nkx2:  (Default value = None)
    :type nkx2:
    :param nky2:  (Default value = None)
    :type nky2:
    :param nkz2:  (Default value = None)
    :type nkz2:
    :param zmsh:  (Default value = None)
    :type zmsh:
    :param nz1:  (Default value = None)
    :type nz1:
    :param nz2:  (Default value = None)
    :type nz2:
    :param nz3:  (Default value = None)
    :type nz3:
    :param nres:  (Default value = None)
    :type nres:
    :param nzd:  (Default value = None)
    :type nzd:
    :param depth:  (Default value = None)
    :type depth:
    :param imagz:  (Default value = None)
    :type imagz:
    :param eps:  (Default value = None)
    :type eps:
    :param elim:  (Default value = None)
    :type elim:
    :param amix:  (Default value = None)
    :type amix:
    :param efmix:  (Default value = None)
    :type efmix:
    :param vmtz:  (Default value = None)
    :type vmtz:
    :param mmom:  (Default value = None)
    :type mmom:
    :param tole:  (Default value = None)
    :type tole:
    :param tolef:  (Default value = None)
    :type tolef:
    :param tolcpa:  (Default value = None)
    :type tolcpa:
    :param tfermi:  (Default value = None)
    :type tfermi:
    :param nsws:  (Default value = None)
    :type nsws:
    :param dsws:  (Default value = None)
    :type dsws:
    :param alpcpa:  (Default value = None)
    :type alpcpa:
    :param efgs:  (Default value = None)
    :type efgs:
    :param hx:  (Default value = None)
    :type hx:
    :param nx:  (Default value = None)
    :type nx:
    :param nz0:  (Default value = None)
    :type nz0:
    :param stmp:  (Default value = None)
    :type stmp:
    :param iex:  (Default value = None)
    :type iex:
    :param dirac_np:  (Default value = None)
    :type dirac_np:
    :param nes:  (Default value = None)
    :type nes:
    :param dirac_niter:  (Default value = None)
    :type dirac_niter:
    :param iwat:  (Default value = None)
    :type iwat:
    :param nprna:  (Default value = None)
    :type nprna:
    :param vmix:  (Default value = None)
    :type vmix:
    :param rwat:  (Default value = None)
    :type rwat:
    :param rmax:  (Default value = None)
    :type rmax:
    :param dx:  (Default value = None)
    :type dx:
    :param dr1:  (Default value = None)
    :type dr1:
    :param test:  (Default value = None)
    :type test:
    :param teste:  (Default value = None)
    :type teste:
    :param testy:  (Default value = None)
    :type testy:
    :param testv:  (Default value = None)
    :type testv:
    :param FOR001:  (Default value = None)
    :type FOR001:
    :param DIR002:  (Default value = None)
    :type DIR002:
    :param DIR003:  (Default value = None)
    :type DIR003:
    :param FOR004:  (Default value = None)
    :type FOR004:
    :param DIR006:  (Default value = None)
    :type DIR006:
    :param DIR009:  (Default value = None)
    :type DIR009:
    :param DIR010:  (Default value = None)
    :type DIR010:
    :param DIR011:  (Default value = None)
    :type DIR011:
    :returns: None
    :rtype: None
    """

    def __init__(self, jobname=None, latname=None, latpath=None, ibz=None, atoms=None,
                 concs=None, iqs=None, its=None, itas=None, qtrs=None, splts=None, fixs=None,
                 sm_ss=None, s_wss=None, ws_wsts=None,
                 atconf=None, sws=None, strt=None, msgl=None, expan=None, fcd=None, func=None,
                 niter=None, nlin=None, nprn=None, ncpa=None, mode=None, frc=None,
                 dos=None, ops=None, afm=None, crt=None, lmaxh=None, lmaxt=None,
                 kgrn_nfi=None, fixg=None, shf=None, sofc=None, kmsh=None, nkx=None,
                 nky=None, nkz=None, fbz=None, kmsh2=None, ibz2=None, nkx2=None,
                 nky2=None, nkz2=None, zmsh=None, nz1=None, nz2=None, nz3=None,
                 nres=None, nzd=None, depth=None, imagz=None, eps=None, elim=None,
                 amix=None, efmix=None, vmtz=None, mmom=None, tole=None, tolef=None,
                 tolcpa=None, tfermi=None, nsws=None, dsws=None, alpcpa=None, efgs=None,
                 hx=None, nx=None, nz0=None, stmp=None, iex=None, dirac_np=None, nes=None,
                 dirac_niter=None, iwat=None, nprna=None, vmix=None, rwat=None, rmax=None,
                 dx=None, dr1=None, test=None, teste=None, testy=None, testv=None,
                 FOR001=None, DIR002=None, DIR003=None, FOR004=None, DIR006=None, DIR009=None,
                 DIR010=None, DIR011=None):

        self.jobname = jobname
        self.latname = latname
        self.ibz = ibz
        self.latpath = latpath
        self.atoms = atoms
        self.concs = concs
        self.iqs = iqs
        self.its = its
        self.itas = itas
        self.qtrs = qtrs
        self.splts = splts
        self.fixs = fixs
        self.sm_ss = sm_ss
        self.s_wss = s_wss
        self.ws_wsts = ws_wsts
        self.atconf = atconf
        self.sws = sws
        self.strt = strt
        self.msgl = msgl
        self.expan = expan
        self.fcd = fcd
        self.func = func
        self.niter = niter
        self.nlin = nlin
        self.nprn = nprn
        self.ncpa = ncpa
        self.mode = mode
        self.frc = frc
        self.dos = dos
        self.ops = ops
        self.afm = afm
        self.crt = crt
        self.lmaxh = lmaxh
        self.lmaxt = lmaxt
        self.kgrn_nfi = kgrn_nfi
        self.fixg = fixg
        self.shf = shf
        self.sofc = sofc
        self.kmsh = kmsh
        self.nkx = nkx
        self.nky = nky
        self.nkz = nkz
        self.fbz = fbz
        self.kmsh2 = kmsh2
        self.ibz2 = ibz2
        self.nkx2 = nkx2
        self.nky2 = nky2
        self.nkz2 = nkz2
        self.zmsh = zmsh
        self.nz1 = nz1
        self.nz2 = nz2
        self.nz3 = nz3
        self.nres = nres
        self.nzd = nzd
        self.depth = depth
        self.imagz = imagz
        self.eps = eps
        self.elim = elim
        self.amix = amix
        self.efmix = efmix
        self.vmtz = vmtz
        self.mmom = mmom
        self.tole = tole
        self.tolef = tolef
        self.tolcpa = tolcpa
        self.tfermi = tfermi
        self.nsws = nsws
        self.dsws = dsws
        self.alpcpa = alpcpa
        self.efgs = efgs
        self.hx = hx
        self.nx = nx
        self.nz0 = nz0
        self.stmp = stmp
        self.iex = iex
        self.dirac_np = dirac_np
        self.nes = nes
        self.dirac_niter = dirac_niter
        self.iwat = iwat
        self.nprna = nprna
        self.vmix = vmix
        self.rwat = rwat
        self.rmax = rmax
        self.dx = dx
        self.dr1 = dr1
        self.test = test
        self.teste = teste
        self.testy = testy
        self.testv = testv
        self.FOR001 = FOR001
        self.DIR002 = DIR002
        self.DIR003 = DIR003
        self.FOR004 = FOR004
        self.DIR006 = DIR006
        self.DIR009 = DIR009
        self.DIR010 = DIR010
        self.DIR011 = DIR011

        if iqs is None:
            self.iqs = np.ones(50, dtype='int32')
        else:
            self.iqs = iqs
        if its is None:
            self.its = np.ones(50, dtype='int32')
        else:
            self.its = its
        if itas is None:
            self.itas = np.arange(1, 50, dtype='int32')
        else:
            self.itas = itas
        if qtrs is None:
            self.qtrs = np.zeros(50)
        else:
            self.qtrs = qtrs
        if splts is None:
            self.splts = np.zeros(50)
        else:
            self.splts = splts
        if fixs is None:
            self.fixs = np.empty([50], dtype='S1')
            self.fixs[:] = 'N'
        else:
            self.fixs = fixs
        if sm_ss is None:
            self.sm_ss = np.ones(50)
        else:
            self.sm_ss = sm_ss
        if s_wss is None:
            self.s_wss = np.ones(50)
        else:
            self.s_wss = s_wss
        if ws_wsts is None:
            self.ws_wsts = np.ones(50)
        else:
            self.ws_wsts = ws_wsts

        self.atconfKeys = ['atoms', 'concs', 'iqs', 'its', 'itas', 'qtrs', 'splts',
                           'fixs', 'sm_ss', 's_wss', 'ws_wsts']

        self.symtonum = {"Va": 0, "Em": 0, "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5,
                         "C": 6, "N": 7, "O": 8, "O-2": 8, "F": 9, "Ne": 10, "Na": 11,
                         "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18,
                         "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25,
                         "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30, "Ga": 31, "Ge": 32,
                         "As": 33, "Se": 34, "Br": 35, "Kr": 36, "Rb": 37, "Sr": 38, "Y": 39,
                         "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46,
                         "Ag": 47, "Cd": 48, "In": 49, "Sn": 50, "Sb": 51, "Te": 52, "I": 53,
                         "Xe": 54, "Cs": 55, "Ba": 56, "La": 57, "Ce": 58, "Dy": 66, "Tm": 69,
                         "Lu": 71, "Hf": 72, "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77,
                         "Pt": 78, "Au": 79, "Hg": 80, "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84,
                         "At": 85, "Rn": 86, "Fr": 87, "Ra": 88, "Ac": 89, "Th": 90, "Pa": 91,
                         "U": 92, "Np": 93, "Pu": 94, "Am": 95, "Cm": 96, "Bk": 97, "Cf": 98}

        self.atomblock = None

        self.Confdtype = [('elem', 'U4'), ('iq', int), ('it', int), ('ita', int),
                          ('conc', float), ('Sm_s', float), ('S_ws',
                                                             float), ('WS_wst', float),
                          ('qtr', float), ('mmom', float), ('fix', 'U1')]

        # Try creating the atomblock
        self.create_atomblock()

    def Atomline(
            self, atom, iq, it, ita, conc, Sm_s, S_ws, WS_wst, qtr, splt, fix):
        """Prints atomic line in kgrn format

        :param atom:
        :type atom:
        :param iq:
        :type iq:
        :param it:
        :type it:
        :param ita:
        :type ita:
        :param conc:
        :type conc:
        :param Sm_s:
        :type Sm_s:
        :param S_ws:
        :type S_ws:
        :param WS_wst:
        :type WS_wst:
        :param qtr:
        :type qtr:
        :param splt:
        :type splt:
        :param fix:
        :type fix:
        :returns:
        :rtype:
        """

        line = "%-2s    %3i %2i %2i  %2i  %5.3f  %5.3f  %5.3f  %5.3f  %3.1f %4.1f  %s" \
            % (atom, iq, it, ita, self.symtonum[atom], conc, Sm_s, S_ws, WS_wst, qtr, splt, fix)

        return line

    def AtomOutput(self):
        """(self) -> str

            Output of atomic lines in kgrn format

        :returns:
        :rtype:
        """

        lines = ""
        for i in range(len(self.atomblock)):
            lines = lines + self.atomblock[i] + "\n"
        return lines

    def aconflines(self, atype):
        """Returns string containing electronic configuration lines

        Contains the database of atomic orbital configuration information
        which is present in atomic block of the KGRN input file.

        :param atype: Name of the element in the periodic table
        :type atype: str
        :returns: A formatted string corresponding to the parameter atype
        :rtype: str
        """

        if atype == "Va":
            results = "Va\n" \
                + "Iz=   0 Norb=  0 Ion=  0 Config= 1s0\n" \
                + "n      1\n" \
                + "Kappa -1\n" \
                + "Occup  0\n" \
                + "Valen  1\n"
        elif atype == "Em":
            results = "Em\n" \
                + "Iz=   0 Norb=  0 Ion=  0 Config= 1s0\n" \
                + "n      1\n" \
                + "Kappa -1\n" \
                + "Occup  0\n" \
                + "Valen  1\n"
        elif atype == "H":
            results = "H\n" \
                + "Iz=   1 Norb=  1 Ion=  0 Config= 1s1\n" \
                + "n      1\n" \
                + "Kappa -1\n" \
                + "Occup  1\n" \
                + "Valen  1\n"
        elif atype == "He":
            results = "He\n" \
                + "Iz=   2 Norb=  1 Ion=  0 Config= 1s2\n" \
                + "n      1\n" \
                + "Kappa -1\n" \
                + "Occup  2\n" \
                + "Valen  1\n"
        elif atype == "Li":
            results = "Li\n" \
                + "Iz=   3 Norb=  2 Ion=  0 Config= 2s1\n" \
                + "n      1  2\n" \
                + "Kappa -1 -1\n" \
                + "Occup  2  1\n" \
                + "Valen  0  1\n"
        elif atype == "Be":
            results = "Be\n" \
                + "Iz=   4 Norb=  2 Ion=  0 Config= 2s2\n" \
                + "n      1  2\n" \
                + "Kappa -1 -1\n" \
                + "Occup  2  2\n" \
                + "Valen  0  1\n"
        elif atype == "B":
            results = "B\n" \
                + "Iz=   5 Norb=  3 Ion=  0 Config= 2s1_2p1\n" \
                + "n      1  2  2\n" \
                + "Kappa -1 -1  1\n" \
                + "Occup  2  2  1\n" \
                + "Valen  0  1  1\n"
        elif atype == "C":
            results = "C\n" \
                + "Iz=   6 Norb=  3 Ion=  0 Config= 2s2 2p2\n" \
                + "n      1  2  2\n" \
                + "Kappa -1 -1  1\n" \
                + "Occup  2  2  2\n" \
                + "Valen  0  1  1\n"
        elif atype == "N":
            results = "N\n" \
                + "Iz=   7 Norb=  4 Ion=  0 Config= 2p3\n" \
                + "n      1  2  2  2\n" \
                + "Kappa -1 -1  1 -2\n" \
                + "Occup  2  2  2  1\n" \
                + "Valen  0  0  1  1\n"
        elif atype == "O":
            results = "O\n" \
                + "Iz=   8 Norb=  4 Ion=  0 Config= 2s2_2p4\n" \
                + "n      1  2  2  2\n" \
                + "Kappa -1 -1  1 -2\n" \
                + "Occup  2  2  2  2\n" \
                + "Valen  0  1  1  1\n"
        elif atype == "O-2":
            results = "O-2\n" \
                + "Iz=   8 Norb=  4 Ion= -2 Config= 2p6\n" \
                + "n      1  2  2  2\n" \
                + "Kappa -1 -1  1 -2\n" \
                + "Occup  2  2  2  4\n" \
                + "Valen  0  0  1  1\n"
        elif atype == "F":
            results = "F\n" \
                + "Iz=   9 Norb=  4 Ion=  0 Config= 2s2_2p5\n" \
                + "n      1  2  2  2\n" \
                + "Kappa -1 -1  1 -2\n" \
                + "Occup  2  2  2  3\n" \
                + "Valen  0  1  1  1\n"
        elif atype == "Ne":
            results = "Ne\n" \
                + "Iz=  10 Norb=  4 Ion=  0 Config= 2s2_2p6\n" \
                + "n      1  2  2  2\n" \
                + "Kappa -1 -1  1 -2\n" \
                + "Occup  2  2  2  4\n" \
                + "Valen  0  1  1  1\n"
        elif atype == "Na":
            results = "Na\n" \
                + "Iz=  11 Norb=  5 Ion=  0 Config= 3s1\n" \
                + "n      1  2  2  2  3\n" \
                + "Kappa -1 -1  1 -2 -1\n" \
                + "Occup  2  2  2  4  1\n" \
                + "Valen  0  0  0  0  1\n"
        elif atype == "Mg":
            results = "Mg\n" \
                + "Iz=  12 Norb=  5 Ion=  0 Config= 3s2\n" \
                + "n      1  2  2  2  3\n" \
                + "Kappa -1 -1  1 -2 -1\n" \
                + "Occup  2  2  2  4  2\n" \
                + "Valen  0  0  0  0  1\n"
        elif atype == "Al":
            results = "Al\n" \
                + "Iz=  13 Norb=  6 Ion=  0 Config= 3s2_3p1\n" \
                + "n      1  2  2  2  3  3\n" \
                + "Kappa -1 -1  1 -2 -1  1\n" \
                + "Occup  2  2  2  4  2  1\n" \
                + "Valen  0  0  0  0  1  1\n"
        elif atype == "Si":
            results = "Si\n" \
                + "Iz=  14 Norb=  6 Ion=  0 Config= 3s2_3p2\n" \
                + "n      1  2  2  2  3  3\n" \
                + "Kappa -1 -1  1 -2 -1  1\n" \
                + "Occup  2  2  2  4  2  2\n" \
                + "Valen  0  0  0  0  1  1\n"
        elif atype == "P":
            results = "P\n" \
                + "Iz=  15 Norb=  7 Ion=  0 Config= 3s2_3p3\n" \
                + "n      1  2  2  2  3  3  3\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2\n" \
                + "Occup  2  2  2  4  2  2  1\n" \
                + "Valen  0  0  0  0  1  1  1\n"
        elif atype == "S":
            results = "S\n" \
                + "Iz=  16 Norb=  7 Ion=  0 Config= 3s2_3p4\n" \
                + "n      1  2  2  2  3  3  3\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2\n" \
                + "Occup  2  2  2  4  2  2  2\n" \
                + "Valen  0  0  0  0  1  1  1\n"
        elif atype == "Cl":
            results = "Cl\n" \
                + "Iz=  17 Norb=  7 Ion=  0 Config= 3s2_3p5\n" \
                + "n      1  2  2  2  3  3  3\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2\n" \
                + "Occup  2  2  2  4  2  2  3\n" \
                + "Valen  0  0  0  0  1  1  1\n"
        elif atype == "Ar":
            results = "Ar\n" \
                + "Iz=  18 Norb=  7 Ion=  0 Config= 3s2_3p6\n" \
                + "n      1  2  2  2  3  3  3\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2\n" \
                + "Occup  2  2  2  4  2  2  4\n" \
                + "Valen  0  0  0  0  1  1  1\n"
        elif atype == "K":
            results = "K\n" \
                + "Iz=  19 Norb=  8 Ion=  0 Config 4s1\n" \
                + "n      1  2  2  2  3  3  3  4\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2 -1\n" \
                + "Occup  2  2  2  4  2  2  4  1\n" \
                + "Valen  0  0  0  0  0  0  0  1\n"
        elif atype == "Ca":
            results = "Ca\n" \
                    + "Iz=  20 Norb=  8 Ion=  0 Config= 4s2\n" \
                    + "n      1  2  2  2  3  3  3  4\n" \
                    + "Kappa -1 -1  1 -2 -1  1 -2 -1\n" \
                    + "Occup  2  2  2  4  2  2  4  2\n" \
                    + "Valen  0  0  0  0  0  0  0  1\n"
        #elif atype == "Ca":
        #    results = "Ca\n" \
        #            + "Iz=  20 Norb=  8 Ion=  0 Config= 3s2_3p6_4s2\n" \
        #            + "n      1  2  2  2  3  3  3  4\n" \
        #            + "Kappa -1 -1  1 -2 -1  1 -2 -1\n" \
        #            + "Occup  2  2  2  4  2  2  4  2\n" \
        #            + "Valen  0  0  0  0  1  1  1  1\n"
        elif atype == "Sc":
            results = "Sc\n" \
                + "Iz=  21 Norb=  9 Ion=  0 Config= 3d1_4s2\n" \
                + "n      1  2  2  2  3  3  3  3  4\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -1\n" \
                + "Occup  2  2  2  4  2  2  4  1  2\n" \
                + "Valen  0  0  0  0  0  0  0  1  1\n"
        elif atype == "Ti":
            results = "Ti\n" \
                + "Iz=  22 Norb=  9 Ion=  0 Config= 3d2_4s2\n" \
                + "n      1  2  2  2  3  3  3  3  4\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -1\n" \
                + "Occup  2  2  2  4  2  2  4  2  2\n" \
                + "Valen  0  0  0  0  0  0  0  1  1\n"
        elif atype == "V":
            results = "V\n" \
                + "Iz=  23 Norb=  9 Ion=  0 Config= 3d3_4s2\n" \
                + "n      1  2  2  2  3  3  3  3  4\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -1\n" \
                + "Occup  2  2  2  4  2  2  4  3  2\n" \
                + "Valen  0  0  0  0  0  0  0  1  1\n"
        elif atype == "Cr":
            results = "Cr\n" \
                + "Iz=  24 Norb=  9 Ion=  0 Config= 3d4_4s2\n" \
                + "n      1  2  2  2  3  3  3  3  4\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  2\n" \
                + "Valen  0  0  0  0  0  0  0  1  1\n"
        elif atype == "Mn":
            results = "Mn\n" \
                + "Iz=  25 Norb= 10 Ion=  0 Config= 3d5_4s2\n" \
                + "n      1  2  2  2  3  3  3  3  3  4\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  1  2\n" \
                + "Valen  0  0  0  0  0  0  0  1  1  1\n"
        elif atype == "Fe":
            results = "Fe\n" \
                + "Iz=  26 Norb= 10 Ion=  0 Config= 3d7_4s1\n" \
                + "n      1  2  2  2  3  3  3  3  3  4\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  3  1\n" \
                + "Valen  0  0  0  0  0  0  0  1  1  1\n"
        elif atype == "Co":
            results = "Co\n" \
                + "Iz=  27 Norb= 10 Ion=  0 Config= 3d7_4s2\n" \
                + "n      1  2  2  2  3  3  3  3  3  4\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  3  2\n" \
                + "Valen  0  0  0  0  0  0  0  1  1  1\n"
        elif atype == "Ni":
            results = "Ni\n" \
                    + "Iz=  28 Norb= 10 Ion=  0 Config= 3d8_4s2\n" \
                    + "n      1  2  2  2  3  3  3  3  3  4\n" \
                    + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1\n" \
                    + "Occup  2  2  2  4  2  2  4  4  4  2\n" \
                    + "Valen  0  0  0  0  0  0  0  1  1  1\n"
        elif atype == "Cu":
            results = "Cu\n" \
                + "Iz=  29 Norb= 10 Ion=  0 Config= 3d10_4s1\n" \
                + "n      1  2  2  2  3  3  3  3  3  4\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  1\n" \
                + "Valen  0  0  0  0  0  0  0  1  1  1\n"
        elif atype == "Zn":
            results = "Zn\n" \
                + "Iz=  30 Norb= 10 Ion=  0 Config= 3d10_4s2\n" \
                + "n      1  2  2  2  3  3  3  3  3  4\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2\n" \
                + "Valen  0  0  0  0  0  0  0  1  1  1\n"
        elif atype == "Ga":
            results = "Ga\n" \
                + "Iz=  31 Norb= 11 Ion=  0 Config= 3d10_4s2_4p1\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  1\n" \
                + "Valen  0  0  0  0  0  0  0  1  1  1  1\n"
        elif atype == "Ge":
            results = "Ge\n" \
                + "Iz=  32 Norb= 11 Ion=  0 Config= 4s2_4p2\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2\n" \
                + "Valen  0  0  0  0  0  0  0  0  0  1  1\n"
        elif atype == "As":
            results = "As\n" \
                + "Iz=  33 Norb= 12 Ion=  0 Config= 3d10_4s2_4p3\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  1\n" \
                + "Valen  0  0  0  0  0  0  0  1  1  1  1  1\n"
        elif atype == "Se":
            results = "Se\n" \
                + "Iz=  34 Norb= 12 Ion=  0 Config= 3d10_4s2_4p4\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  2\n" \
                + "Valen  0  0  0  0  0  0  0  1  1  1  1  1\n"
        elif atype == "Br":
            results = "Br\n" \
                + "Iz=  35 Norb= 12 Ion=  0 Config= 3d10_4s2_4p5\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  3\n" \
                + "Valen  0  0  0  0  0  0  0  1  1  1  1  1\n"
        elif atype == "Kr":
            results = "Kr\n" \
                + "Iz=  36 Norb= 12 Ion=  0 Config= 3d10_4s2_4p6\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4\n" \
                + "Valen  0  0  0  0  0  0  0  1  1  1  1  1\n"
        elif atype == "Rb":
            results = "Rb\n" \
                + "Iz=  37 Norb= 13 Ion=  0 Config= 5s1\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  5\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  1\n" \
                + "Valen  0  0  0  0  0  0  0  0  0  0  0  0  1\n"

        elif atype == "Sr":
            results = "Sr\n" \
                + "Iz=  38 Norb= 13 Ion=  0 Config= 5s2\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  5\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  2\n" \
                + "Valen  0  0  0  0  0  0  0  0  0  0  0  0  1\n"
        #elif atype == "Sr":
        #    results = "Sr\n" \
        #            + "Iz=  38 Norb= 13 Ion=  0 Config= 4p6_5s2\n" \
        #            + "n      1  2  2  2  3  3  3  3  3  4  4  4  5\n" \
        #            + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2 -1\n" \
        #            + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  2\n" \
        #            + "Valen  0  0  0  0  0  0  0  0  0  0  1  1  1\n"
        elif atype == "Y":
            results = "Y\n" \
                + "Iz=  39 Norb= 14 Ion=  0 Config= 4d1_5s2\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  5\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  1  2\n" \
                + "Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1\n"
        elif atype == "Zr":
            results = "Zr\n" \
                + "Iz=  40 Norb= 14 Ion=  0 Config= 4d2_5s2\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  5\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  2  2\n" \
                + "Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1\n"
        elif atype == "Nb":
            results = "Nb\n" \
                + "Iz=  41 Norb= 14 Ion=  0 Config= 4d3_5s2\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  5\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  3  2\n" \
                + "Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1\n"
        elif atype == "Mo":
            results = "Mo\n" \
                + "Iz=  42 Norb= 14 Ion=  0 Config= 4d4_5s2\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  5\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  2\n" \
                + "Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1\n"
        elif atype == "Tc":
            results = "Tc\n" \
                + "Iz=  43 Norb= 15 Ion=  0 Config= 4d5_5s2\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  5\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  1  2\n" \
                + "Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1\n"
        elif atype == "Ru":
            results = "Ru\n" \
                + "Iz=  44 Norb= 15 Ion=  0 Config= 4d6_5s2\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  5\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  2  2\n" \
                + "Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1\n"
        elif atype == "Rh":
            results = "Rh\n" \
                + "Iz=  45 Norb= 15 Ion=  0 Config= 4d7_5s2\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  5\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  3  2\n" \
                + "Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1\n"
        elif atype == "Pd":
            results = "Pd\n" \
                + "Iz=  46 Norb= 15 Ion=  0 Config= 4d8_5s2\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  5\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  4  2\n" \
                + "Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1\n"
        elif atype == "Ag":
            results = "Ag\n" \
                + "Iz=  47 Norb= 15 Ion=  0 Config= 4d10_5s1\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  5\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  1\n" \
                + "Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1\n"
        elif atype == "Cd":
            results = "Cd\n" \
                + "Iz=  48 Norb= 15 Ion=  0 Config= 4d10_5s2\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  5\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  2\n" \
                + "Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1\n"
        elif atype == "In":
            results = "In\n" \
                + "Iz=  49 Norb= 16 Ion=  0 Config= 4d10_5s2_5p1\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  5  5\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3 -1  1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  2  1\n" \
                + "Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1\n"
        elif atype == "Sn":
            results = "Sn\n" \
                + "Iz=  50 Norb= 16 Ion=  0 Config= 4d10_5s2_5p2\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  5  5\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3 -1  1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  2  2\n" \
                + "Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1\n"
        elif atype == "Sb":
            results = "Sb\n" \
                + "Iz=  51 Norb= 17 Ion=  0 Config= 4d10_5s2_5p3\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  5  5  5\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3 -1  1 -2\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  2  2  1\n" \
                + "Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1\n"
        elif atype == "Te":
            results = "Te\n" \
                + "Iz=  52 Norb= 17 Ion=  0 Config= 4d10_5s2_5p4\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  5  5  5\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3 -1  1 -2\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  2  2  2\n" \
                + "Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1\n"
        elif atype == "I":
            results = "I\n" \
                + "Iz=  53 Norb= 17 Ion=  0 Config= 4d10_5s2_5p5\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  5  5  5\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3 -1  1 -2\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  2  2  3\n" \
                + "Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1\n"
        elif atype == "Xe":
            results = "Xe\n" \
                + "Iz=  54 Norb= 17 Ion=  0 Config= 4d10_5s2_5p6\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  5  5  5\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3 -1  1 -2\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  2  2  4\n" \
                + "Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1\n"
        elif atype == "Cs":
            results = "Cs\n" \
                + "Iz=  55 Norb= 18 Ion=  0 Config= 6s1\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  5  5  5  6\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3 -1  1 -2 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  2  2  4  1\n" \
                + \
                "Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1\n"
        elif atype == "Ba":
            results = "Ba\n" \
                + "Iz=  56 Norb= 18 Ion=  0 Config= 6s2\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  5  5  5  6\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3 -1  1 -2 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  2  2  4  2\n" \
                + \
                "Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1\n"
        elif atype == "La":
            results = "La\n" \
                + "Iz=  57 Norb= 19 Ion=  0 Config= 5d2_6s1\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  5  5  5  5  6\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3 -1  1 -2  2 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  2  2  4  1  2\n" \
                + \
                "Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1\n"
        elif atype == "Ce":
            results = "Ce\n" \
                + "Iz=  58 Norb= 20 Ion=  0 Config= 4f1_5d2_6s1\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  5  5  5  5  6\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -1  1 -2  2 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  1  2  2  4  2  1\n" \
                + \
                "Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  1  1\n"
        elif atype == "Dy":
            results = "Dy\n" \
                + "Iz=  66 Norb= 21 Ion=  0 Config= 4f9 5d1_6s2\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  6\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  3  2  2  4  1  2\n" \
                + \
                "Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1\n"
        elif atype == "Tm":
            results = "Tm\n" \
                + "Iz=  69 Norb= 21 Ion=  0 Config= 4f12 5d1_6s2\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  6\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  6  2  2  4  1  2\n" \
                + \
                "Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1\n"
        elif atype == "Lu":
            results = "Lu\n" \
                + "Iz=  71 Norb= 21 Ion=  0 Config= 5d1_6s2\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  6\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  1  2\n" \
                + \
                "Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1\n"
        elif atype == "Hf":
            results = "Hf\n" \
                + "Iz=  72 Norb= 21 Ion=  0 Config= 5d2_6s2\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  6\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  2  2\n" \
                + \
                "Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1\n"
        elif atype == "Ta":
            results = "Ta\n" \
                + "Iz=  73 Norb= 21 Ion=  0 Config= 5d3_6s2\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  6\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  3  2\n" \
                + \
                "Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1\n"
        elif atype == "W":
            results = "W\n" \
                + "Iz=  74 Norb= 21 Ion=  0 Config= 5d4_6s2\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  6\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  2\n" \
                + \
                "Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1\n"
        elif atype == "Re":
            results = "Re\n" \
                + "Iz=  75 Norb= 22 Ion=  0 Config= 5d5_6s2\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  6\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  1  2\n" \
                + \
                "Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1\n"
        elif atype == "Os":
            results = "Os\n" \
                + "Iz=  76 Norb= 22 Ion=  0 Config= 5d6_6s2\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  6\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  2  2\n" \
                + \
                "Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1\n"
        elif atype == "Ir":
            results = "Ir\n" \
                + "Iz=  77 Norb= 22 Ion=  0 Config= 5d7_6s2\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  6\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  3  2\n" \
                + \
                "Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1\n"
        elif atype == "Pt":
            results = "Pt\n" \
                + "Iz=  78 Norb= 22 Ion=  0 Config= 5d8_6s2\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  6\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  4  2\n" \
                + \
                "Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1\n"
        elif atype == "Au":
            results = "Au\n" \
                + "Iz=  79 Norb= 22 Ion=  0 Config= 5d10_6s1\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  6\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  6  1\n" \
                + \
                "Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1\n"
        elif atype == "Hg":
            results = "Hg\n" \
                + "Iz=  80 Norb= 22 Ion=  0 Config= 5d10_6s2\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  6\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  6  2\n" \
                + \
                "Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1\n"
        elif atype == "Tl":
            results = "Tl\n" \
                + "Iz=  81 Norb= 23 Ion=  0 Config= 5d10_6s2_6p1\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  6  6\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3 -1  1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  6  2  1\n" \
                + \
                "Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1\n"
        elif atype == "Pb":
            results = "Pb\n" \
                + "Iz=  82 Norb= 23 Ion=  0 Config= 6s2_6p2\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  6  6\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3 -1  1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  6  2  2\n" \
                + \
                "Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1\n"
        elif atype == "Bi":
            results = "Bi\n" \
                + "Iz=  83 Norb= 24 Ion=  0 Config= 5d10_6s2_6p3\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  6  6  6\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3 -1  1 -2\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  6  2  2  1\n" \
                + \
                "Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1\n"
        elif atype == "Po":
            results = "Po\n" \
                + "Iz=  84 Norb= 24 Ion=  0 Config= 5d10_6s2_6p4\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  6  6  6\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3 -1  1 -2\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  6  2  2  2\n" \
                + \
                "Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1\n"
        elif atype == "At":
            results = "At\n" \
                + "Iz=  85 Norb= 24 Ion=  0 Config= 5d10_6s2_6p5\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  6  6  6\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3 -1  1 -2\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  6  2  2  3\n" \
                + \
                "Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1\n"
        elif atype == "Rn":
            results = "Rn\n" \
                + "Iz=  86 Norb= 24 Ion=  0 Config= 5d10_6s2_6p6\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  6  6  6\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3 -1  1 -2\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  6  2  2  4\n" \
                + \
                "Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1\n"
        elif atype == "Fr":
            results = "Fr\n" \
                + "Iz=  87 Norb= 25 Ion=  0 Config= 7s1\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  6  6  6  7\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3 -1  1 -2 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  6  2  2  4  1\n" \
                + \
                "Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1\n"
        elif atype == "Ra":
            results = "Ra\n" \
                + "Iz=  88 Norb= 25 Ion=  0 Config= 7s2\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  6  6  6  7\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3 -1  1 -2 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  6  2  2  4  2\n" \
                + \
                "Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1\n"
        elif atype == "Ac":
            results = "Ac\n" \
                + "Iz=  89 Norb= 26 Ion=  0 Config= 6d2_7s1\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  6  6  6  6  7\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3 -1  1 -2  2 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  6  2  2  4  2  1\n" \
                + \
                "Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1\n"
        elif atype == "Th":
            results = "Th\n" \
                + "Iz=  90 Norb= 26 Ion=  0 Config= 6d3_7s1\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  6  6  6  6  7\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3 -1  1 -2  2 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  6  2  2  4  3  1\n" \
                + \
                "Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1\n"
        elif atype == "Pa":
            results = "Pa\n" \
                + "Iz=  91 Norb= 27 Ion=  0 Config= 5f1_6d3_7s1\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  5  6  6  6  6  7\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3  3 -1  1 -2  2 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  6  1  2  2  4  3  1\n" \
                + \
                "Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  1  1\n"
        elif atype == "U":
            results = "U\n" \
                + "Iz=  92 Norb= 27 Ion=  0 Config= 5f2_6d3_7s1\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  5  6  6  6  6  7\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3  3 -1  1 -2  2 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  6  2  2  2  4  3  1\n" \
                + \
                "Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  1  1\n"
        elif atype == "Np":
            results = "Np\n" \
                + "Iz=  93 Norb= 27 Ion=  0 Config= 5f4_6d2_7s1\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  5  6  6  6  6  7\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3  3 -1  1 -2  2 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  6  4  2  2  4  2  1\n" \
                + \
                "Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  1  1\n"
        elif atype == "Pu":
            results = "Pu\n" \
                + "Iz=  94 Norb= 27 Ion=  0 Config= 5f5_6d2_7s1\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  5  6  6  6  6  7\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3  3 -1  1 -2  2 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  6  5  2  2  4  2  1\n" \
                + \
                "Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  1  1\n"
        elif atype == "Am":
            results = "Am\n" \
                + "Iz=  95 Norb= 27 Ion=  0 Config= 5f6_6d2_7s1\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  5  6  6  6  6  7\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3  3 -1  1 -2  2 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  6  6  2  2  4  2  1\n" \
                + \
                "Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  1  1\n"
        elif atype == "Cm":
            results = "Cm\n" \
                + "Iz=  96 Norb= 28 Ion=  0 Config= 5f7_6d2_7s1\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  5  5  6  6  6  6  7\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  6  6  1  2  2  4  2  1\n" \
                + \
                "Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  0  0  0  1  1\n"
        elif atype == "Bk":
            results = "Bk\n" \
                + "Iz=  97 Norb= 28 Ion=  0 Config= 5f8_6d2_7s1\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  5  5  6  6  6  6  7\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  6  6  2  2  2  4  2  1\n" \
                + \
                "Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  0  0  0  1  1\n"
        elif atype == "Cf":
            results = "Cf\n" \
                + "Iz=  98 Norb= 28 Ion=  0 Config= 5f9_6d2_7s1\n" \
                + "n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  5  5  6  6  6  6  7\n" \
                + "Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -1\n" \
                + "Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  6  6  3  2  2  4  2  1\n" \
                + \
                "Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  0  0  0  1  1\n"
        else:
            sys.exit(
                "ERROR: Atom configuration for {0} not implemented!".format(atype))

        return results

    def output(self):
        """(self) -> (str)

            Output first part of the kgrn input file in formated string

        :returns:
        :rtype:
        """

        now = datetime.datetime.now()
        line = "KGRN                                              "\
            + str(now.day) + "." + str(now.month) + "." + str(now.year) + "\n"
        line = line + "JOBNAM=" + self.jobname + "\n"
        line = line + "STRT..=  " + self.strt + " MSGL.=  " + str(self.msgl) \
            + " EXPAN.= " + self.expan + " FCD..=  " + self.fcd \
            + " FUNC..= " + self.func + "\n"
        line = line + "FOR001=" + self.FOR001 + "\n"
        line = line + "FOR001=" + self.FOR001_2 + "\n"
        line = line + "DIR002=" + self.DIR002 + "\n"
        line = line + "DIR003=" + self.DIR003 + "\n"
        line = line + "FOR004=" + self.FOR004 + "\n"
        line = line + "DIR006=" + self.DIR006 + "\n"
        line = line + "DIR009=" + self.DIR009 + "\n"
        line = line + "DIR010=" + self.DIR010 + "\n"
        line = line + "DIR011=" + self.DIR011 + "\n"
        line = line + "Self-consistent KKR Calc for %s.\n" % (self.jobname)
        line = line + "Band: 10 lines\n"
        line = line + "NITER.=%3i NLIN.=%3i NPRN.=  %1i NCPA.=%3i " \
            % (self.niter, self.nlin, self.nprn, self.ncpa) \
            + "NT...= %2i MNTA.= %2i\n" % (self.nt, self.mnta)
        line = line + "MODE..= %2s FRC..=  %1s DOS..=  %1s OPS..=  %1s " \
            % (self.mode, self.frc, self.dos, self.ops) + \
            "AFM..=  " + self.afm + " CRT..=  " + self.crt + "\n"
        line = line + "Lmaxh.= %2i Lmaxt= %2i NFI..=%3i " \
            % (self.lmaxh, self.lmaxt, self.kgrn_nfi) \
            + "FIXG.= %2i SHF..=  %1i SOFC.=  %1s" \
            % (self.fixg, self.shf, self.sofc) + "\n"
        line = line + "KMSH...= %1s IBZ..= %2i NKX..= %2i " \
            % (self.kmsh, self.ibz, self.nkx) \
            + "NKY..= %2i NKZ..= %2i FBZ..=  %1s\n" \
            % (self.nky, self.nkz, self.fbz)
        line = line + "KMSH2..= %1s IBZ2.=%3i NKX2.=%3i NKY2.=%3i NKZ2.=%3i" \
            % (self.kmsh2, self.ibz2, self.nkx2, self.nky2, self.nkz2) + "\n"
        line = line + "ZMSH...= %1s NZ1..= %2i " % (self.zmsh, self.nz1) \
            + "NZ2..=%3i NZ3..=%3i NRES.=%3i NZD.=%4i" \
            % (self.nz2, self.nz3, self.nres, self.nzd) + "\n"
        line = line + "DEPTH..= %6.3f IMAGZ.= %6.3f " \
            % (self.depth, self.imagz) \
            + "EPS...=  %5.3f ELIM..= %6.3f" \
            % (self.eps, self.elim) + "\n"
        line = line + "AMIX...= %6.3f EFMIX.= %6.3f " % (self.amix, self.efmix) \
            + "VMTZ..=%7.3f MMOM..=%7.3f" \
            % (self.vmtz, self.mmom) + "\n"
        line = line + "TOLE...=%7.1e TOLEF.=%7.1e " % (self.tole, self.tolef) \
            + "TOLCPA=%7.1e TFERMI= %6.1f\n" % (self.tolcpa, self.tfermi)
        line = line + "SWS......=%10.7f NSWS.=%3i " % (self.sws, self.nsws) \
            + "DSWS..=   %4.2f ALPCPA= %6.4f\n" % (self.dsws, self.alpcpa)
        line = line + "Setup: 2 + NQ*NS lines\n"
        line = line + "EFGS...= %6.3f HX....= %6.3f " % (self.efgs, self.hx) \
            + "NX...=  %1i NZ0..= %2i " % (self.nx, self.nz0) \
            + "STMP..= " + self.stmp + "\n"
        line = line + "Symb   IQ IT ITA NZ  CONC   "\
            + "Sm(s)  S(ws) WS(wst) QTR SPLT Fix\n"
        line = line + self.AtomOutput()
        line = line + "Atom:  4 lines + NT*6 lines\n"
        line = line + "IEX...= %2i NP..=%4i " % (self.iex, self.dirac_np) \
            + "NES..=%3i NITER=%3i IWAT.=%3i NPRNA=%3i" \
            % (self.nes, self.dirac_niter, self.iwat, self.nprna) + "\n"
        line = line + "VMIX.....=  %8.6f RWAT....=  %8.6f RMAX....=%10.6f" \
            % (self.vmix, self.rwat, self.rmax) + "\n"
        line = line + "DX.......=  %8.6f DR1.....=  %8.6f TEST....=  %8.2e" \
            % (self.dx, self.dr1, self.test) + "\n"
        line = line + "TESTE....=  %8.2e TESTY...=  %8.2e TESTV...=  %8.2e" \
            % (self.teste, self.testy, self.testv) + "\n"

        for i in self.confidx:
            line = line + self.aconflines(self.conforder[i][0])

        return(line)

    def write_input_file(self, folder=None):
        """(self,str) ->(None)

            Save KGRN input data to file named filename

        :param folder:  (Default value = None)
        :type folder:
        :returns:
        :rtype:
        """

        # Check data integrity before anything is written on disk or run
        self.check_input_file()

        if folder is None:
            #sys.exit('Kgrn.create_input_file: \'folder\' has to be given!')
            folder = './'
        else:
            common.check_folders(folder)

        fl = open(folder + '/{0}.kgrn'.format(self.jobname), "w")
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
        if not hasattr(self, key):
            print('WARNING: Kgrn() class has no attribute \'{0}\''.format(key))
            return

        if key in self.atconfKeys:
            setattr(self, key, value)
            # Make sure the element symbols are in the proper format
            if key == 'atoms':
                for i in range(len(self.atoms)):
                    self.atoms[i] = self.atoms[i].lower().capitalize()
            # Try creating the atomblock
            self.create_atomblock()

        # Lattice name or path has changed => we have to update the FOR and DIR
        # information
        elif key == 'latname' or key == 'latpath':
            setattr(self, key, value)
            if self.latname is not None and self.latpath is not None:
                self.FOR001 = self.latpath + '/kstr/' + self.latname + '.tfh'
                self.FOR001_2 = self.latpath + \
                    '/kstr/' + self.latname + '2.tfh'
                self.FOR001 = common.cleanup_path(self.FOR001)
                self.FOR001_2 = common.cleanup_path(self.FOR001_2)
                self.FOR004 = self.latpath + '/bmdl/' + self.latname + '.mdl'
                self.FOR004 = common.cleanup_path(self.FOR004)

        else:
            setattr(self, key, value)

        return

    def create_atconf(self):
        """Constructs the self.atconf list out of the atomic information.

        :returns: None
        :rtype: None
        """

        self.atconf = []
        for i in range(len(self.atoms)):
            self.atconf.append([self.atoms[i], self.iqs[i], self.its[i], self.itas[i],
                                self.concs[i], self.sm_ss[
                                    i], self.s_wss[i], self.ws_wsts[i],
                                self.qtrs[i], self.splts[i], self.fixs[i]])
        return

    def create_atomblock(self):
        """Constructs the KGRN input file atomblock if all
        the necessary parameters are present.

        :returns: None
        :rtype: None
        """

        for key in self.atconfKeys:
            if getattr(self, key) == None:
                # Some information needed to create the atomblock is still missing,
                # postpone further processing.
                return

        # Every attribute necessary has some meaningful value, atomblock can therefore
        # be constructed.

        self.atomblock = []

        # First we must create the atconf-array because the code below constructs
        # the atomblock using the atconf-array.

        self.create_atconf()
        ln = len(self.atconf)
        self.confidx = []

        self.aconfig = np.empty(ln, dtype=self.Confdtype)

        for i in range(ln):
            self.aconfig[i] = (self.atconf[i][0], self.atconf[i][1], self.atconf[i][2],
                               self.atconf[i][3], self.atconf[
                                   i][4], self.atconf[i][5],
                               self.atconf[i][6], self.atconf[
                                   i][7], self.atconf[i][8],
                               self.atconf[i][9], self.atconf[i][10])

        # Sort list so that figuring out the first appearances
        # of unique [it,ita]-pairs is easy.
        self.conforder = np.sort(self.aconfig, order=['it', 'ita', 'iq'])

        self.its_itas = []
        for i in range(ln):
            it = self.conforder[i][2]
            ita = self.conforder[i][3]
            it_ita = [it, ita]
            if it_ita not in self.its_itas:
                self.its_itas.append(it_ita)
                self.confidx.append(i)

        # Calculate NT and MNTA so that they don't have to be changed manually
        self.its_itas = np.array(self.its_itas)
        self.nt = np.max(self.its_itas[:, 0])
        self.mnta = np.max(self.its_itas[:, 1])

        # Sort list in terms of iq, which will be the order
        # in which the lines appear in the final input file.
        self.aconfig = np.sort(self.aconfig, order=['iq', 'it', 'ita'])
        for i in range(ln):
            self.atomblock.append(
                self.Atomline(self.aconfig[i][0], self.aconfig[i][1], self.aconfig[i][2],
                              self.aconfig[i][3], self.aconfig[
                                  i][4], self.aconfig[i][5],
                              self.aconfig[i][6], self.aconfig[
                                  i][7], self.aconfig[i][8],
                              self.aconfig[i][9], self.aconfig[i][10]))
        return

    def check_input_file(self):
        """Perform various checks on the class data

        Makes sure that all necessary data exists
        before we attempt to write the input file to disk.

        :returns: None
        :rtype: None
        """

        # Mission critical arguments
        if self.jobname is None:
            sys.exit('Kgrn: \'jobname\' has to be given!')
        if self.latname is None:
            sys.exit('Kgrn: \'latname\' has to be given!')
        if self.ibz is None:
            try:
                self.ibz = common.lat_to_ibz(self.latname)
            except KeyError:
                sys.exit(
                    'Kgrn: \'ibz\' has to be given for \'latname\'={0}!'.format(self.latname))
        if self.latpath is None:
            self.latpath = './'
        if self.sws is None or self.sws == 0.0:
            sys.exit('Kgrn: \'sws\' has to be given!')
        if self.atomblock is None:
            sys.exit('Kgrn: Some atomic information is missing =>' +
                     ' self.atomblock has not been constructed!')
        if self.atconf is None:
            sys.exit('Kgrn: Some atomic information is missing =>' +
                     ' self.atconf has not been constructed!')
        else:
            if isinstance(self.atconf[0], list):
                pass
            else:
                self.atconf = [self.atconf]

        # Arguments for which reasonable default values exist
        if self.strt is None:
            self.strt = 'A'
        if self.msgl is None:
            self.msgl = 0
        if self.expan is None:
            self.expan = 'S'
        if self.fcd is None:
            self.fcd = 'Y'
        if self.func is None:
            self.func = 'SCA'
        if self.niter is None:
            self.niter = 100
        if self.nlin is None:
            self.nlin = 31
        if self.nprn is None:
            self.nprn = 0
        if self.ncpa is None:
            self.ncpa = 7
        if self.mode is None:
            self.mode = '3D'
        if self.frc is None:
            self.frc = 'N'
        if self.dos is None:
            self.dos = 'N'
        if self.ops is None:
            self.ops = 'N'
        if self.afm is None:
            self.afm = 'P'
        if self.crt is None:
            self.crt = 'M'
        if self.lmaxh is None:
            self.lmaxh = 8
        if self.lmaxt is None:
            self.lmaxt = 4
        if self.kgrn_nfi is None:
            self.kgrn_nfi = 31
        if self.fixg is None:
            self.fixg = 2
        if self.shf is None:
            self.shf = 0
        if self.sofc is None:
            self.sofc = 'N'
        if self.kmsh is None:
            self.kmsh = 'G'
        if self.nkx is None:
            self.nkx = 0
        if self.nky is None:
            if self.ibz == 2:
                self.nky = 13
            elif self.ibz == 3:
                self.nky = 17
            elif self.ibz == 4:
                self.nky = 13
                self.nkz = 9
            else:
                self.nky = 13
        if self.nkz is None:
            self.nkz = 0
        if self.fbz is None:
            self.fbz = 'N'
        if self.kmsh2 is None:
            self.kmsh2 = 'G'
        if self.ibz2 is None:
            self.ibz2 = 1
        if self.nkx2 is None:
            self.nkx2 = 4
        if self.nky2 is None:
            self.nky2 = 0
        if self.nkz2 is None:
            self.nkz2 = 51
        if self.zmsh is None:
            self.zmsh = 'C'
        if self.nz1 is None:
            self.nz1 = 16
        if self.nz2 is None:
            self.nz2 = 8
        if self.nz3 is None:
            self.nz3 = 8
        if self.nres is None:
            self.nres = 4
        if self.nzd is None:
            self.nzd = 1500
        if self.depth is None:
            self.depth = 1.0
        if self.imagz is None:
            self.imagz = 0.020
        if self.eps is None:
            self.eps = 0.2
        if self.elim is None:
            self.elim = -1.0
        if self.amix is None:
            self.amix = 0.05
        if self.efmix is None:
            self.efmix = 1.0
        if self.vmtz is None:
            self.vmtz = 0.0
        if self.mmom is None:
            self.mmom = 0.0
        if self.tole is None:
            self.tole = 1.0e-7
        if self.tolef is None:
            self.tolef = 1.0e-7
        if self.tolcpa is None:
            self.tolcpa = 1.0e-6
        if self.tfermi is None:
            self.tfermi = 500.0
        if self.nsws is None:
            self.nsws = 1
        if self.dsws is None:
            self.dsws = 0.05
        if self.alpcpa is None:
            self.alpcpa = 0.6020
        if self.efgs is None:
            self.efgs = 0.0
        if self.hx is None:
            self.hx = 0.1
        if self.nx is None:
            self.nx = 5
        if self.nz0 is None:
            self.nz0 = 6
        if self.stmp is None:
            self.stmp = 'Y'
        if self.iex is None:
            self.iex = 4
        if self.dirac_np is None:
            self.dirac_np = 251
        if self.nes is None:
            self.nes = 15
        if self.dirac_niter is None:
            self.dirac_niter = 100
        if self.iwat is None:
            self.iwat = 0
        if self.nprna is None:
            self.nprna = 0
        if self.vmix is None:
            self.vmix = 0.3
        if self.rwat is None:
            self.rwat = 3.5
        if self.rmax is None:
            self.rmax = 20.0
        if self.dx is None:
            self.dx = 0.03
        if self.dr1 is None:
            self.dr1 = 0.002
        if self.test is None:
            self.test = 1.0e-12
        if self.teste is None:
            self.teste = 1.0e-12
        if self.testy is None:
            self.testy = 1.0e-12
        if self.testv is None:
            self.testv = 1.0e-12
        if self.FOR001 is None:
            self.FOR001 = self.latpath + '/kstr/' + self.latname + '.tfh'
            self.FOR001_2 = self.latpath + '/kstr/' + self.latname + '2.tfh'
            self.FOR001 = common.cleanup_path(self.FOR001)
            self.FOR001_2 = common.cleanup_path(self.FOR001_2)
        if self.DIR002 is None:
            self.DIR002 = 'kgrn/'
        if self.DIR003 is None:
            self.DIR003 = 'kgrn/'
        if self.FOR004 is None:
            self.FOR004 = self.latpath + '/bmdl/' + self.latname + '.mdl'
            self.FOR004 = common.cleanup_path(self.FOR004)
        if self.DIR006 is None:
            self.DIR006 = 'kgrn/'
        if self.DIR009 is None:
            self.DIR009 = 'kgrn/'
        if self.DIR010 is None:
            self.DIR010 = 'kgrn/'
        if self.DIR011 is None:
            self.DIR011 = 'kgrn/tmp/'
        return

