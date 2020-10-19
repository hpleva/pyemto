# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 14:48:25 2014

@author: Matti Ropo
@author: Henrik Levämäki

"""

import sys
import os
import datetime
import pickle
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
    :param fxms:  (Default value = None)
    :type fxms:
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
    :param FOR002:  (Default value = None)
    :type FOR002:
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
                 concs=None, iqs=None, its=None, itas=None, qtrs=None, splts=None, fxms=None,
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
                 FOR001=None, DIR002=None, DIR003=None, FOR002=None, DIR006=None, DIR009=None,
                 DIR010=None, DIR011=None,ncpu=None,CQNA=None,KGRN_file_type=None,
                 setups=None, gpm='N', fsm='N', DIR021='', FOR098='./ATOM.cfg',
                 DIR022=None, vmixatm=None, kpole=None,
                 nrms=None, a_scrs=None, b_scrs=None, tetas=None, phis=None,
                 qx=None, qy=None, qz=None):

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
        self.fxms = fxms
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
        self.FOR002 = FOR002
        self.DIR002 = DIR002
        self.DIR003 = DIR003
        self.DIR006 = DIR006
        self.DIR009 = DIR009
        self.DIR010 = DIR010
        self.DIR011 = DIR011
        self.ncpu = ncpu
        self.CQNA = CQNA
        self.KGRN_file_type = KGRN_file_type
        self.setups = setups
        self.gpm = gpm
        self.fsm = fsm
        self.DIR021 = DIR021
        self.DIR022 = DIR022
        self.FOR098 = FOR098
        self.vmixatm = vmixatm
        self.kpole = kpole
        self.qx = qx
        self.qy = qy
        self.qz = qz

        if a_scrs is None:
            self.a_scrs = np.ones(50) * 0.67
        else:
            self.a_scrs = a_scrs
        if b_scrs is None:
            self.b_scrs = np.ones(50) * 1.05
        else:
            self.b_scrs = b_scrs
        if tetas is None:
            self.tetas = np.zeros(50)
        else:
            self.tetas = tetas
        if phis is None:
            self.phis = np.zeros(50)
        else:
            self.phis = phis
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
        if nrms is None:
            self.nrms = np.ones(50)
        else:
            self.nrms = nrms
        if qtrs is None:
            self.qtrs = np.zeros(50)
        else:
            self.qtrs = qtrs
        if splts is None:
            self.splts = np.zeros(50)
        else:
            self.splts = splts
        if fxms is None:
            self.fxms = np.empty([50], dtype='S1')
            self.fxms[:] = 'N'
        else:
            self.fxms = fxms
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
        if KGRN_file_type is None:
            self.KGRN_file_type = 'kgrn'


        self.atconfKeys = ['atoms', 'concs', 'iqs', 'its', 'itas', 'qtrs', 'splts',
                           'fxms', 'sm_ss', 's_wss', 'ws_wsts']

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

        # self.Confdtype = [('elem', 'U4'), ('iq', int), ('it', int), ('ita', int),
                          # ('conc', float), ('Sm_s', float), ('S_ws',float), ('WS_wst', float),
                          # ('qtr', float), ('mmom', float), ('fix', 'U1')]

        self.Confdtype = [('elem', 'U4'), ('iq', int), ('it', int),
                ('ita', int), ('nrm', int), ('conc', float),
                ('a_scr', float), ('b_scr', float), ('teta', float),
                ('phi', float), ('fxm', 'U1'), ('m_split', float)]

        self.default_setups = {
            "Va": "1s", "Em": "1s", "H": "1s", "He": "1s", "Li": "2s",
            "Be": "2s", "B": "2s2p", "C": "2s2p", "N": "2p", "O": "2s2p",
            "O-2": "2p", "F": "2s2p", "Ne": "2s2p", "Na": "3s", "Mg": "3s",
            "Al": "3s3p", "Si": "3s3p", "P": "3s3p", "S": "3s3p", "Cl": "3s3p",
            "Ar": "3s3p", "K": "4s", "Ca": "4s", "Sc": "3d4s", "Ti": "3d4s",
            "V": "3d4s", "Cr": "3d4s", "Mn": "3d4s", "Fe": "3d4s", "Co": "3d4s",
            "Ni": "3d4s", "Cu": "3d4s", "Zn": "3d4s", "Ga": "3d4s4p",
            "Ge": "4s4p", "As": "3d4s4p", "Se": "3d4s4p", "Br": "3d4s4p",
            "Kr": "3d4s4p", "Rb": "5s", "Sr": "5s", "Y": "4d5s", "Zr": "4d5s",
            "Nb": "4d5s", "Mo": "4d5s", "Tc": "4d5s", "Ru": "4d5s", "Rh": "4d5s",
            "Pd": "4d5s", "Ag": "4d5s", "Cd": "4d5s", "In": "4d5s5p", "Sn": "4d5s5p",
            "Sb": "4d5s5p", "Te": "4d5s5p", "I": "4d5s5p", "Xe": "4d5s5p",
            "Cs": "6s", "Ba": "6s", "La": "5d6s", "Ce": "4f5d6s",
            "Dy": "5d6s", "Tm": "5d6s", "Lu": "5d6s", "Hf": "5d6s", "Ta": "5d6s",
            "W": "5d6s", "Re": "5d6s", "Os": "5d6s", "Ir": "5d6s", "Pt": "5d6s",
            "Au": "5d6s", "Hg": "5d6s", "Tl": "5d6s6p", "Pb": "6s6p",
            "Bi": "5d6s6p", "Po": "5d6s6p", "At": "5d6s6p", "Rn": "5d6s6p",
            "Fr": "7s", "Ra": "7s", "Ac": "6d7s", "Th": "6d7s", "Pa": "5f6d7s",
            "U": "5f6d7s", "Np": "5f6d7s", "Pu": "5f6d7s", "Am": "5f6d7s",
            "Cm": "5f6d7s", "Bk": "5f6d7s", "Cf": "5f6d7s",
        }

        with open(os.path.join(os.path.dirname(__file__), "aconf.pckl"), "rb") as f:
            self.aconf_dict = pickle.load(f)

        # Try creating the atomblock
        self.create_atomblock()

    def Atomline(self, inp):
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
        atom, iq, it, ita, nrm, conc, a_scr, b_scr, teta, phi, fxm, m_split = inp

        if np.round(conc, decimals=3) >= 10:
            line = "%-2s    %2i  %2i  %2i   %1i  %8.5f  %5.3f %5.3f  %6.4f  %6.4f  %s   %6.4f" \
                   % (atom, iq, it, ita, nrm, conc, a_scr, b_scr, teta, phi,
                      fxm, m_split)
        else:
            line = "%-2s    %2i  %2i  %2i   %1i  %8.6f  %5.3f %5.3f  %6.4f  %6.4f  %s   %6.4f" \
                   % (atom, iq, it, ita, nrm, conc, a_scr, b_scr, teta, phi,
                      fxm, m_split)

        return line

    def AtomOutput(self):
        """(self) -> str

            Output of atomic lines in kgrn format

        :returns: A atomic lines for kgrn input
        :rtype: str
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

        if self.setups is not None:
            if atype in self.setups.keys():
                key = atype + "_" + self.setups[atype]
            else:
                key = atype + "_" + self.default_setups[atype]
        else:
            key = atype + "_" + self.default_setups[atype]

        results = self.aconf_dict[key]
        return results

    def output(self):
        """(self) -> (str)

            Output first part of the kgrn input file in formated string

        :returns: A first part of the kgrn input file
        :rtype: str
        """

        now = datetime.datetime.now()
        line = "KGRN      HP..= 0                                 "\
            + str(now.day) + "." + str(now.month) + "." + str(now.year) + "\n"
        line += "JOBNAM...=" + self.jobname + "\n"
        line += "MSGL.=" + str(self.msgl) + " STRT.=" + self.strt  \
            + " FUNC.=" + self.func + " EXPAN=" + str(self.expan) + \
            " FCD.=" + self.fcd + " GPM.=" + self.gpm + \
            " FSM.=" + self.fsm + "\n"
        line += "FOR001=" + self.FOR001 + "\n"
        line += "FOR002=" + self.FOR002 + "\n"
        line += "DIR003=" + self.DIR003 + "\n"
        line += "DIR006=" + self.DIR006 + "\n"
        line += "DIR010=" + self.DIR010 + "\n"
        line += "DIR011=" + self.DIR011 + "\n"
        line += "DIR021=" + self.DIR021 + "\n"
        line += "DIR022=" + self.DIR022 + "\n"
        line += "FOR098=" + self.FOR098 + "\n"
        line += "Self-consistent KKR calculation for: %s\n" % (self.jobname)
        line += "**********************************************************************\n"
        line += "SCFP:  information for self-consistency procedure:                   *\n"
        line += "**********************************************************************\n"
        line += "NITER.=%3i NLIN.=%3i NCPA.=%3i NPRN....=%08i\n" \
            % (self.niter, self.nlin, self.ncpa, self.nprn )
        # line += "NT...= %2i MNTA.= %2i\n" % (self.nt, self.mnta)
        line += "FRC...=  %1s DOS..=  %1s OPS..=  %1s " \
            % (self.frc, self.dos, self.ops) + \
            "AFM..=  " + self.afm + " CRT..=  " + self.crt + " STMP..= " + self.stmp + "\n"
        line += "Lmaxh.= %2i Lmaxt= %2i NFI..=%3i " \
            % (self.lmaxh, self.lmaxt, self.kgrn_nfi) \
            + "FIXG.= %2i SHF..=  %1i SOFC.=  %1s" \
            % (self.fixg, self.shf, self.sofc) + "\n"
        line += "KMSH...= %1s IBZ..= %2i NKX..= %2i " \
            % (self.kmsh, self.ibz, self.nkx) \
            + "NKY..= %2i NKZ..= %2i FBZ..=  %1s\n" \
            % (self.nky, self.nkz, self.fbz)
        #line += "KMSH2..= %1s IBZ2.=%3i NKX2.=%3i NKY2.=%3i NKZ2.=%3i NCPU.=%03i"\
        #    % (self.kmsh2, self.ibz2, self.nkx2, self.nky2, self.nkz2, self.ncpu) + "\n"
        line += "ZMSH...= %1s NZ1..= %2i " % (self.zmsh, self.nz1) \
            + "NZ2..=%3i NZ3..=%3i NRES.=%3i NZD..=%3i" \
            % (self.nz2, self.nz3, self.nres, self.nzd) + "\n"
        line += "DEPTH..= %6.3f IMAGZ.= %6.4f " \
            % (self.depth, self.imagz) \
            + "EPS...=  %5.3f ELIM..= %6.3f" \
            % (self.eps, self.elim) + "\n"
        line += "AMIX...= %6.3f VMIX..= %6.4f EFMIX.= %6.3f " % (self.amix,
                self.vmix, self.efmix) + "VMTZ..=%7.3f" \
            % (self.vmtz) + "\n"
        line += "TOLE...=%7.1e TOLEF.=%7.1e " % (self.tole, self.tolef) \
            + "TOLCPA=%7.1e TFERMI= %6.1f (K)\n" % (self.tolcpa, self.tfermi)
        line += "SWS....= %6.4f MMOM..= %6.3f\n" % (self.sws, self.mmom)
            # NSWS.=%3i " % (self.sws, self.nsws) \
            # + "DSWS..=   %4.2f ALPCPA= %6.4f\n" % (self.dsws, self.alpcpa)
        # line += "Setup: 2 + NQ*NS lines\n"
        line += "EFGS...= %6.3f HX....= %6.3f " % (self.efgs, self.hx) \
            + "NX...= %2i NZ0..= %2i" % (self.nx, self.nz0) \
            + " KPOLE=  " + str(self.kpole) + "\n"
        line += "**********************************************************************\n"
        line += "Sort:  information for alloy:                                        *\n"
        line += "******************************SS-screening**|***Magnetic structure ***\n"
        line += "Symb  IQ  IT ITA NRM  CONC      a_scr b_scr |Teta    Phi    FXM  m(split)\n"
        line += self.AtomOutput()
        line += "**********************************************************************\n"
        line += "Spin-spiral wave vector:\n"
        line += "qx....={0:9.6f} qy....={1:9.6f} qz....={2:9.6f}\n".format(
                self.qx, self.qy, self.qz)
        line += "**********************************************************************\n"
        line += "Atom:  information for atomic calculation:                           *\n"
        line += "**********************************************************************\n"
        line += "IEX...= %2i " % (self.iex) \
            + "NES..=%3i NITER=%3i IWAT.=%3i NPRNA=%3i" \
            % (self.nes, self.dirac_niter, self.iwat, self.nprna) + "\n"
        line += "VMIXATM..=  %8.6f RWAT....=  %8.6f RMAX....=%10.6f" \
            % (self.vmixatm, self.rwat, self.rmax) + "\n"
        line += "DX.......=  %8.6f DR1.....=  %8.6f TEST....=  %8.2e" \
            % (self.dx, self.dr1, self.test) + "\n"
        line += "TESTE....=  %8.2e TESTY...=  %8.2e TESTV...=  %8.2e" \
            % (self.teste, self.testy, self.testv) + "\n"

        # for i in self.confidx:
            # line += self.aconflines(self.conforder[i][0])

        return(line)

    def write_input_file(self, folder=None):
        """(self,str) ->(None)

            Save KGRN input data to file named filename

        :param folder:  (Default value = None)
        :type folder:
        :returns:  None
        :rtype:
        """

        # Check data integrity before anything is written on disk or run
        self.check_input_file()

        if folder is None:
            #sys.exit('Kgrn.create_input_file: \'folder\' has to be given!')
            folder = './'
        else:
            common.check_folders(folder)

        fl = open(folder + '/{0}.{1}'.format(self.jobname,self.KGRN_file_type), "w")
        fl.write(self.output())
        fl.close()

    def set_values(self, key, value):
        """


        :type key:
        :param value:, so we have to use .any()
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
                    '/kstr/' + self.latname + 'M.tfh'
                self.FOR001 = common.cleanup_path(self.FOR001)
                self.FOR001_2 = common.cleanup_path(self.FOR001_2)
                self.FOR002 = self.latpath + '/kstr/' + self.latname + '.mdl'
                self.FOR002 = common.cleanup_path(self.FOR002)

        else:
            setattr(self, key, value)

        return

    def create_atconf(self):
        """Constructs the self.atconf list out of the atomic information.

        :returns: None
        :rtype: None
        """

        self.atconf = []
        try:
            for i in range(len(self.atoms)):
                # self.atconf.append([self.atoms[i], self.iqs[i], self.its[i], self.itas[i],
                                    # self.concs[i], self.sm_ss[i],
                                    # self.s_wss[i], self.ws_wsts[i],
                                    # self.qtrs[i], self.splts[i], self.fxms[i]])
                self.atconf.append([self.atoms[i], self.iqs[i], self.its[i],
                    self.itas[i], self.nrms[i], self.concs[i],
                    self.a_scrs[i], self.b_scrs[i], self.tetas[i],
                    self.phis[i], self.fxms[i], self.splts[i]])
        except:
            pass


    def create_atomblock(self):
        """Constructs the KGRN input file atomblock if all
        the necessary parameters are present.

        :returns: None
        :rtype: None
        """

        for key in self.atconfKeys:
            # Check if we have a numpy array and do nothing
            if type(getattr(self, key)) == type(np.array([0.1,0.1])):
                continue

            # Some information needed to create the atomblock is still missing,
            # postpone further processing.
            elif getattr(self, key) == None:
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
            self.aconfig[i] = (self.atconf[i][0], self.atconf[i][1],
                                self.atconf[i][2], self.atconf[i][3],
                                self.atconf[i][4], self.atconf[i][5],
                               self.atconf[i][6], self.atconf[i][7],
                               self.atconf[i][8], self.atconf[i][9],
                               self.atconf[i][10], self.atconf[i][11])
            # self.aconfig[i] = self.atconf[i]

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
            self.atomblock.append(self.Atomline(self.aconfig[i]))
                # self.Atomline(self.aconfig[i][0], self.aconfig[i][1], self.aconfig[i][2],
                              # self.aconfig[i][3], self.aconfig[
                                  # i][4], self.aconfig[i][5],
                              # self.aconfig[i][6], self.aconfig[
                                  # i][7], self.aconfig[i][8],
                              # self.aconfig[i][9], self.aconfig[i][10],
                              # self.aconfig[i][11]))

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
            self.msgl = 1
        if self.expan is None:
            self.expan = 1
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
            for i in range(len(self.splts)):
                if not self.splts[i] == 0.0:
                    self.afm = 'F'
                    break
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
            self.kmsh = 'S'
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
            self.nzd = 100
        if self.depth is None:
            self.depth = 1.0
        if self.imagz is None:
            self.imagz = 0.020
        if self.eps is None:
            self.eps = 0.2
        if self.elim is None:
            self.elim = -1.0
        if self.amix is None:
            self.amix = 0.02
        if self.efmix is None:
            self.efmix = 0.9
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
            self.nz0 = self.nz1
        if self.stmp is None:
            self.stmp = 'N'
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
            self.vmix = 0.7
        if self.vmixatm is None:
            self.vmixatm = 0.3
        if self.rwat is None:
            self.rwat = 3.5
        if self.rmax is None:
            self.rmax = 20.0
        if self.dx is None:
            self.dx = 0.03
        if self.dr1 is None:
            self.dr1 = 0.001
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
        if self.FOR002 is None:
            self.FOR002 = self.latpath + '/bmdl/' + self.latname + '.mdl'
            self.FOR002 = common.cleanup_path(self.FOR002)
        if self.DIR006 is None:
            self.DIR006 = 'kgrn/'
        if self.DIR009 is None:
            self.DIR009 = 'kgrn/'
        if self.DIR010 is None:
            self.DIR010 = 'kgrn/'
        if self.DIR011 is None:
            self.DIR011 = 'kgrn/tmp/'
        if self.CQNA is None:
            self.CQNA = 'N'
        if self.ncpu is None:
            self.ncpu = 1
        if self.DIR022 is None:
            self.DIR022 = self.latpath + '/shape/' + self.latname + '.shp'
            self.DIR022 = common.cleanup_path(self.DIR022)
        if self.kpole is None:
            self.kpole = 0
        if self.qx is None:
            self.qx = 0.0
        if self.qy is None:
            self.qy = 0.0
        if self.qz is None:
            self.qz = 0.0
