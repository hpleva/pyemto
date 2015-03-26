# -*- coding: utf-8 -*-
"""

Created on Wed Dec  3 14:25:06 2014

@author: Matti Ropo
@author: Henrik Levämäki

"""

import os
import subprocess
import sys
import numpy as np


def run_emto(name, folder="./"):
    """Submits a batch script to the queue.

    Finds all the files in a folder that start with the name


    :param name:
    :type name:
    :param folder:  (Default value = "./")
    :type folder:
    :returns:
    :rtype:
    """

    # filetype = 'sh'
    filetype = 'cmd'
    namelist = []
    jobids = []
    for fltest in os.listdir(folder):
        if fltest.endswith('.' + filetype) and fltest.startswith(name):
            namelist.append(os.path.splitext(fltest)[0])
    for jobname in namelist:
        run_job = 'cd {0};sbatch {1}.{2}'.format(folder, jobname, filetype)
        jobid = 0
        jobid_raw = run_bash(run_job)
        print(jobid_raw)
        for t in jobid_raw.split():
            try:
                jobid = int(t)
            except ValueError:
                pass
        jobids.append(jobid)
    return jobids


def run_bash(cmd):
    """

    :param cmd:
    :type cmd:
    :returns:
    :rtype:
    """

    # Runs a bash command and returns the output
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    out = p.stdout.read().strip()
    return out  # This is the stdout from the shell command


def write_batch(folder, jobname):
    """

    :param folder:
    :type folder:
    :param jobname:
    :type jobname:
    :returns:
    :rtype:
    """

    # Writes a SLURM batch script for the job
    line = "#!/bin/bash" + "\n"
    line = line + "\n"
    line = line + "#SBATCH -J {0}".format(jobname) + "\n"
    line = line + "#SBATCH -t 01:00:00" + "\n"
    line = line + "#SBATCH -o {0}.output".format(jobname) + "\n"
    line = line + "#SBATCH -e {0}.error".format(jobname) + "\n"
    line = line + "\n"
    line = line + "python -u {0}.py".format(jobname) + "\n"

    fln = open("{0}/{1}.cmd".format(folder, jobname), "w")
    fln.write(line)
    fln.close()

    return


def write_job(folder, jobname, atom, conc, sws):
    """

    :param folder:
    :type folder:
    :param jobname:
    :type jobname:
    :param atom:
    :type atom:
    :param conc:
    :type conc:
    :param sws:
    :type sws:
    :returns:
    :rtype:
    """

    line = "import pyemto" + "\n"
    line += "\n"
    line += "test=pyemto.System()" + "\n"
    line += "test.bulk(latpath=\'../\',lat=\'fcc\',atoms={0}".format(atom)
    line += ",concs=[{1},{2}],sws={0},amix=0.02,efmix=0.9".format(sws, *conc)
    line += ",expan=\'M\',sofc=\'Y\')" + "\n"
    line += "test.latticeconstant(sws0={0})".format(sws) + "\n"

    fln = open("{0}/{1}.py".format(folder, jobname), "w")
    fln.write(line)
    fln.close()

    return


def submit_to_batch(folder, jobname, system='slurm'):
    """

    :param folder:
    :type folder:
    :param jobname:
    :type jobname:
    :param system:  (Default value = 'slurm')
    :type system:
    :returns:
    :rtype:
    """

    if system == 'slurm':
        command = "cd {0};sbatch {0}/{1}.cmd".format(folder, jobname)
        os.system(command)
    else:
        sys.exit('Only SLURM has been implemented so far! Exiting.')

    return


def wsrad_to_latparam(sws, lat, ca=None, c=None):
    """

    :param sws:
    :type sws:
    :param lat:
    :type lat:
    :param ca:  (Default value = None)
    :type ca:
    :param c:  (Default value = None)
    :type c:
    :returns:
    :rtype:
    """

    # Converts Wigner-Seitz radius (in bohr) to
    # a lattice parameter a (in Angstroms) (and
    # also c/a for hcp)
    bohr2angstrom = 0.52917721092
    if lat == 'sc':
        tmp = 4.0 * np.pi / 3.0
        a = tmp**(1.0 / 3.0) * sws * bohr2angstrom
        return a
    if lat == 'bcc':
        tmp = 4.0 * np.pi * 2.0 / 3.0
        a = tmp**(1.0 / 3.0) * sws * bohr2angstrom
        return a
    if lat == 'fcc':
        tmp = 4.0 * np.pi * 4.0 / 3.0
        a = tmp**(1.0 / 3.0) * sws * bohr2angstrom
        return a
    if lat == 'hcp':
        if c is None and ca is None:
            sys.exit('wsrad_to_latparam: Either \'c\' (in angstrom) or ' +\
                     '\'ca\' (c/a) has to be given for hcp structure!')
        else:
            if c is not None:
                a = np.sqrt(2.0 *
                            6.0 *
                            4.0 *
                            np.pi /
                            3.0 /
                            3.0 /
                            np.sqrt(3.0) /
                            c *
                            (sws *
                             bohr2angstrom)**3)
                return a
            else:
                vol = 4 * np.pi / 3.0 * sws**3 * 2
                a = (vol / ca / 0.8660254)**(1.0 / 3.0) * bohr2angstrom
                return a


def latparam_to_wsrad(a=None, lat=None, ca=None, c=None):
    """

    :param a:  (Default value = None)
    :type a:
    :param lat:  (Default value = None)
    :type lat:
    :param ca:  (Default value = None)
    :type ca:
    :param c:  (Default value = None)
    :type c:
    :returns:
    :rtype:
    """

    # Converts lattice parameters (in Angstroms) to Wigner-Seitz radius (in
    # bohr).
    bohr2angstrom = 0.52917721092
    if lat == 'sc':
        tmp = 4.0 * np.pi
        wsradius = (3.0 / tmp)**(1.0 / 3.0) * a / bohr2angstrom
        return wsradius
    if lat == 'bcc':
        tmp = 4.0 * np.pi * 2.0
        wsradius = (3.0 / tmp)**(1.0 / 3.0) * a / bohr2angstrom
        return wsradius
    if lat == 'fcc':
        tmp = 4.0 * np.pi * 4.0
        wsradius = (3.0 / tmp)**(1.0 / 3.0) * a / bohr2angstrom
        return wsradius
    if lat == 'hcp':
        if c is None and ca is None:
            sys.exit('latparam_to_wsrad: Either \'c\' (in angstrom) or ,' +\
                     '\'ca\' (c/a) has to be given for hcp structure!')
        else:
            if c is not None:
                wsradius = (3.0 * 3.0 * np.sqrt(3) * a**2 * c /
                            (2.0 * 6.0 * 4.0 * np.pi))**(1.0 / 3.0) / bohr2angstrom
                return wsradius
            else:
                wsradius = (3.0 * 3.0 * np.sqrt(3) * a**3 * ca /
                            (2.0 * 6.0 * 4.0 * np.pi))**(1.0 / 3.0) / bohr2angstrom
                return wsradius


def extrapolate_0k(
        a=None,
        Ta=None,
        bmod=None,
        TB=None,
        lat=None,
        TD=None,
        alpha=None,
        Talpha=None,
        B1=None,
        T0=None,
        ca=None,
        ZP=True):
    """Extrapolates experimental data to zero Kelvin.

    A function to extrapolate experimental room temperature lattice constants and
    bulk moduli to zero Kelvin. Zero point effects are also taken into account.

    +-----------------------+----------------------------------------------------------+
    | **Input parameter**   | **Description**                                          |
    +=======================+==========================================================+
    | a                     | Experimental lattice constant measured at                |
    |                       | temperature Ta (in Angstroms)                            |
    +-----------------------+----------------------------------------------------------+
    | Ta                    | The temperature at which a was measured (in Kelvin)      |
    +-----------------------+----------------------------------------------------------+
    | bmod                  | Experimental bulk modulus measured at temperature        |
    |                       | TB (in GPa)                                              |
    +-----------------------+----------------------------------------------------------+
    | TB                    | The temperature at which bmod was measured (in Kelvin)   |
    +-----------------------+----------------------------------------------------------+
    | lat                   | Lattice structure                                        |
    +-----------------------+----------------------------------------------------------+
    | TD                    | Experimental Debye temperature of the substance          |
    +-----------------------+----------------------------------------------------------+
    | alpha                 | VOLUMETRIC!!! thermal expansion coefficient at           |
    |                       | temperature T (in 10⁻6 1/K)                              |
    +-----------------------+----------------------------------------------------------+
    | Talpha                | The temperature at which alpha was measured (in Kelvin)  |
    +-----------------------+----------------------------------------------------------+
    | B1                    | Pressure derivative of the bulk modulus                  |
    +-----------------------+----------------------------------------------------------+
    | T0                    | The temperature to which we want to extrapolate          |
    |                       | (Optional argument, if not given assume 0K)              |
    +-----------------------+----------------------------------------------------------+
    | ca                    | c/a lattice parameter for hcp structures                 |
    +-----------------------+----------------------------------------------------------+
    | ZP                    | Whether or not zero-point (ZP) corrections are performed.|
    |                       | If one only wants to calculate the fractional volume     |
    |                       | change due to dropping temperature, ZP=False.            |
    +-----------------------+----------------------------------------------------------+

    :param a:  (Default value = None)
    :type a:
    :param Ta:  (Default value = None)
    :type Ta:
    :param bmod:  (Default value = None)
    :type bmod:
    :param TB:  (Default value = None)
    :type TB:
    :param lat:  (Default value = None)
    :type lat:
    :param TD:  (Default value = None)
    :type TD:
    :param alpha:  (Default value = None)
    :type alpha:
    :param Talpha:  (Default value = None)
    :type Talpha:
    :param B1:  (Default value = None)
    :type B1:
    :param T0:  (Default value = None)
    :type T0:
    :param ca:  (Default value = None)
    :type ca:
    :param ZP:  (Default value = True)
    :type ZP:
    :returns: ZPAE corrected 0K lattice constant (in Angstrom),
              ZPPE corrected 0K bulk modulus (in GPa)
    :rtype: float,float
    """

    import sys
    import numpy as np
    from scipy.integrate import quad, dblquad

    def debye_integrand(x):
        """Returns the value of the Debye function at a
        point x.

        :param x: Point on the x-axis where the function should
                  be evaluated
        :type x: float
        :returns: value of the function at point x
        :rtype: float
        """
        return (x**4 * np.e**x) / (np.e**x - 1)**2

    bohr = 0.52917721092  # Bohr to Angstrom conversion factor
    kb = 1.3806488E-23   # Boltzmann constant

    # Make sure all the input parameters are floats. This is because
    # if e.g. T and TD are given as integers you get integer division
    # (Python 2.7) which ruins the results: TD/Ta = 470/300 = 1.
    a = float(a)
    bmod = float(bmod)
    Ta = float(Ta)
    TB = float(TB)
    TD = float(TD)
    alpha = float(alpha)
    Talpha = float(Talpha)
    B1 = float(B1)
    if ca is not None:
        ca = float(ca)

    if lat is None:
        sys.exit('extrapolate_0K: \'lat\' has to be given!')
    if T0 is None:
        T0 = 2.0
    elif T0 < 2.0:
        # Stop a bit above 0K to prevent an overflow in the second integral.
        T0 = 2.0
    else:
        T0 = float(T0)
    if Ta is None:
        sys.exit('extrapolate_0K: \'Ta\' has to be given!')
    elif Ta < 2.0:
        # Stop a bit above 0K to prevent an overflow in the second integral.
        Ta = 2.0
    if TB is None:
        sys.exit('extrapolate_0K: \'TB\' has to be given!')
    elif TB < 2.0:
        # Stop a bit above 0K to prevent an overflow in the second integral.
        TB = 2.0
    #
    #
    ##########################################################
    #                                                        #
    # Compute the fractional volume change                   #
    # when temperature drops from T1 to T2                   #
    # as explained in:                                       #
    #                                                        #
    # Calphad                                                #
    # Volume 30, Issue 3, September 2006, Pages 354–356      #
    #                                                        #
    ##########################################################
    #
    #
    # Calculate the unit cell volume
    if lat in ['sc', 'bcc', 'fcc']:
        vol = a**3  # Holds true for cubic lattices

    elif lat == 'hcp':
        if ca is None:
            sys.exit(
                'extrapolate_0K: \'ca\' (c/a) has to be given for hcp structures!')
        else:
            vol = 3.0 * np.sqrt(3) * ca / 2 * (a * 1.0E-10)**3

    # Fractional Expansion Coefficient = FEC

    FECa = alpha * 1.0E-6 / Talpha**3 / quad(debye_integrand, 0,
                                             TD / Talpha)[0] * dblquad(lambda x,
                                                                       T: T**3 * (x**4 * np.e**x) /
                                                                       (np.e**x - 1)**2,
                                                                       Ta, T0,
                                                                       lambda x: 0,
                                                                       lambda x: TD / x)[0]

    FECB = alpha * 1.0E-6 / Talpha**3 / quad(debye_integrand, 0,
                                             TD / Talpha)[0] * dblquad(lambda x,
                                                                       T: T**3 * (x**4 * np.e**x) /
                                                                       (np.e**x - 1)**2,
                                                                       TB, T0,
                                                                       lambda x: 0,
                                                                       lambda x: TD / x)[0]

    # Thermal change in lattice parameter
    vol0 = (1.0 + FECa) * vol
    # Thermal change in bulk modulus
    P_T = bmod * FECB
    dB_thermal = B1 * P_T
    bmod0 = bmod - dB_thermal

    if lat in ['sc', 'bcc', 'fcc']:
        a0 = vol0**(1.0 / 3.0)

    elif lat == 'hcp':
        a0 = (2.0 * vol0 / (3.0 * np.sqrt(3.0) * ca))**(1.0 / 3.0)

    if not ZP:
        return a0, bmod0

    elif ZP:
        # Volume per atom is needed for the ZPAE !!!!!!
        if lat == 'sc':
            unit_vol0 = vol0
        elif lat == 'bcc':
            unit_vol0 = vol0 / 2.0
        elif lat == 'fcc':
            unit_vol0 = vol0 / 4.0
        elif lat == 'hcp':
            unit_vol0 = vol0 / 6.0
        #
        #
        #########################################################################
        #                                                                       #
        # Next, compute ZPAE for the lattice constant and ZPPE for bulk modulus #
        # as explained in PHYSICAL REVIEW B 66, 052104 (2002)                   #
        #                                                                       #
        #########################################################################
        #
        #
        factor = 9.0 / 8.0 * kb
        zeta = factor * TD

        P_Z = -1.0 / 6.0 * zeta * (B1 - 0.0) / unit_vol0
        dB_ZPPE = B1 * P_Z * 1.0E21
        bmod0_ZPPE = bmod0 - dB_ZPPE
        FECa_ZPAE = 1.0 / 6.0 * zeta * (B1 - 1.0) / (bmod0 * unit_vol0) * 1.0E21
        a0_ZPAE = (1.0 - FECa_ZPAE) * a0

        return a0_ZPAE, bmod0_ZPPE
