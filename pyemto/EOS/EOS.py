# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 10:51:00 2014

@author: Matti Ropo
@author: Henrik Levämäki

"""

from __future__ import print_function
import sys
import os
import numpy as np
import datetime
from scipy.optimize import leastsq


def _general_function(params, xdata, ydata, function):
    """This part comes from:

    http://projects.scipy.org/scipy/browser/trunk/scipy/optimize/minpack.py

    :param params:
    :type params:
    :param xdata:
    :type xdata:
    :param ydata:
    :type ydata:
    :param function:
    :type function:
    :returns:
    :rtype:
    """

    return function(xdata, *params) - ydata


def curve_fit(f, x, y, p0):
    """Fits an arbitraty function to data (x,y).

    This part comes from:
    http://projects.scipy.org/scipy/browser/trunk/scipy/optimize/minpack.py

    :param f:
    :type f:
    :param x:
    :type x:
    :param y:
    :type y:
    :param p0:
    :type p0:
    :returns:
    :rtype:
    """

    func = _general_function
    args = (x, y, f)
    popt, pcov, infodict, mesg, ier = leastsq(
        func, p0, args=args, full_output=1, maxfev=200000, ftol=1.0e-15, xtol=1.0e-15)

    if ier not in [1, 2, 3, 4]:
        raise RuntimeError("Optimal parameters not found: " + mesg)
    # end of this part
    return popt, pcov, infodict, mesg, ier


class EOS:
    """Fit equation of state for bulk systems.

    The following equations are available:

    ================ ==============================================================
    morse            PRB 37, 790 (1988)
    sjeos            PRB 63, 224115 (2001)
    taylor           A third order Taylor series expansion about the minimum volume
    murnaghan        PRB 28, 5480 (1983)
    birch            Intermetallic compounds: Principles and Practice,
                     Vol I: Principles. pages 195-210
    birchmurnaghan   PRB 70, 224107
    pouriertarantola PRB 70, 224107
    vinet            PRB 70, 224107
    antonschmidt     Intermetallics 11, 23-32 (2003)
    oldpoly          A third order polynomial fit (alternative implementation)
    ================ ==============================================================

    Use::

       eos = EquationOfState(volumes, energies, eos='murnaghan')
       v0, e0, B = eos.fit()
       eos.plot()

    :param name:
    :type name:
    :param xc:  (Default value = 'PBE')
    :type xc:
    :param method:  (Default value = 'morse')
    :type method:
    :param units:  (Default value = 'bohr')
    :type units:
    :returns:
    :rtype:
    """

    def __init__(self, name, xc='PBE', method='morse', units='bohr'):

        self.name = name
        self.xc = xc
        self.eos_string = method
        self.units = units
        self.error = 0.0
        self.relerror = 0.0
        self.residuals = []
        self.chisqr = 0.0
        self.redchi = 0.0
        self.rsquared = 0.0
        self.e = np.zeros(1)
        self.v = np.zeros(1)
        self.wsrad = 1.0
        self.nq = 1
        self.e0 = 0.0
        self.v0 = 0.0
        self.B = 0.0
        self.grun = 0.0
        self.bohr2a = 0.52917721092
        self.ry2ev = 13.605692533
        self.ev2gpa = 160.21773
        self.eMin = 0.0
        self.eos_parameters = None
        self.sjeosfit0 = None

    def vol2wsrad(self, V):
        """

        :param V:
        :type V:
        :returns:
        :rtype:
        """
        return (3 * V / (4 * np.pi))**(1.0 / 3.0)

    def wsrad2vol(self, WSrad):
        """

        :param WSrad:
        :type WSrad:
        :returns:
        :rtype:
        """
        #print('debug:', WSrad)
        return 4.0 * np.pi / 3.0 * WSrad**3

    def angstrom2bohr(self, V):
        """

        :param V:
        :type V:
        :returns:
        :rtype:
        """
        return (3 * V / (4 * np.pi))**(1.0 / 3.0) / 0.52917721092

    def bohr2angstrom(self, WSrad):
        """

        :param WSrad:
        :type WSrad:
        :returns:
        :rtype:
        """
        return 4.0 * np.pi / 3.0 * (WSrad * 0.52917721092)**3

    def morse(self, w, a0, b0, c0, l0):
        """

        :param w:
        :type w:
        :param a0:
        :type a0:
        :param b0:
        :type b0:
        :param c0:
        :type c0:
        :param l0:
        :type l0:

        """
        E = a0 + b0 * np.exp(-l0 * w) + c0 * np.exp(-2.0 * l0 * w)
        return E

    def taylor(self, V, E0, beta, alpha, V0):
        """Taylor Expansion up to 3rd order about V0

        :param V:
        :type V:
        :param E0:
        :type E0:
        :param beta:
        :type beta:
        :param alpha:
        :type alpha:
        :param V0:
        :type V0:
        :returns:
        :rtype:
        """

        E = E0 + beta / 2. * (V - V0)**2 / V0 + alpha / 6. * (V - V0)**3 / V0
        return E

    def murnaghan(self, V, E0, B0, BP, V0):
        """From PRB 28,5480 (1983)

        :param V:
        :type V:
        :param E0:
        :type E0:
        :param B0:
        :type B0:
        :param BP:
        :type BP:
        :param V0:
        :type V0:
        :returns:
        :rtype:
        """

        E = E0 + B0 * V / BP * \
            (((V0 / V)**BP) / (BP - 1) + 1) - V0 * B0 / (BP - 1)
        return E

    def birch(self, V, E0, B0, BP, V0):
        """From Intermetallic compounds: Principles and Practice, Vol. I: Principles
        Chapter 9 pages 195-210 by M. Mehl. B. Klein, D. Papaconstantopoulos
        paper downloaded from Web

        case where n=0

        :param V:
        :type V:
        :param E0:
        :type E0:
        :param B0:
        :type B0:
        :param BP:
        :type BP:
        :param V0:
        :type V0:
        :returns:
        :rtype:
        """

        E = (E0 + 9.0 / 8.0 * B0 * V0 * ((V0 / V)**(2.0 / 3.0) - 1.0)**2 +
             9.0 / 16.0 * B0 * V0 * (BP - 4.) * ((V0 / V)**(2.0 / 3.0) - 1.0)**3)
        return E

    def birchmurnaghan(self, V, E0, B0, BP, V0):
        """
        BirchMurnaghan equation from PRB 70, 224107

        :param V:
        :type V:
        :param E0:
        :type E0:
        :param B0:
        :type B0:
        :param BP:
        :type BP:
        :param V0:
        :type V0:
        :returns:
        :rtype:
        """

        eta = (V / V0)**(1. / 3.)
        E = E0 + 9. * B0 * V0 / 16. * \
            (eta**2 - 1)**2 * (6 + BP * (eta**2 - 1.) - 4. * eta**2)
        return E

    def pouriertarantola(self, V, E0, B0, BP, V0):
        """
        Pourier-Tarantola equation from PRB 70, 224107

        :param V:
        :type V:
        :param E0:
        :type E0:
        :param B0:
        :type B0:
        :param BP:
        :type BP:
        :param V0:
        :type V0:
        :returns:
        :rtype:
        """

        eta = (V / V0)**(1. / 3.)
        squiggle = -3. * np.log(eta)

        E = E0 + B0 * V0 * squiggle**2 / 6. * (3. + squiggle * (BP - 2))
        return E

    def vinet(self, V, E0, B0, BP, V0):
        """
        Vinet equation from PRB 70, 224107

        :param V:
        :type V:
        :param E0:
        :type E0:
        :param B0:
        :type B0:
        :param BP:
        :type BP:
        :param V0:
        :type V0:
        :returns:
        :rtype:
        """

        eta = (V / V0)**(1. / 3.)

        E = (E0 + 2. * B0 * V0 / (BP - 1.)**2 * (
            2. - (5. + 3. * BP * (eta - 1.) - 3. * eta) * np.exp(-3. * (BP - 1.) * (eta - 1.) / 2.))
            )
        return E

    def antonschmidt(self, V, Einf, B, n, V0):
        """From Intermetallics 11, 23-32 (2003)

        Einf should be E_infinity, i.e. infinite separation, but
        according to the paper it does not provide a good estimate
        of the cohesive energy. They derive this equation from an
        empirical formula for the volume dependence of pressure,

        E(vol) = E_inf + int(P dV) from V=vol to V=infinity

        but the equation breaks down at large volumes, so E_inf
        is not that meaningful

        n should be about -2 according to the paper.

        I find this equation does not fit volumetric data as well
        as the other equtions do.

        :param V:
        :type V:
        :param Einf:
        :type Einf:
        :param B:
        :type B:
        :param n:
        :type n:
        :param V0:
        :type V0:
        :returns: Energy
        :rtype: float
        """

        E = B * V0 / (n + 1.) * (V / V0)**(n + 1.) * \
            (np.log(V / V0) - (1. / (n + 1.))) + Einf
        return E

    def oldpoly(self, V, c0, c1, c2, c3):
        """polynomial fit, 3rd order

        :param V:
        :type V:
        :param c0:
        :type c0:
        :param c1:
        :type c1:
        :param c2:
        :type c2:
        :param c3:
        :type c3:
        :returns: Energy
        :rtype: float
        """

        E = c0 + c1 * V + c2 * V**2 + c3 * V**3
        return E

    def parabola(self, x, a, b, c):
        """Parabola polynomial function

        This function is used to fit the data to get good guesses for
        the equation of state fits.

        A 4th order polynomial fit to get good guesses for
        was not a good idea because for noisy data the fit is too wiggly
        2nd order seems to be sufficient, and guarentees a single minimum.

        :param x:
        :type x:
        :param a:
        :type a:
        :param b:
        :type b:
        :param c:
        :type c:
        :returns: Energy
        :rtype: float
        """

        return a + b * x + c * x**2

    def predicted(self):
        """Evaluates the EOS function using the calc. EOS params.

        :returns: EOS function
        :rtype: func
        """

        fity = eval('self.{0}'.format(self.eos_string))(
            self.v,
            self.eos_parameters[0],
            self.eos_parameters[1],
            self.eos_parameters[2],
            self.eos_parameters[3],)
        return fity

    def prepareData(self):
        """

        :returns:
        :rtype:
        """

        fitpts = open('fit/' + self.name + '.fit', 'w')
        now = datetime.datetime.now()
        timeline = str(now.day) + "." + str(now.month) + "." + str(now.year) +\
            " -- " + "{0:02}".format(now.hour) + ":" +\
            "{0:02}".format(now.minute) + ":" + \
            "{0:02}".format(now.second) + "\n"
        fitpts.write(
            timeline +
            "JOBNAM = " +
            self.name +
            ' -- ' +
            self.xc +
            "\n")

        namelist = []
        for file in os.listdir("kfcd/"):
            if self.name in file:
                namelist.append(file)
        data = np.array([[0, 0]])
        enTag = 'TOT-' + self.xc
        volTag = 'VOL'
        nqTag = 'NQ'
        nq = 1
        for n in namelist:
            nqRead = False
            lineIndex = 0
            fl = open("kfcd/" + n)
            lines = fl.readlines()
            for line in lines:
                if not nqRead:
                    if nqTag in line:
                        nq = int(line.split()[2])
                        nqRead = True
                """
                if volTag in line:
                    # cellVol = volume of the cell in unit
                    # of cellA**3
                    cellVol = float(line.split()[5])
                    #cellA = length of the lattice parameter A
                    # in unit of Bohr
                    cellA = float(lines[lineIndex+2].split()[2])
                    #
                    cellVol = cellVol*cellA**3
                    cellWSRad = (3.0*cellVol/(4*np.pi))**(1.0/3.0)
                """
                if enTag in line:
                    linesplit = line.split()
                    data = np.append(data, [[float(linesplit[6]),
                                             float(linesplit[3])]], axis=0)
                lineIndex += 1
            fl.close()

        data = np.delete(data, 0, 0)
        # Make sure datapoints are in ascending order
        data = data[data[:, 0].argsort()]

        fitpts.write('NQ = {0}'.format(nq) + "\n")
        fitpts.write("Sws.......Etot..........." + "\n")
        for i in range(len(data[:, 0])):
            fitpts.write('{0:.6f} {1:.6f}\n'.format(data[i, 0], data[i, 1]))
        fitpts.close()

    def fit(self, swses, energies):
        """Calculate volume, energy, and bulk modulus.

        A fitting input file is created in the '/fit' folder.
        Returns the optimal volume, the minumum energy, and the bulk
        modulus.  Notice that the units for the bulk modulus is
        eV/Angstrom^3 - to get the value in GPa, do this::

        v0, e0, B = eos.fit()
        print B*160.21773, 'GPa'

        :param swses: List of WS-radii
        :type swses: list(float)
        :param energies: List of energies
        :type energies: list(float)
        :returns: Eq. WS-rad, energy, bulk modulus, Gruneisen parameter
        :rtype: float, float, float, float
        """

        self.v = np.asarray(swses)
        self.e = np.asarray(energies)
        # Remove duplicates and make sure the arrays are sorted
        vLen0 = len(self.v)
        self.v, uniqueInd = np.unique(self.v, return_index=True)
        self.e = np.array(self.e)[uniqueInd]
        if len(self.v) < vLen0:
            print('EOS.fit(): WARNING! Duplicate volumes have been detected.' + '\n' +\
                  'Duplicates have been removed in no particular fashion.' +
                  '\n' +\
                  'Check the integrity of the data!' + '\n')
        if self.eos_string != 'morse':
            # Convert WS-radius to Angstrom^3
            self.v[:] = 4.0 * np.pi / 3.0 * (self.bohr2a * self.v[:])**3
            # Convert Rydberg to eV (No need to do this unless one wants to
            # print out in eVs)
            # self.e[:] *= self.ry2ev
        # Translate energies so that smallest value becomes zero
        self.eMin = np.min(self.e)
        self.e -= self.eMin

        if self.eos_string == 'sjeos':
            fit0 = np.poly1d(np.polyfit(self.v**-(1.0 / 3), self.e, 3))
            fit1 = np.polyder(fit0, 1)
            fit2 = np.polyder(fit1, 1)
            self.v0 = None
            for t in np.roots(fit1):
                if isinstance(t, float) and t > 0 and fit2(t) > 0:
                    self.v0 = t**-3
                    break
            if self.v0 is None:
                raise ValueError('SJEOS: No minimum!')
            self.e0 = fit0(t) + self.eMin
            self.B = t**5 * fit2(t) / 9
            self.sjeosfit0 = fit0
            # Convert B to GPa
            self.B *= self.ry2ev * self.ev2gpa

        else:
            p0 = [min(self.e), 1, 1]
            if self.eos_string == 'morse':
                popt, pcov, infodict0, mesg0, ier0 = curve_fit(
                    self.parabola, self.wsrad2vol(
                        self.v), self.e, p0)
            else:
                popt, pcov, infodict0, mesg0, ier0 = curve_fit(
                    self.parabola, self.v, self.e, p0)

            parabola_parameters = popt
            # Here I just make sure the minimum is bracketed by the volumes
            # this if for the solver
            minvol = min(self.v)
            maxvol = max(self.v)

            # the minimum of the parabola is at dE/dV = 0, or 2*c V +b =0
            c = parabola_parameters[2]
            b = parabola_parameters[1]
            a = parabola_parameters[0]
            parabola_vmin = -b / 2 / c
            if self.eos_string == 'morse':
                parabola_vmin = self.vol2wsrad(parabola_vmin)
            # print('parabola_vmin=',parabola_vmin)

            if not (minvol < parabola_vmin and parabola_vmin < maxvol):
                print('EOS.fit(): Warning the minimum volume of a fitted' +\
                      ' parabola is not in your volumes.\n' +\
                      '           You may not have a minimum in your dataset.')

            # evaluate the parabola at the minimum to estimate the groundstate
            # energy
            E0 = self.parabola(parabola_vmin, a, b, c)
            # print('parabola_E0=',E0)
            # estimate the bulk modulus from Vo*E''.  E'' = 2*c
            if self.eos_string == 'morse':
                B0 = 2 * c * self.wsrad2vol(parabola_vmin)
            else:
                B0 = 2 * c * parabola_vmin
            # print('parabola_B0=',B0)

            if self.eos_string == 'antonschmidt':
                BP = -2
            else:
                BP = 4

            if self.eos_string == 'morse':
                w0 = parabola_vmin
                l0 = 2.0 * (2.0 - 1.0 / 3.0) / w0
                x0 = np.exp(-l0 * w0)
                c0 = -6.0 * np.pi * np.log(x0) * B0 / (x0**2 * l0**3)
                b0 = -2.0 * c0 * x0
                a0 = E0 - b0 * x0 - c0 * x0**2
                initial_guess = [a0, b0, c0, l0]
                # print('a0=',a0,'b0=',b0,'c0=',c0,'lambda0=',l0)
            else:
                initial_guess = [E0, B0, BP, parabola_vmin]

            # now fit the equation of state
            p0 = initial_guess
            popt, pcov, infodict, mesg, ier = curve_fit(
                eval('self.{0}'.format(self.eos_string)), self.v, self.e, p0)

            self.eos_parameters = popt

            self.e += self.eMin

            if self.eos_string == 'morse':
                ma, mb, mc, ml = self.eos_parameters
                # print('morse_a=',ma,'morse_b=',mb,'morse_c=',mc,'morse_lambda=',ml)
                mx0 = -0.5 * mb / mc
                # print('morse_x0=',mx0)
                self.v0 = -np.log(mx0) / ml
                self.e0 = ma + mb * mx0 + mc * mx0**2 + self.eMin
                self.B = -(mc * mx0**2 * ml**3) / (6 * np.pi * np.log(mx0))
                # Convert B to GPa
                self.B = self.ev2gpa * self.ry2ev / self.bohr2a**3 * self.B
                self.grun = 0.5 * ml * self.v0
            elif self.eos_string == 'oldpoly':
                c0, c1, c2, c3 = self.eos_parameters
                # find minimum E in E = c0 + c1*V + c2*V**2 + c3*V**3
                # dE/dV = c1+ 2*c2*V + 3*c3*V**2 = 0
                # solve by quadratic formula with the positive root

                a = 3 * c3
                b = 2 * c2
                c = c1

                self.v0 = (-b + np.sqrt(b**2 - 4 * a * c)) / (2 * a)
                self.e0 = self.oldpoly(self.v0, c0, c1, c2, c3) + self.eMin
                # Convert e0 to Rydbergs
                #self.e0 /= self.ry2ev
                self.B = (2 * c2 + 6 * c3 * self.v0) * self.v0
                # Convert B to GPa
                self.B *= self.ry2ev * self.ev2gpa
                # Gruneisen parameter for oldpoly not implemented
                self.grun = 0.0
            else:
                self.e0 = self.eos_parameters[0] + self.eMin
                self.v0 = self.eos_parameters[3]
                self.B = self.eos_parameters[1]
                self.grun = -0.5 + self.eos_parameters[2] / 2.0
                # Convert e0 to Rydbergs
                #self.e0 /= self.ry2ev
                # Convert B to GPa
                self.B *= self.ry2ev * self.ev2gpa

        if self.units == 'bohr':
            if self.eos_string == 'morse':
                pass
            else:
                self.v0 = self.angstrom2bohr(self.v0)

        if self.units == 'angstrom':
            if self.eos_string == 'morse':
                self.v0 = self.bohr2angstrom(self.v0)
            else:
                pass

        if self.eos_string == 'sjeos':
            sjeospredicted = self.sjeosfit0(self.v**-(1.0 / 3.0))
            self.residuals = self.e - sjeospredicted
            self.chisqr = np.sum(self.residuals**2)
            self.redchi = self.chisqr / (len(self.e) - 4)
            ss_tot = np.sum((self.e - np.mean(self.e))**2)
            self.rsquared = 1.0 - self.chisqr / ss_tot
            fitE = sjeospredicted + self.eMin
            self.e += self.eMin
        else:
            self.residuals = infodict['fvec']
            self.chisqr = np.sum(self.residuals**2)
            self.redchi = self.chisqr / (len(self.e) - 4)
            ss_tot = np.sum((self.e - np.mean(self.e))**2)
            self.rsquared = 1.0 - self.chisqr / ss_tot
            fitE = self.predicted() + self.eMin

        # Calculate WS-radius
        self.wsrad = self.v0

        # Printouts
        now = datetime.datetime.now()
        timeline = str(now.day) + "." + str(now.month) + "." + str(now.year) +\
            " -- " + "{0:02}".format(now.hour) + ":" +\
            "{0:02}".format(now.minute) + ":" + \
            "{0:02}".format(now.second) + "\n"
        print(timeline + "JOBNAM = " + self.name + ' -- ' + self.xc + "\n")
        print('Using ' + self.eos_string + ' function' + "\n")
        print('Chi squared         = ' + str(self.chisqr))
        print('Reduced Chi squared = ' + str(self.redchi))
        print('R squared           = ' + str(self.rsquared) + "\n")
        print(self.eos_string + ' parameters:' + "\n")
        if self.eos_string == 'morse':
            print('a      = {0:12.6f}'.format(self.eos_parameters[0]))
            print('b      = {0:12.6f}'.format(self.eos_parameters[1]))
            print('c      = {0:12.6f}'.format(self.eos_parameters[2]))
            print('lambda = {0:12.6f}'.format(self.eos_parameters[3]) + "\n")
        elif self.eos_string == 'sjeos':
            print('a = {0:12.6f}'.format(self.sjeosfit0.c[0]))
            print('b = {0:12.6f}'.format(self.sjeosfit0.c[1]))
            print('c = {0:12.6f}'.format(self.sjeosfit0.c[2]))
            print('d = {0:12.6f}'.format(self.sjeosfit0.c[3]) + "\n")
        elif self.eos_string == 'oldpoly':
            print('c0 = {0:12.6f}'.format(self.eos_parameters[0]))
            print('c1 = {0:12.6f}'.format(self.eos_parameters[1]))
            print('c2 = {0:12.6f}'.format(self.eos_parameters[2]))
            print('c3 = {0:12.6f}'.format(self.eos_parameters[3]) + "\n")
        else:
            print('E0    = {0:12.6f}'.format(self.eos_parameters[0]))
            print('Bmod  = {0:12.6f}'.format(self.eos_parameters[1]))
            print('Bmod\' = {0:12.6f}'.format(self.eos_parameters[2]))
            print('V0    = {0:12.6f}'.format(self.eos_parameters[3]) + "\n")
        print('Ground state parameters:' + "\n")
        if self.units == 'bohr':
            print('V0           = {0:12.6f} Bohr^3 (unit cell volume)'
                  .format(self.v0))
            print('             = {0:12.6f} Bohr   (WS-radius)'
                  .format(self.wsrad))
        elif self.units == 'angstrom':
            print('V0           = {0:12.6f} Angstrom^3'.format(self.v0))
        print('E0           = {0:12.6f} Ry'.format(self.e0 * self.nq))
        print('Bmod         = {0:12.6f} GPa'.format(self.B))
        print('Grun. param. = {0:12.6f}'.format(self.grun) + "\n")
        print('sws            Einp           Eout          Residual       ' +
              'err (% * 10**6)')
        for i in range(len(self.e)):
            print('{0:.6f} {1:13.6f} {2:13.6f} {3:13.6f} {4:15.6f}'
                  .format(self.v[i], self.e[i], fitE[i], self.residuals[i],
                          self.residuals[i] / np.abs(self.e[i]) * 1.0E6))
        print('\n')

        # Save an image for quality checking
        #self.plot(filename=self.name, show=None)

        return self.wsrad, self.e0, self.B, self.grun  # , self.error

    def fit2file(self):
        """Calculate volume, energy, and bulk modulus.

        Returns the optimal volume, the minumum energy, and the bulk
        modulus.  Notice that the units for the bulk modulus is
        eV/Angstrom^3 - to get the value in GPa, do this::

          v0, e0, B = eos.fit()
          print B*160.21773, 'GPa'

        :returns:
        :rtype:
        """

        # Read data from the fit folder
        file = open('fit/' + self.name + '.fit', mode='r')
        lines = file.readlines()
        self.nq = int(lines[2].split()[2])
        file.close()
        data = np.loadtxt('fit/' + self.name + '.fit', skiprows=4)
        self.v = data[:, 0]
        self.e = data[:, 1]
        if self.eos_string != 'morse':
            # Convert WS-radius to Angstrom^3
            self.v[:] = 4.0 * np.pi / 3.0 * (self.bohr2a * self.v[:])**3
            # Convert Rydberg to eV (No need to do this unless one wants to
            # print out in eVs)
            # self.e[:] *= self.ry2ev
        # Translate energies so that smallest value becomes zero
        self.eMin = np.min(self.e)
        self.e -= self.eMin

        if self.eos_string == 'sjeos':
            fit0 = np.poly1d(np.polyfit(self.v**-(1.0 / 3), self.e, 3))
            fit1 = np.polyder(fit0, 1)
            fit2 = np.polyder(fit1, 1)
            self.v0 = None
            for t in np.roots(fit1):
                if isinstance(t, float) and t > 0 and fit2(t) > 0:
                    self.v0 = t**-3
                    break
            if self.v0 is None:
                raise ValueError('sjeos: No minimum!')
            self.e0 = fit0(t) + self.eMin
            self.B = t**5 * fit2(t) / 9
            self.sjeosfit0 = fit0
            # Convert B to GPa
            self.B *= self.ry2ev * self.ev2gpa

        else:
            p0 = [min(self.e), 1, 1]
            if self.eos_string == 'morse':
                popt, pcov, infodict0, mesg0, ier0 = curve_fit(
                    self.parabola, self.wsrad2vol(
                        self.v), self.e, p0)
            else:
                popt, pcov, infodict0, mesg0, ier0 = curve_fit(
                    self.parabola, self.v, self.e, p0)

            parabola_parameters = popt
            # Here I just make sure the minimum is bracketed by the volumes
            # this if for the solver
            minvol = min(self.v)
            maxvol = max(self.v)

            # the minimum of the parabola is at dE/dV = 0, or 2*c V +b =0
            c = parabola_parameters[2]
            b = parabola_parameters[1]
            a = parabola_parameters[0]
            parabola_vmin = -b / 2 / c
            if self.eos_string == 'morse':
                parabola_vmin = self.vol2wsrad(parabola_vmin)
            # print('parabola_vmin=',parabola_vmin)

            if not (minvol < parabola_vmin and parabola_vmin < maxvol):
                print('EOS.fit2file(): Warning the minimum volume of a fitted' +\
                      ' parabola is not in your volumes.\n' +\
                      '                You may not have a minimum in your dataset.')

            # evaluate the parabola at the minimum to estimate the groundstate
            # energy
            E0 = self.parabola(parabola_vmin, a, b, c)
            # print('parabola_E0=',E0)
            # estimate the bulk modulus from Vo*E''.  E'' = 2*c
            if self.eos_string == 'morse':
                B0 = 2 * c * self.wsrad2vol(parabola_vmin)
            else:
                B0 = 2 * c * parabola_vmin
            # print('parabola_B0=',B0)

            if self.eos_string == 'antonschmidt':
                BP = -2
            else:
                BP = 4

            if self.eos_string == 'morse':
                w0 = parabola_vmin
                l0 = 2.0 * (2.0 - 1.0 / 3.0) / w0
                x0 = np.exp(-l0 * w0)
                c0 = -6.0 * np.pi * np.log(x0) * B0 / (x0**2 * l0**3)
                b0 = -2.0 * c0 * x0
                a0 = E0 - b0 * x0 - c0 * x0**2
                initial_guess = [a0, b0, c0, l0]
                # print('a0=',a0,'b0=',b0,'c0=',c0,'lambda0=',l0)
            else:
                initial_guess = [E0, B0, BP, parabola_vmin]

            # now fit the equation of state
            p0 = initial_guess
            popt, pcov, infodict, mesg, ier = curve_fit(
                eval(
                    'self.{0}'.format(
                        self.eos_string)), self.v, self.e, p0)

            self.eos_parameters = popt

            self.e += self.eMin

            if self.eos_string == 'morse':
                ma, mb, mc, ml = self.eos_parameters
                # print('morse_a=',ma,'morse_b=',mb,'morse_c=',mc,'morse_lambda=',ml)
                mx0 = -0.5 * mb / mc
                # print('morse_x0=',mx0)
                self.v0 = -np.log(mx0) / ml
                self.e0 = ma + mb * mx0 + mc * mx0**2 + self.eMin
                self.B = -(mc * mx0**2 * ml**3) / (6 * np.pi * np.log(mx0))
                # Convert B to GPa
                self.B = self.ev2gpa * self.ry2ev / self.bohr2a**3 * self.B
                self.grun = 0.5 * ml * self.v0
            elif self.eos_string == 'oldpoly':
                c0, c1, c2, c3 = self.eos_parameters
                # find minimum E in E = c0 + c1*V + c2*V**2 + c3*V**3
                # dE/dV = c1+ 2*c2*V + 3*c3*V**2 = 0
                # solve by quadratic formula with the positive root

                a = 3 * c3
                b = 2 * c2
                c = c1

                self.v0 = (-b + np.sqrt(b**2 - 4 * a * c)) / (2 * a)
                self.e0 = self.oldpoly(self.v0, c0, c1, c2, c3) + self.eMin
                # Convert e0 to Rydbergs
                #self.e0 /= self.ry2ev
                self.B = (2 * c2 + 6 * c3 * self.v0) * self.v0
                # Convert B to GPa
                self.B *= self.ry2ev * self.ev2gpa
                # Gruneisen parameter for oldpoly not implemented
                self.grun = 0.0
            else:
                self.e0 = self.eos_parameters[0] + self.eMin
                self.v0 = self.eos_parameters[3]
                self.B = self.eos_parameters[1]
                self.grun = -0.5 + self.eos_parameters[2] / 2.0
                # Convert e0 to Rydbergs
                #self.e0 /= self.ry2ev
                # Convert B to GPa
                self.B *= self.ry2ev * self.ev2gpa

        if self.units == 'bohr':
            if self.eos_string == 'morse':
                pass
            else:
                self.v0 = self.angstrom2bohr(self.v0)

        if self.units == 'angstrom':
            if self.eos_string == 'morse':
                self.v0 = self.bohr2angstrom(self.v0)
            else:
                pass

        if self.eos_string == 'sjeos':
            sjeospredicted = self.sjeosfit0(self.v**-(1.0 / 3.0))
            self.residuals = self.e - sjeospredicted
            self.chisqr = np.sum(self.residuals**2)
            self.redchi = self.chisqr / (len(self.e) - 4)
            ss_tot = np.sum((self.e - np.mean(self.e))**2)
            self.rsquared = 1.0 - self.chisqr / ss_tot
            fitE = sjeospredicted + self.eMin
            self.e += self.eMin
        else:
            self.residuals = infodict['fvec']
            self.chisqr = np.sum(self.residuals**2)
            self.redchi = self.chisqr / (len(self.e) - 4)
            ss_tot = np.sum((self.e - np.mean(self.e))**2)
            self.rsquared = 1.0 - self.chisqr / ss_tot
            fitE = self.predicted() + self.eMin

        # Calculate WS-radius
        self.wsrad = self.v0
        # Calculate unit cell volume
        self.v0 = 4.0 / 3.0 * np.pi * self.v0**3.0 * self.nq

        # Write results to a text file
        fitresult = open('fit/' + self.name + '.prn', 'w')
        now = datetime.datetime.now()
        timeline = str(now.day) + "." + str(now.month) + "." + str(now.year) +\
            " -- " + "{0:02}".format(now.hour) + ":" +\
            "{0:02}".format(now.minute) + ":" + \
            "{0:02}".format(now.second) + "\n"
        fitresult.write(
            timeline +
            "JOBNAM = " +
            self.name +
            ' -- ' +
            self.xc +
            "\n" +
            "\n")
        fitresult.write('Using ' + self.eos_string + ' function' + "\n" + "\n")
        fitresult.write('Chi squared         = ' + str(self.chisqr) + "\n")
        fitresult.write('Reduced Chi squared = ' + str(self.redchi) + "\n")
        fitresult.write('R squared           = ' +
                        str(self.rsquared) +
                        "\n" +
                        "\n")
        fitresult.write(self.eos_string + ' parameters:' + "\n" + "\n")
        if self.eos_string == 'morse':
            fitresult.write(
                'a      = {0:12.6f}'.format(
                    self.eos_parameters[0]) +
                "\n")
            fitresult.write(
                'b      = {0:12.6f}'.format(
                    self.eos_parameters[1]) +
                "\n")
            fitresult.write(
                'c      = {0:12.6f}'.format(
                    self.eos_parameters[2]) +
                "\n")
            fitresult.write(
                'lambda = {0:12.6f}'.format(
                    self.eos_parameters[3]) +
                "\n" +
                "\n")
        elif self.eos_string == 'sjeos':
            fitresult.write('a = {0:12.6f}'.format(self.sjeosfit0.c[0]) + "\n")
            fitresult.write('b = {0:12.6f}'.format(self.sjeosfit0.c[1]) + "\n")
            fitresult.write('c = {0:12.6f}'.format(self.sjeosfit0.c[2]) + "\n")
            fitresult.write(
                'd = {0:12.6f}'.format(
                    self.sjeosfit0.c[3]) +
                "\n" +
                "\n")
        elif self.eos_string == 'oldpoly':
            fitresult.write(
                'c0 = {0:12.6f}'.format(
                    self.eos_parameters[0]) +
                "\n")
            fitresult.write(
                'c1 = {0:12.6f}'.format(
                    self.eos_parameters[1]) +
                "\n")
            fitresult.write(
                'c2 = {0:12.6f}'.format(
                    self.eos_parameters[2]) +
                "\n")
            fitresult.write(
                'c3 = {0:12.6f}'.format(
                    self.eos_parameters[3]) +
                "\n" +
                "\n")
        else:
            fitresult.write(
                'E0    = {0:12.6f}'.format(
                    self.eos_parameters[0]) +
                "\n")
            fitresult.write(
                'Bmod  = {0:12.6f}'.format(
                    self.eos_parameters[1]) +
                "\n")
            fitresult.write(
                'Bmod\' = {0:12.6f}'.format(
                    self.eos_parameters[2]) +
                "\n")
            fitresult.write(
                'V0    = {0:12.6f}'.format(
                    self.eos_parameters[3]) +
                "\n" +
                "\n")
        fitresult.write('Ground state parameters:' + "\n" + "\n")
        if self.units == 'bohr':
            fitresult.write(
                'V0           = {0:12.6f} Bohr^3 (unit cell volume)' .format(
                    self.v0) +
                "\n")
            fitresult.write('             = {0:12.6f} Bohr   (WS-radius)'
                            .format(self.wsrad) + "\n")
        elif self.units == 'angstrom':
            fitresult.write(
                'V0           = {0:12.6f} Angstrom^3'.format(
                    self.v0) +
                "\n")
        fitresult.write(
            'E0           = {0:12.6f} Ry'.format(
                self.e0 *
                self.nq) +
            "\n")
        fitresult.write('Bmod         = {0:12.6f} GPa'.format(self.B) + "\n")
        fitresult.write(
            'Grun. param. = {0:12.6f}'.format(
                self.grun) +
            "\n" +
            "\n")
        fitresult.write(
            'sws            Einp           Eout          Residual       ' +
            'err (% * 10**6)' +
            "\n")
        for i in range(len(self.e)):
            fitresult.write(
                '{0:.6f} {1:13.6f} {2:13.6f} {3:13.6f} {4:15.6f}' .format(
                    self.v[i],
                    self.e[i],
                    fitE[i],
                    self.residuals[i],
                    self.residuals[i] /
                    np.abs(
                        self.e[i]) *
                    1.0E6) +
                "\n")
        fitresult.close()

        # Save an image for quality checking
        #self.plot(filename=self.name, show=None)

        return self.wsrad, self.e0, self.B, self.grun  # , self.error

    def plot(self, filename=None, show=None):
        """Plot fitted energy curve.

        Uses Matplotlib to plot the energy curve.  Use *show=True* to
        show the figure and *filename='abc.png'* or
        *filename='abc.eps'* to save the figure to a file.

        :param filename:  (Default value = None)
        :type filename:
        :param show:  (Default value = None)
        :type show:
        :returns:
        :rtype:
        """

        import matplotlib as mpl
        mpl.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
        #import pylab as plt

        if self.v0 is None:
            sys.exit('plot(): self.v0 is None!')
            # self.fit()

        if filename is None and show is None:
            show = True

        fig = plt.figure(figsize=(10, 7))
        ax1 = fig.add_subplot(111)
        fig.subplots_adjust(left=0.12, right=0.9, top=0.9, bottom=0.15)
        if self.units == 'bohr':
            if self.eos_string == 'morse':
                plt.plot(self.v, self.e, 'o')
                x = np.linspace(min(self.v), max(self.v), 100)
                y = eval('self.{0}'.format(self.eos_string))(
                    x,
                    self.eos_parameters[0] +
                    self.eMin,
                    self.eos_parameters[1],
                    self.eos_parameters[2],
                    self.eos_parameters[3],)
            elif self.eos_string == 'sjeos':
                plt.plot(self.angstrom2bohr(self.v), self.e / self.ry2ev, 'o')
                x = np.linspace(min(self.v), max(self.v), 100)
                y = self.sjeosfit0(x**-(1.0 / 3.0))
                x = self.angstrom2bohr(x)
                y /= self.ry2ev
            else:
                plt.plot(self.angstrom2bohr(self.v), self.e / self.ry2ev, 'o')
                x = np.linspace(min(self.v), max(self.v), 100)
                y = eval('self.{0}'.format(self.eos_string))(
                    x,
                    self.eos_parameters[0] +
                    self.eMin,
                    self.eos_parameters[1],
                    self.eos_parameters[2],
                    self.eos_parameters[3],)
                x = self.angstrom2bohr(x)
                y /= self.ry2ev

        if self.units == 'angstrom':
            if self.eos_string == 'morse':
                plt.plot(self.bohr2angstrom(self.v), self.e, 'o')
                x = np.linspace(min(self.v), max(self.v), 100)
                y = eval('self.{0}'.format(self.eos_string))(
                    x,
                    self.eos_parameters[0] +
                    self.eMin,
                    self.eos_parameters[1],
                    self.eos_parameters[2],
                    self.eos_parameters[3],)
                x = self.bohr2angstrom(x)
            else:
                plt.plot(self.v, self.e / self.ry2ev, 'o')
                x = np.linspace(min(self.v), max(self.v), 100)
                y = eval('self.{0}'.format(self.eos_string))(
                    x,
                    self.eos_parameters[0] +
                    self.eMin,
                    self.eos_parameters[1],
                    self.eos_parameters[2],
                    self.eos_parameters[3],)
                y /= self.ry2ev

        plt.plot(x, y, '-r')
        plt.ylabel('energy [Ry]')
        ax1.yaxis.set_major_formatter(FormatStrFormatter('%0.6f'))
        ax = plt.gca()
        plt.text(0.2, 0.8, filename, transform=ax.transAxes)
        plt.text(
            0.2,
            0.75,
            'Chi^2 = {0}'.format(
                self.chisqr),
            transform=ax.transAxes)
        plt.text(
            0.2,
            0.7,
            'Red. Chi^2 = {0}'.format(
                self.redchi),
            transform=ax.transAxes)
        plt.text(
            0.2,
            0.65,
            'R^2 = {0}'.format(
                self.rsquared),
            transform=ax.transAxes)
        if self.units == 'bohr':
            plt.xlabel('volume [sws (in Bohr)]')
            plt.title(
                '{0}: E0: {1:.3f} Ry, w0: {2:.3f} Bohr, B0: {3:3.1f} GPa'.format(
                    self.eos_string,
                    self.e0,
                    self.v0,
                    self.B))
        if self.units == 'angstrom':
            plt.xlabel('volume [Angstrom^3]')
            plt.title(
                '{0}: E0: {1:.3f} Ry, V0: {2:.3f} Angstrom^3, B0: {3:3.1f} GPa'.format(
                    self.eos_string,
                    self.e0,
                    self.v0,
                    self.B))

        if show:
            # plt.tight_layout()
            plt.show()
        if filename is not None:
            fig.savefig("fit/{0}.png".format(filename))

        # return f

    def fit_eval(self, sws):
        """Evaluate the fitting function at given points

        :param sws:
        :type sws:
        :returns:
        :rtype:
        """

        y = eval('self.{0}'.format(self.eos_string))(
            sws,
            self.eos_parameters[0] +
            self.eMin,
            self.eos_parameters[1],
            self.eos_parameters[2],
            self.eos_parameters[3],)
        return y

    def ca_fit(self, x, y, n):
        """Fits a polynomial to x vs. y data and calculates
        xmin and ymin from the curve.

        :param x:
        :type x:
        :param y:
        :type y:
        :param n:
        :type n:
        :returns:
        :rtype:
        """

        z = np.polyfit(x, y, n)
        fit = np.poly1d(z)
        xmin = -z[1] / (2.0 * z[0])
        ymin = fit(xmin)
        return xmin, ymin

    def distortion_poly1(self, x, a2, a1):
        """One variable form of the 2nd order polynomial
           which is used to fitting distortion vs. energy
           data in order to find elastic constants.

        :param x:
        :type x:
        :param a2:
        :type a2:
        :param a1:
        :type a1:
        :returns:
        :rtype:
        """

        E = a2 * x**2
        return E

    def distortion_poly2(self, x, a2, a0):
        """Two variable form of the 2nd order polynomial
           which is used to fitting distortion vs. energy
           data in order to find elastic constants.

        :param x:
        :type x:
        :param a2:
        :type a2:
        :param a0:
        :type a0:
        :returns:
        :rtype:
        """

        E = a2 * x**2 + a0
        return E

    def distortion_poly3(self, x, a2, a1, a0):
        """Three variable form of the 2nd order polynomial
           which is used to fitting distortion vs. energy
           data in order to find elastic constants.

        :param x:
        :type x:
        :param a2:
        :type a2:
        :param a1:
        :type a1:
        :param a0:
        :type a0:
        :returns:
        :rtype:
        """

        E = a2 * x**2 + a1 * x + a0
        return E

    def distortion_fit(self,x,y,num=1,title='',ascii_art=True):
        """Fits the distortion_poly function to
        the distortion data.

        num : number of variables in the fitting function
          1 : E=a2*x**2
          2 : E=a2*x**2 + a0
          3 : E=a2*x**2 + a1*x + a0

        The fit coefficient(s) and r-squared describing the accuracy of the fit are returned.

        :param x:
        :type x:
        :param y:
        :type y:
        :param num:  (Default value = 2)
        :type num:
		:param title:  (Default value = '')
		:type title: str
		:param ascii_art: True, if ascii figure of the fit should be drawn.
		:type ascii_art: Boolean
        :returns:
        :rtype:
        """

        if num == 1:
            ydata = y-y[0]
        else:
            ydata = y
        
        # Starting estimates for the actual fit
        z = np.polyfit(x,ydata,2)

		"""
        # TEST: simulate symmetric high-quality 6th order polynomial fit
        eps = 1.0E-6
        zero_included = False
        for i in range(len(x)):
            if x[i]<eps:
                zero_included = True
        if zero_included == True:
            xsim = np.zeros(2*len(x)-1)
            ysim = np.zeros(2*len(x)-1)
        else:
            xsim = np.zeros(2*len(x))
            ysim = np.zeros(2*len(x))
        if zero_included == True:
            xsim[:len(x)-1] = -x[:0:-1]
            xsim[len(x)-1:] = x[:]
            ysim[:len(x)-1] = y[:0:-1]
            ysim[len(x)-1:] = y[:]
        else:
            xsim[:len(x)] = -x[::-1]
            xsim[len(x):] = x[:]
            ysim[:len(x)] = y[::-1]
            ysim[len(x):] = y[:]
                                        
        print(xsim)
        print(ysim)

        print(np.polyfit(xsim,ysim,2)[-3])
        print(np.polyfit(xsim,ysim,3)[-3])
        print(np.polyfit(xsim,ysim,4)[-3])
        print(np.polyfit(xsim,ysim,5)[-3])
        print(np.polyfit(xsim,ysim,6)[-3])
        """

        if num == 1:
            p0 = [z[0]]
        elif num == 2:
            p0 = [z[0],z[2]]
        elif num == 3:
            p0 = [z[0],z[1],z[2]]

        popt,pcov,infodict,mesg,ier=curve_fit(eval('self.distortion_poly{0}'.format(num)),x,ydata,p0)

        # Compute error estimate for the fit
        residuals = infodict['fvec']
        chisqr = np.sum(residuals**2)
        redchi = chisqr/(len(ydata)-4)
        ss_tot = np.sum((ydata-np.mean(ydata))**2)
        rsquared = 1.0-chisqr/ss_tot
        #fitE = self.predicted() + self.eMin

        # Make sure the return variable is always an array
        if type(popt) is type(np.sqrt(1.0)):
            popt = np.array([popt])

        print(popt)

        # Draw the curve
        if ascii_art:
            import time
            
            if num == 1:
                fitti = eval('self.distortion_poly{0}'.format(num))(x,popt[0])
            elif num == 2:
                fitti = eval('self.distortion_poly{0}'.format(num))(x,popt[0],popt[1])
            elif num == 3:
                fitti = eval('self.distortion_poly{0}'.format(num))(x,popt[0],popt[1],popt[2])
            self.ascii_plot(x,ydata,fitti,title)
			
            """
            # Draw simulated high-quality fit:
            fitti = np.poly1d(np.polyfit(xsim,ysim,6))(xsim)
            self.ascii_plot(xsim,ysim,fitti,title)
            """
            # Allow some time to pass so that the gnuplot process finishes before
            # python continues. Otherwise plots might get drawn into the same canvas.
            time.sleep(0.1)
            
        # TEST
        # TEST
        # TEST: use symmetric fit as output
        #return [np.polyfit(xsim,ysim,6)[-3]],rsquared
        # TEST
        # TEST

        return popt,rsquared

    def ascii_plot(self,x,y,z,title=''):
        from pyemto.utilities.utils import run_bash
        import numpy as np
        import subprocess

        gnuplot = subprocess.Popen(["/usr/bin/gnuplot"],stdin=subprocess.PIPE)

        gnuplot.stdin.write("set term dumb 70 20\n")

        # Adjust the legend location
        gnuplot.stdin.write("set key left top\n")

        #gnuplot.stdin.write("plot '-' using 1:2 title '{0}' with lines\n".format(title))
        #gnuplot.stdin.write("plot '-' using 1:2 title 'fit' with points\n")
        gnuplot.stdin.write("plot '-' using 1:2 title '{0}' with points, '-' using 1:2 title 'fit' with lines\n".format(title))

        for i,j in zip(x,y):
            gnuplot.stdin.write("%f %f\n" % (i,j))
        gnuplot.stdin.write("e\n")

        for i,j in zip(x,z):
            gnuplot.stdin.write("%f %f\n" % (i,j))
        gnuplot.stdin.write("e\n")

        #for i,j,k in zip(x,y,z):
        #    gnuplot.stdin.write("%f %f %f\n" % (i,j,k))
        #gnuplot.stdin.write("e\n")

        gnuplot.stdin.write("quit")
        gnuplot.stdin.flush()

        return
