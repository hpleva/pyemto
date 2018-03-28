"""
Python port of the Asturfit program written in Octave:

@article{Otero-De-La-Roza2011c,
author = {Otero-de-la-Roza, A. and Lua{\~{n}}a, V{\'{i}}ctor},
doi = {10.1016/j.cpc.2011.04.016},
journal = {Computer Physics Communications},
number = {8},
pages = {1708--1720},
publisher = {Elsevier B.V.},
title = {{Gibbs2: A new version of the quasi-harmonic model code. I. Robust treatment of the static data}},
url = {http://dx.doi.org/10.1016/j.cpc.2011.04.016 http://linkinghub.elsevier.com/retrieve/pii/S0010465511001470},
volume = {182},
year = {2011}
}

"""

from __future__ import division
import numpy as np
import numpy.polynomial.polynomial as poly

# global constants
hybohr3togpa = 2*14710.50498740275538944426

class polynomial:
    """

    """
    def __init__(self):
        self.c = None
        self.SSerr = None
        self.order = None
        self.data = None
        self.smin = None
        self.w = None
    def __repr__(self):
        return "c    = {0}\nSSerr = {1}\norder = {2}\ndata  = {3}\nw     = {4}\nsmin = {5}\n"\
                .format(self.c,self.SSerr,self.order,self.data,self.w,self.smin)

class minimum:
    """
    Class to store minimums
    """
    def __init__(self):
        self.Vmin = None
        self.fmin = None
        self.Emin = None
        self.err  = 0
    def __repr__(self):
        return "Vmin = {0}\nfmin = {1}\nEmin = {2}\nerr  = {3}\n"\
                .format(self.Vmin,self.fmin,self.Emin,self.err)

class straineval_class:
    """
    Class to store strains
    """
    def __init__(self):
        self.E   = None
        self.E1v = None
        self.E2v = None
        self.E3v = None
        self.E4v = None
        self.E5v = None
        self.p   = None
        self.B   = None
        self.B1p = None
        self.B2p = None
        self.B3p = None
    def __repr__(self):
        string = "E   = {0}\nE1v = {1}\nE2v = {2}\nE3v = {3}\nE4v = {4}\nE5v = {5}\n" +\
                 "p   = {6}\nB   = {7}\nB1p = {8}\nB1p = {9}\nB3p = {10}\n"
        return string.format(self.E,self.E1v,self.E2v,self.E3v,self.E4v,self.E5v,self.p,self.B,self.B1p,self.B2p,self.B3p)

class savg_class:
    def __init__(self):
        self.eqmean = None
        self.eqstd = None
        self.R2 = None
        self.Efit = None
        self.Estd = None
        self.pmean = None
        self.pstd = None
        self.Bmean = None
        self.Bstd = None
        self.B1pmean = None
        self.B1pstd = None
        self.B2pmean = None
        self.B2pstd = None
        self.B3pmean = None
        self.B3pstd = None

class sb_class:
    def __init__(self):
        self.mean = None
        self.std = None
        self.outliers = None
        self.Efit = None
        self.R2 = None

def volume2strain(V, V0, strain='eulerian'):
    """"""
    # Determine the strain:
    if strain == 'eulerian':
        f = ((V/V0)**(-2./3)-1)/2
        return f
    else:
        import sys
        sys.exit('volume2strain: Requested strain form \'{0}\' is unknown!'.format(strain))

def strain2volume(f, V0, strain='eulerian'):
    """"""
    # Determine the strain:
    if strain == 'eulerian':
        V = (f*2+1)**(-3./2) * V0
        return V
    else:
        import sys
        sys.exit('strain2volume: Requested strain form \'{0}\' is unknown!'.format(strain))

def check_noise(vols,ens,show_plot=False):
    """Detect outliers/phase transitions in E vs. V data."""
    status = 0
    vref = np.median(vols)
    rk = [1,3,4,5]
    strain = 'eulerian'
    cf,sf = avgstrainfit(vols,ens,vref,nmax=10,strain=strain,LOG=1,nargout=2)

    # Is this a noisy dataset?
    fiterr = np.abs(sf.eqstd/sf.eqmean)
    worse = np.max(fiterr)

    # TEST:
    #worse = 0.5

    print('################################################')
    print('Noise level: ~0.0 => Smooth, >1.0 => Very noisy.')
    print('################################################')
    if worse < 0.01:
        status = 0
        print("Datafile appears very smooth: Noise level = {0}\nCheck the E(V) plot, anyway.\n".format(worse))
        return vols, ens
    elif worse > 1.0:
        status += 20
        print("Noise is a serious problem in your data: Noise level = {0}\n".format(worse))
        sample = 10000
    else:
        status += 10
        print("Some level of noise exists: Noise level = {0}\nCheck carefully the next results.\n".format(worse))
        sample = 1000

    #"""
    # Try first a bootstrap fit to detect possible outliers:
    print("\nA bootstrap fit will be tried to detect possible outliers:\n")
    print("Bootstrap sample size: {0}\n".format(sample))
    cb,sb = strainbootstrap(vols, ens, vref, ndeg=12, nsample=sample, strain=strain, LOG=1, nargout=2)
    bad = sb.outliers
    good = [i for i in range(len(vols))]
    for i in range(len(bad)):
        good.remove(bad[i])
    #fileplot = sprintf("%s-evn.eps", rootname)
    #plot(v,sb.Efit,'-r;BS fit;', v(good),e(good),'ob;good data;' \
    #    , v(bad),e(bad),'om;bad data;', 'markersize', 12)
    #xlabel('V (bohr^3)'); ylabel('E (hartree)'); grid('on')
    #print(fileplot, '-FHelvetica:24', '-depsc')
    #printf('See the E(V) curve in the file: %s\n', fileplot)
    if (len(bad) > 0):
        status += 1
    #"""

    #"""
    # Remove the detected outliers:
    vg = vols[good]; eg = ens[good]
    print("\nRemove the outliers and keep {0} points.\n".format(len(vg)))
    #"""
    #
    if show_plot == True:
        import matplotlib.pyplot as plt
        from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
        fig = plt.figure(figsize=(10, 7))
        ax1 = fig.add_subplot(111)
        fig.subplots_adjust(left=0.12, right=0.9, top=0.9, bottom=0.15)
        #
        plt.plot(vols,sb.Efit,color='black')
        plt.plot(vg, eg, 'o', color='green')
        #plt.plot(vols[bad], ens[bad], 'o', mfc='none', markersize=12, color='red', linewidth='3')
        plt.scatter(vols[bad], ens[bad], s=80, facecolors='none', edgecolors='r', linewidth='3')
        plt.tight_layout()
        plt.show()
        print(sb.Efit)

    return vg,eg

def avgstrainfit(V, E, V0, nmax=16, MODE=1, strain='eulerian', LOG=0, nargout=1):
    """avgstrainfit - Fit to an average strain polynomials."""

    import sys
    if (len(V) != len(E)):
        sys.exit('avgstrainfit: V and E must be vectors of the same length!')
    elif (len(V) < 7):
        sys.exit('avgstrainfit: dataset must have at least 7 points!')
    elif MODE != 1 and MODE != 2:
        sys.exit('avgstrainfit: weighting mode must be 1 or 2!')

    ndata = len(V)
    Vrange = [np.amin(V), np.amax(V)]

    # Determine the maximum degree of the polynomials:
    if (nmax < 0):
        MaxDegree = int(np.amin([ndata-5, np.trunc(ndata/2.)]))
    else:
        MaxDegree = int(np.amin([nmax,np.amin([ndata-5, np.trunc(ndata/2.)])]))

    # Some statistics of the averaging proccess:
    Morder = 0
    npol = 0
    pol = []
    s2 = []
    for n in range(2,MaxDegree+1):
        c = strainfit(V, E, V0, n, strain, 0, nargout=1)
        sm = strainmin(c, V0, Vrange, strain)
        Efit = strainevalE(c, V0, V, strain)
        SSerr = np.sum((E-Efit)**2)
        pol.append(polynomial())
        pol[npol].c = c
        pol[npol].SSerr = SSerr
        s2.append(SSerr)
        pol[npol].order = n
        pol[npol].data = ndata
        pol[npol].smin = sm
        Morder = np.max([n,Morder])
        npol += 1
        #print(c)
        #print('n,Morder,npol = ',n,Morder,npol)
        #print(sm)
        #print(Efit)
        #print(SSerr)

    #print('s2 = ',s2)
    #print(pol)

    # Get the polynomial weights:
    SSmin = np.amin(s2)
    Q = 0
    if (MODE == 1):
        for k in range(npol):
            w = (pol[k].SSerr/SSmin) * (pol[k].order / pol[k].data)
            pol[k].w = np.exp(-w)
            Q += pol[k].w
    elif (MODE == 2):
        for k in range(npol):
            w = (pol[k].SSerr/SSmin) * (pol[k].order / pol[k].data)
            pol[k].w = np.exp(-w*w)
            Q += pol[k].w
    ##ww = zeros(1,1:npol);
    ww = np.zeros(npol)
    for k in range(npol):
        ww[k] = pol[k].w / Q
        pol[k].w = pol[k].w / Q

    #print(pol)

    # Form the average polynomial:
    cavg = np.zeros(Morder+1)
    for k in range(npol):
        n = pol[k].order
        cavg[Morder-n:] += pol[k].c * pol[k].w

    #"""
    # Extra output and optional LOG record:
    # If the average of polynomials has a minimum, analyze the equilibrium
    # geometry.
    # Otherwise, analyze the reference volume.
    if (LOG > 0 or nargout > 1):
        avgmin = strainmin(cavg, V0, Vrange, strain)
        if (avgmin.err == 0):
            # Analyze and report the equilibrium geometry:
            if (LOG > 0):
                print('\n\navgpolyfit: AVERAGE OF {0} STRAIN POLYNOMIALS'.format(strain))
                print('Volume reference (V0): {0:.6f}'.format(V0))
                print('Range of degrees: 2--{0}'.format(MaxDegree))
                print('Number of polynomials: {0}'.format(npol))
                print('\nProperties at the minimum of each polynomial:')
                print('--i-- npol data order --SSerr-- ---w---- ----Vmin--- -----Emin---  ----Bmin--- ---B1min--- ---B2min---   ---B3min---')
            isrt = np.argsort(ww)[::-1]
            srt = ww[isrt]
            pmean = np.zeros(6); pstd = np.zeros(6)
            for k in range(npol):
                k1 = isrt[k]
                ###[smin] = strainmin(pol{k1}.c, V0, Vrange, strain);
                smin = pol[k1].smin
                prop = straineval(pol[k1].c, V0, smin.Vmin, strain)
                if (LOG > 0 and k <= 25):
                    # In the conversion of the bulk modulus and derivatives we
                    # assume the units: volume (bohr^3), energy (Hy).
                    output_string  = "{0:4} {1:4} {2:4} {3:4}".format(k,k1+1,pol[k1].data,pol[k1].order)
                    output_string += "   {0:9.2e} {1:8.6f}".format(pol[k1].SSerr, pol[k1].w)
                    output_string += " {0:10.6f}".format(smin.Vmin)
                    output_string += " {0:13.6f}".format(prop.E)
                    output_string += "  {0:10.6f}".format(prop.B * hybohr3togpa)
                    output_string += "  {0:9.6f}".format(prop.B1p)
                    output_string += "  {0:13.9f}".format(prop.B2p / hybohr3togpa)
                    output_string += " {0:13.9f}".format(prop.B3p / hybohr3togpa**2)
                    print(output_string)
                pr = np.array([smin.Vmin, prop.E, prop.B, prop.B1p, prop.B2p, prop.B3p])
                pmean = pmean + pr * pol[k1].w
                pstd = pstd + (pr**2) * pol[k1].w
            pstd = np.lib.scimath.sqrt(pstd - pmean**2)
            if (LOG > 0):
                print('\nAverage properties (weighted polynomials):')
                print('------ ---volume-- ---energy--   --B-(GPa)-- ----B1p---- B2p-(1/GPa) B3p--(1/GPa^2)')
                print('-mean- {0:11.6f} {1:11.6f} {2:11.6f} {3:11.6f} {4:11.6f} {5:14.9f}'.format(pmean[0]
                      , pmean[1], pmean[2]*hybohr3togpa, pmean[3], pmean[4]/hybohr3togpa
                      , pmean[5]/hybohr3togpa**2))
                print('stdvev {0:11.6f} {1:11.6f} {2:11.6f} {3:11.6f} {4:11.6f} {5:14.9f}\n'.format(float(np.real(pstd[0]))
                      , float(np.real(pstd[1])), float(np.real(pstd[2]))*hybohr3togpa, float(np.real(pstd[3]))
                      , float(np.real(pstd[4]))/hybohr3togpa, float(np.real(pstd[5]))/hybohr3togpa**2))
        else:
            # Analyze and report the reference geometry:
            if (LOG > 0):
                print('\n\navgpolyfit: AVERAGE OF {0} STRAIN POLYNOMIALS'.format(strain))
                print('Volume reference (V0): {0:.6f}'.format(V0))
                print('Range of degrees: 2--{0}'.format(MaxDegree))
                print('Number of polynomials: {0}'.format(npol))
                print('\nProperties at the reference volume: {0}'.format(V0))
                print('--i-- npol data order --SSerr-- ---w---- -----Eref---  ----Bref--- ---B1ref--- ---B2ref---   ---B3ref---')
            isrt = np.argsort(ww)[::-1]
            srt = ww[isrt]
            pmean = np.zeros(6); pstd = np.zeros(6)
            for k in range(npol):
                k1 = isrt[k]
                prop = straineval(pol[k1].c, V0, V0, strain)
                if (LOG > 0 and k <= 25):
                    # In the conversion of the bulk modulus and derivatives we
                    # assume the units: volume (bohr^3), energy (Hy).
                    output_string  = "{0:4} {1:4} {2:4} {3:4}".format(k,k1+1,pol[k1].data,pol[k1].order)
                    output_string += "   {0:9.2e} {1:8.6f}".format(pol[k1].SSerr, pol[k1].w)
                    output_string += " {0:13.6f}".format(prop.E)
                    output_string += "  {0:10.6f}".format(prop.B * hybohr3togpa)
                    output_string += "  {0:9.6f}".format(prop.B1p)
                    output_string += "  {0:13.9f}".format(prop.B2p / hybohr3togpa)
                    output_string += " {0:13.9f}".format(prop.B3p / hybohr3togpa**2)
                    print(output_string)
                pr = np.array([V0, prop.E, prop.B, prop.B1p, prop.B2p, prop.B3p])
                pmean = pmean + pr * pol[k1].w
                pstd = pstd + (pr**2) * pol[k1].w
            pstd = np.lib.scimath.sqrt(pstd - pmean**2)
            if (LOG > 0):
                print('\nAverage properties at the ref. volume: {0}'.format(V0))
                print('------ ---energy--   --B-(GPa)-- ----B1p---- B2p-(1/GPa) B3p--(1/GPa^2)')
                print('-mean- {0:11.6f} {1:11.6f} {2:11.6f} {3:11.6f} {4:14.9f}'.format(
                        pmean[1], pmean[2]*hybohr3togpa, pmean[3], pmean[4]/hybohr3togpa
                      , pmean[5]/hybohr3togpa**2))
                print('stdvev {0:11.6f} {1:11.6f} {2:11.6f} {3:11.6f} {4:14.9f}\n'.format(
                        float(np.real(pstd[1])), float(np.real(pstd[2]))*hybohr3togpa, float(np.real(pstd[3]))
                      , float(np.real(pstd[4]))/hybohr3togpa, float(np.real(pstd[5]))/hybohr3togpa**2))

        Efit = strainevalE(cavg, V0, V, strain)
        SSerr = np.sum((E-Efit)**2)
        SStot = np.sum((E-np.mean(E))**2)
        R2 = 1. - SSerr / SStot
        if (nargout > 1):
            savg = savg_class()
            savg.eqmean = pmean
            savg.eqstd = pstd
            savg.R2 = R2
            savg.Efit = Efit
    #"""
    if (nargout > 1):
        # Putting error bars to the pressure, bulk modulus, etc at all volumes
        savg.Emean, savg.Estd = np.zeros(len(V)), np.zeros(len(V))
        savg.pmean, savg.pstd = np.zeros(len(V)), np.zeros(len(V))
        savg.Bmean, savg.Bstd = np.zeros(len(V)), np.zeros(len(V))
        savg.B1pmean, savg.B1pstd = np.zeros(len(V)), np.zeros(len(V))
        savg.B2pmean, savg.B2pstd = np.zeros(len(V)), np.zeros(len(V))
        savg.B3pmean, savg.B3pstd = np.zeros(len(V)), np.zeros(len(V))
        for k in range(npol):
            k1 = isrt[k]
            prop = straineval(pol[k1].c, V0, V, strain)
            savg.Emean = savg.Emean + prop.E * pol[k1].w
            savg.Estd  = savg.Estd  + (prop.E**2) * pol[k1].w
            savg.pmean = savg.pmean + prop.p * pol[k1].w
            savg.pstd  = savg.pstd  + (prop.p**2) * pol[k1].w
            savg.Bmean = savg.Bmean + prop.B * pol[k1].w
            savg.Bstd  = savg.Bstd  + (prop.B**2) * pol[k1].w
            savg.B1pmean = savg.B1pmean + prop.B1p * pol[k1].w
            savg.B1pstd  = savg.B1pstd  + (prop.B1p**2) * pol[k1].w
            savg.B2pmean = savg.B2pmean + prop.B2p * pol[k1].w
            savg.B2pstd  = savg.B2pstd  + (prop.B2p**2) * pol[k1].w
            savg.B3pmean = savg.B3pmean + prop.B3p * pol[k1].w
            savg.B3pstd  = savg.B3pstd  + (prop.B3p**2) * pol[k1].w
        savg.Estd = np.sqrt(savg.Estd - savg.Emean**2)
        savg.pstd = np.sqrt(savg.pstd - savg.pmean**2)
        savg.Bstd = np.sqrt(savg.Bstd - savg.Bmean**2)
        savg.B1pstd = np.sqrt(savg.B1pstd - savg.B1pmean**2)
        savg.B2pstd = np.sqrt(savg.B2pstd - savg.B2pmean**2)
        savg.B3pstd = np.sqrt(savg.B3pstd - savg.B3pmean**2)

    if nargout == 1:
        return cavg
    elif nargout == 2:
        return cavg,savg
    #"""

def strainfit(V, E, V0, ndeg=4, strain='eulerian', LOG=0, nargout=1):
    """"""
    # Determine the strain:
    if strain == 'eulerian':
        f = ((V/V0)**(-2./3)-1)/2
        f0 = 0
        if (LOG > 0 or nargout > 1):
            fpv = -(f + f + 1)**(5/2) / (3*V0)
            fppv = (f + f + 1)**4 * (5/(3*V0)**2)
    else:
        import sys
        sys.exit('strainfit: Requested strain form \'{0}\' is unknown!'.format(strain))

    #c = np.polyfit(f, E, ndeg)
    c = poly.polyfit(f, E, ndeg)[::-1]


    #TODO
    """
   if (LOG > 0 | nargout > 1)
      # Evaluate the fitting polynomial and get the determination coefficient:
      Efit = polyval(c, f);
      SSerr = sum((E-Efit).^2);
      SStot = sum((E-mean(E)).^2);
      R2 = 1 - SSerr / SStot;
      # Find the polynomial minimum and its properties:
      # 1. evaluate the polynomial in a fine grid and get the best grid point
      # 2. get the roots of the polynomial derivative
      # 3. select the real roots that give a positive second derivative
      # 4. get the candidate from (3) that is closest to the best grid point
      ff = linspace(min(f), max(f), 51);
      ee = polyval(c, ff);
      [emin,imin] = min(ee);
      ffmin = ff(imin);
      c1 = polyder(c);
      c2 = polyder(c1);
      rr = roots(c1);
      ipos = find(abs(imag(rr)) <= 1e-15 & polyval(c2,rr) > 0);
      if (length(ipos) < 1)
         printf('strainfit algorithmic error (1):\n');
         printf('strain & ndeg: %s %d\n', strain, ndeg);
         printf('coefs. of polynomial and derivative:\n');
         printf('-icoef----real(c)---imag(c)---real(c1)---imag(c1)---\n');
         for i = 1:length(c)-1
            printf('%6d %22.15e %22.15e %22.15e %22.15e\n', i, real(c(i)), imag(c(i)), real(c1(i)), imag(c1(i)));
         endfor
         printf('%6d %22.15e %22.15e\n', i, real(c(end)), imag(c(end)));
         printf('-iroot---real---imag---deriv2---\n');
         for i = 1:length(rr)
             printf('%6d %13.6e %13.6e %13.6e\n', i, real(rr(i)), imag(rr(i)), polyval(c2,rr(i)));
         endfor
         ###error('strainfit: fitted polynomial has no minima!');
         s.err = 1;
         return
      elseif (length(ipos) == 1)
         fmin = rr(ipos);
      else
         rr = rr(ipos);
         rx = abs(rr .- ffmin).**2;
         [rxmin,irxmin] = min(rx);
         fmin = rr(irxmin);
      endif
      Emin = polyval(c, fmin);
      E1min = polyval(c1, fmin);
      c2 = polyder(c1); E2min = polyval(c2, fmin);
      c3 = polyder(c2); E3min = polyval(c3, fmin);
      if (strcmp(strain, 'eulerian'))
         Vmin = V0 * (1+2*fmin)^(-3/2);
         fpvmin = -(2*fmin+1)^(5/2)/(3*V0);
         fppvmin = (2*fmin+1)^4 * (5/(3*V0)**2);
         fpppvmin = -(2*fmin+1)^(11/2) * (40/(3*V0)**3);
      Bmin = Vmin * (E2min * fpvmin^2 + E1min * fppvmin);
      Bfmin = E3min*Vmin*fpvmin^2 + E2min*(fpvmin+3*fppvmin*Vmin) \
            + E1min*(fppvmin+Vmin*fpppvmin)/fpvmin;
      pfmin = -(E2min*fpvmin-E1min*fppvmin/fpvmin);
      Bpmin = Bfmin/pfmin;
   endif

   if (LOG > 0)
      printf('\n\n*************\n')
      printf('* strainfit * Polynomial fitting to the %s strain\n', strain)
      printf('*************\n')
      printf('Reference volume (V0):  %.6f\n', V0);
      printf('Polynomial degree :     %d\n', ndeg);
      printf('Determin. coef. (R2):   %.12f\n', R2);
      printf('Vol. at min. (Vmin):    %.6f\n', Vmin);
      printf('Energy at min. (Emin):  %.9g\n', Emin);
      printf('B at min.      (Bmin):  %.9g\n', Bmin);
      printf('Bp at min.    (Bpmin):  %.9g\n', Bpmin);
      printf('All values in the input units\n');
   endif

   if (nargout > 1)
      s.V0 = Vmin;
      s.E0 = Emin;
      s.B0 = [Bmin, Bpmin];
      s.R2 = R2;
      s.S2 = SSerr;
      s.Efit = Efit;
      s.f = f;
      s.err = 0;
   endif

    """
    return c

def strainmin(c, V0, Vrange, strain='eulerian'):
    """"""
    LOG = 0
    frange = volume2strain(Vrange, V0, strain)
    #
    # Get the roots of the polynomial derivative:
    # Choose only those that are real and produce a positive second derivative.
    #
    grid_done = False
    s = minimum()
    c1 = np.polyder(c)
    c2 = np.polyder(c1)
    rr = np.roots(c1)
    ipos = []
    cond1 = np.abs(rr.imag) <= 1e-15
    cond2 = np.polyval(c2,rr) > 0
    for i in range(len(rr)):
        if cond1[i] == True and cond2[i] == True:
            ipos.append(i)
    #
    if len(ipos) < 1:
        #TODO
        """
        #
        # There are no left roots. An error in the roots() routine?
        #
        if (LOG)
         printf('polymin algorithmic error (1):\n');
         printf('coefs. of polynomial and derivative:\n');
         printf('-i----real(c)---imag(c)---real(c1)---imag(c1)---\n');
         for i = 1:length(c)-1
            printf('%3d %12.5e %12.5e %12.5e %12.5e\n', \
                   i, real(c(i)), imag(c(i)), real(c1(i)), imag(c1(i)));
         endfor
         printf('%3d %12.5e %12.5e\n', length(c), real(c(end)), imag(c(end)));
         printf('-iroot---real---imag---deriv2---\n');
         for i = 1:length(rr)
            printf('%6d %13.6e %13.6e %13.6e\n', \
                   i, real(rr(i)), imag(rr(i)), polyval(c2,rr(i)));
         endfor
        endif
        """
        s.err = 1
    elif len(ipos) == 1:
        #
        # There is a single root. Check if it is inside range.
        #
        s.fmin = np.real(rr[ipos[0]])
        if s.fmin < np.amin(frange)*0.98:
            s.err = 1
        elif s.fmin > np.amax(frange)*1.02:
            s.err = 1
        else:
            s.err = 0
    else:
        #
        # There are several candidates. Evaluate the polynomial in a grid
        # of points within the range. Get the position of the minimum in the
        # grid and choost the root closest to that point.
        # Check again if it is inside range.
        #
        ff = np.linspace(np.amin(frange), np.amax(frange), 51)
        ee = np.polyval(c,ff)
        eemin = np.amin(ee); ieemin = np.argmin(ee)
        ffmin = ff[ieemin]
        grid_done = True
        rr = np.real(rr[ipos])
        rx = np.abs(rr - ffmin)**2
        rxmin = np.amin(rx); irxmin = np.argmin(rx)
        s.fmin = rr[irxmin]
        if s.fmin < np.amin(frange)*0.98:
            s.err = 1
        elif s.fmin > np.amax(frange)*1.02:
            s.err = 1
        else:
            s.err = 0

    if (s.err != 0):
        print('Running Newton-Raphson')
        #
        # The above method has failed. Try Newton-Raphson.
        #
        # Evaluate the grid points if that has not been done before:
        if grid_done == False:
            ff = np.linspace(np.amin(frange), np.amax(frange), 51)
            ee = np.polyval(c,ff)
            eemin = np.amin(ee); ieemin = np.argmin(ee)
            ffmin = ff[ieemin]
        CONV = 1e-8; MAXIT = 30
        fnew = ffmin; delta = 1e30; it = 0
        for i in range(MAXIT):
            fold = fnew; dold = delta
            fnew = fold - np.polyval(c1,fold) / np.polyval(c2,fold)
            delta = np.abs(fnew-fold)
            it += 1
            if delta < CONV:
                break
        s.fmin = fnew
        # Check
        if delta > CONV:
            s.err = 1
        elif s.fmin < np.amin(frange)*0.98:
            s.err = 1
        elif s.fmin > np.amax(frange)*1.02:
            s.err = 1
        else:
            s.err = 0

    s.Vmin = strain2volume(s.fmin, V0, strain)
    s.Emin = np.polyval(c, s.fmin)
    return s

def strainevalE(c, V0, V, strain='eulerian'):
    """"""
    #if (nargin < 3 | nargin > 4)
    #    print_usage ();
    #endif

    # Determine the strain:
    if strain == 'eulerian':
        f = ((V/V0)**(-2./3)-1)/2
    elif strain == 'natural':
        f = log(V/V0)/3
    elif strain == 'lagrangian':
        f = ((V/V0)**(2./3)-1)/2
    elif strain == 'infinitesimal':
        f = -(V/V0)**(-1./3)+1
    elif strain == 'quotient' or strain == 'x1':
        f = V/V0
    elif strain == 'x3':
        f = (V/V0)**(1./3)
    elif strain == 'xinv3':
        f = (V/V0)**(-1./3)
    elif strain == 'V':
        f = V
    else:
        import sys
        sys.exit('strainevalE: Requested strain form \'{0}\' is unknown!'.format(strain))

    # Evaluate E(f):
    E = np.polyval(c, f)
    return E

def straineval(c, V0, V, strain='eulerian'):
    """"""
    # Determine the strain:
    if strain == 'eulerian':
        f = ((V/V0)**(-2./3)-1)/2
        f2 = f + f + 1
        ss = -f2**(3./2)/(3*V0)
        f1v = f2**(5./2) * (-1./(3*V0))
        f2v = f1v * ss * 5
        f3v = f2v * ss * 8
        f4v = f3v * ss * 11
        f5v = f4v * ss * 14
        f6v = f5v * ss * 17
        ###f2v = f2.^4 * (5/(3*V0)^2)
        ###f3v = f2.^(11/2) * (-40/(3*V0)^3)
        ###f4v = f2.^7 * (440/(3*V0)^4)
        ###f5v = f2.^(17/2) * (-6160/(3*V0)^5)
    else:
        import sys
        sys.exit('straineval: Requested strain form \'{0}\' is unknown!'.format(strain))

    s = straineval_class()
    # Evaluate E(f) and the derivatives of E versus f:
    s.E = np.polyval(c, f)
    #
    c1 = np.polyder(c);  E1f = np.polyval(c1, f)
    c2 = np.polyder(c1); E2f = np.polyval(c2, f)
    c3 = np.polyder(c2); E3f = np.polyval(c3, f)
    c4 = np.polyder(c3); E4f = np.polyval(c4, f)
    c5 = np.polyder(c4); E5f = np.polyval(c5, f)

    # Get the derivatives of E versus V:
    s.E1v = E1f*f1v
    s.E2v = E2f*f1v**2 + E1f*f2v
    s.E3v = E3f*f1v**3 + E2f*f1v*f2v*3 + E1f*f3v
    s.E4v = E4f*f1v**4 + E3f*f1v**2*f2v*6 + \
            E2f*(f1v*f3v*4+f2v**2*3) + E1f*f4v
    s.E5v = E5f*f1v**5 + E4f*f1v**3*f2v*10 + \
            E3f*(f1v**2.*f3v*10+f1v*f2v**2*15) + \
            E2f*(f1v*f4v*5+f2v*f3v*10) + E1f*f5v

    # Get the pressure, the bulk modulus, and its derivatives versus pressure:
    s.p = -s.E1v
    s.B = V * s.E2v
    s.B1p = -(V * s.E3v + s.E2v) / s.E2v
    s.B2p = ((s.E4v*s.E2v - s.E3v**2) * V + s.E3v*s.E2v) / s.E2v**3
    s.B3p = -((s.E2v**2*s.E5v-s.E2v*s.E3v*s.E4v*4+s.E3v**3*3)*V \
            + s.E2v**2*s.E4v*2 - s.E2v*s.E3v**2*3) / s.E2v**5
    return s


def strainbootstrap(V, E, V0, ndeg=12, nsample=100, strain='eulerian', LOG=0, nargout=1):
    """"""
    n = len(V)
    Vrange = [np.amin(V), np.amax(V)]
    nreject = 0
    npol = 0
    #s2 = np.zeros(nsample)
    pol = []
    s2 = []
    Morder = 0
    for k in range(nsample):
        rr = np.random.rand(n)
        r,  = np.where(rr >= 0.5)
        nr = len(r)
        #print(r)
        for nd in range(2,ndeg+1):
            if (nr < nd + 5):
                # Too few points for a meaningful fit
                nreject += 1
            else:
                c = strainfit(V[r], E[r], V0, nd, strain, 0)
                Morder = np.amax([Morder,len(c)-1])
                sm = strainmin(c, V0, Vrange, strain)
                Efit = strainevalE(c, V0, V[r], strain)
                SSerr = np.sum((E[r]-Efit)**2)
                pol.append(polynomial())
                pol[npol].c = c
                #s2[npol] = SSerr
                s2.append(SSerr)
                pol[npol].SSerr = SSerr
                pol[npol].deg  = nd
                pol[npol].data = nr
                pol[npol].V0 = sm.Vmin
                pol[npol].smin = sm
                npol +=1

    # Get the polynomial weights:
    #print('nreject = ',nreject)
    #print('s2 =',s2)
    SSmin = np.amin(s2)
    Q = 0
    for k in range(npol):
        w = (pol[k].SSerr/SSmin) * (pol[k].deg / pol[k].data)
        pol[k].w = np.exp(-w)
        Q += pol[k].w

    ww = np.zeros(npol)
    for k in range(npol):
        ww[k] = pol[k].w / Q
        pol[k].w = pol[k].w / Q

    # Form the average polynomial:
    cb = np.zeros(Morder+1)
    for k in range(npol):
        n = len(pol[k].c) - 1
        cb[Morder-n:] += pol[k].c * pol[k].w

    # Extra output and optional LOG record:
    # If the average of polynomials has a minimum, analize the equilibrium
    # geometry.
    # Otherwise analyze the reference volume.
    if (LOG>0 or nargout > 1):
        avgmin = strainmin(cb, V0, Vrange, strain)
        if (avgmin.err == 0):
            if (LOG>0):
                print('\n\nsb2: BOOTSTRAP ANALYSIS OF {0} STRAIN POLYNOMIALS'.format(strain))
                print('Volume reference (V0): {0:.6f}'.format(V0))
                print('Fitting mode: average of polynomials up to {0} degree'.format(ndeg))
                print('Samples tried, Pol. rejected/accepted:: {0} {1}/{2}'.format(nsample, nreject, npol))
                print('\nProperties at the minimum of each polynomial:')
                print('--i--    npol data order  --SSerr-- ---w---- ----Vmin--- -----Emin---  ----Bmin--- ---B1min--- ---B2min---   ---B3min---')
            isrt = np.argsort(ww)[::-1]
            srt = ww[isrt]
            pmean  = np.zeros(6); pstd  = np.zeros(6)
            pmean2 = np.zeros(6); pstd2 = np.zeros(6)
            for k in range(npol):
                k1 = isrt[k]
                smin = pol[k1].smin
                Vmin = smin.Vmin
                prop = straineval(pol[k1].c, V0, Vmin, strain)
                if (LOG>0 and k <= 25):
                    # In the conversion of the bulk modulus and derivatives we
                    # assume the units: volume (bohr^3), energy (Hy).
                    output_string  = "{0:4} {1:8} {2:4} {3:4}".format(k,k1+1,pol[k1].data,pol[k1].deg)
                    output_string += "   {0:9.2e} {1:8.6f}".format(pol[k1].SSerr, pol[k1].w)
                    output_string += " {0:10.6f}".format(Vmin)
                    output_string += " {0:13.6f}".format(prop.E)
                    output_string += "  {0:10.6f}".format(prop.B * hybohr3togpa)
                    output_string += "  {0:9.6f}".format(prop.B1p)
                    output_string += "  {0:13.9f}".format(prop.B2p / hybohr3togpa)
                    output_string += " {0:13.9f}".format(prop.B3p / hybohr3togpa**2)
                    print(output_string)
                pr = np.array([Vmin, prop.E, prop.B, prop.B1p, prop.B2p, prop.B3p])
                pmean = pmean + pr * pol[k1].w
                pstd = pstd + (pr**2) * pol[k1].w
                pmean2 = pmean2 + pr * (1./npol)
                pstd2 = pstd2 + (pr**2) * (1./npol)
            pstd  = np.lib.scimath.sqrt(pstd - pmean**2)
            pstd2 = np.lib.scimath.sqrt(pstd2 - pmean2**2)
            if (LOG > 0):
                print('\nAverage properties (weighted polynomials):')
                print('------ ---volume-- ---energy--   --B-(GPa)-- ----B1p---- B2p-(1/GPa) B3p--(1/GPa^2)')
                print('-mean- {0:11.6f} {1:11.6f} {2:11.6f} {3:11.6f} {4:11.6f} {5:14.9f}'.format(pmean[0]
                      , pmean[1], pmean[2]*hybohr3togpa, pmean[3], pmean[4]/hybohr3togpa
                      , pmean[5]/hybohr3togpa**2))
                print('stdvev {0:11.6f} {1:11.6f} {2:11.6f} {3:11.6f} {4:11.6f} {5:14.9f}\n'.format(float(np.real(pstd[0]))
                      , float(np.real(pstd[1])), float(np.real(pstd[2]))*hybohr3togpa, float(np.real(pstd[3]))
                      , float(np.real(pstd[4]))/hybohr3togpa, float(np.real(pstd[5]))/hybohr3togpa**2))
                print('\nAverage properties (equally weighted polynomials):')
                print('------ ---volume-- ---energy--   --B-(GPa)-- ----B1p---- B2p-(1/GPa) B3p--(1/GPa^2)')
                print('-mean- {0:11.6f} {1:11.6f} {2:11.6f} {3:11.6f} {4:11.6f} {5:14.9f}'.format(pmean2[0]
                      , pmean2[1], pmean2[2]*hybohr3togpa, pmean2[3], pmean2[4]/hybohr3togpa
                      , pmean2[5]/hybohr3togpa**2))
                print('stdvev {0:11.6f} {1:11.6f} {2:11.6f} {3:11.6f} {4:11.6f} {5:14.9f}\n'.format(float(np.real(pstd2[0]))
                      , float(np.real(pstd2[1])), float(np.real(pstd2[2]))*hybohr3togpa, float(np.real(pstd2[3]))
                      , float(np.real(pstd2[4]))/hybohr3togpa, float(np.real(pstd2[5]))/hybohr3togpa**2))
                print('\nFor a smooth dataset weighted and equally weighted -means- and stdvevs should be very similar!\n')
        #
        Efit = strainevalE(cb, V0, V, strain)
        SSerr = np.sum((E-Efit)**2)
        SStot = np.sum((E-np.mean(E))**2)
        R2 = 1. - SSerr / SStot
        if (nargout > 1):
            sb = sb_class()
            sb.mean = pmean
            sb.std = pstd
            sb.R2 = R2
            sb.Efit = Efit

        if (LOG>0):
            print('\n\nProperties of the average polynomial:')
            #print('Vmin, Emin: {0:.6f} {1:.6f}'.format(Vmin, prop.E))
            print('Vmin, Emin: {0:.6f} {1:.6f}'.format(avgmin.Vmin, avgmin.Emin))
            print('SSerr, SStot, R2, 1-R2: {0:.9e} {1:.9e} {2:.12f} {3:.2e}'
                  .format(SSerr, SStot, R2, 1-R2))
            nd = np.zeros(npol)
            for i in range(npol):
                nd[i] = pol[i].data
            print('Data (min/max/mean/std): {0} {1} {2:.4f} {3:.4f}'
                  .format(np.amin(nd), np.amax(nd), np.mean(nd), np.std(nd)))

        # To detect outliers we check the residuals for each point, i.e.
        # the difference between the actual energy and the value predicted
        # by the average polynomial. A robust measurement of what is normal
        # for a residual is provided by the median (but not the mean, which
        # too sensitive to outliear with large residuals):
        if (nargout > 1 or LOG > 0):
            alambda = 10.0
            Edif = np.abs(Efit - E)
            Ecentral = np.median(Edif)
            Edisp = np.median(np.abs(Edif-Ecentral))
            iout = np.where(Edif > Ecentral+alambda*Edisp)
            sb.outliers, = iout
            #print('sb.outliers = ',sb.outliers)
            if (LOG>0):
                print('\nDetect outlier points (lambda={0:.6f}):'.format(alambda))
                print('Ecentral (median): {0:.6e}'.format(Ecentral))
                print('Edisp    (median): {0:.6e}'.format(Edisp))
                print('-npt- remove ---volume-- ---energy-- ---error---')
                output_string = ''
                for i in range(len(V)):
                    #print('{0:5} '.format(i))
                    #print('{0} '.format(Edif[i] > Ecentral+alambda*Edisp))
                    #print('{0:11.6f} {1:11.6f} '.format(V[i], E[i]))
                    #print('{0:11.4e}'.format(Edif[i]))
                    output_string += '{0:5} '.format(i)
                    output_string += '{0} '.format(Edif[i] > Ecentral+alambda*Edisp)
                    output_string += '{0:11.6f} {1:11.6f} '.format(V[i], E[i])
                    output_string += '{0:11.4e}\n'.format(Edif[i])
                print(output_string)
                output_string = '\nOutliers ({0}):'.format(len(iout))
                for i in range(len(iout)):
                    output_string += ' {0},'.format(iout[i])
                print(output_string)

    if nargout == 1:
        return cb
    elif nargout == 2:
        return cb,sb


def guess_jumps(vols, ens, deltaE=0.0, LOG=0):
    """"""
    n = len(ens)

    # Linear interpolation
    E1 = np.zeros(n)
    pass
