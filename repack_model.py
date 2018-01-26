import numpy as np
from math import sqrt, pi, log10, e, exp, log
from sumdigit import create_sum
from scipy.optimize import broyden1
import pickle
import matplotlib.pyplot as plt
from multiprocessing import Pool, Lock, Manager
from time import time

import matplotlib.pyplot  as plt
from scipy.interpolate import interp1d

LAME = 1.7 #[cm] = mu_e \tau_e E; where mu_e - electron mobility, tau_e - electron life time; E - electric field stremgth
LAMH = 0.185 #[cm] same as above for the holes #!!!!0.185 - looks very good
d = D = 0.1 #[cm] - detector spatial size
CdTeden = 5.85 #[g/cm^3]
Eth = 5.

data = np.loadtxt("CdTe_csection.txt", skiprows = 3)
tabs = interp1d(np.log(data[:,0]), np.log(data[:,-1]*CdTeden))
l60 = 1./exp(tabs(log(60.e-3)))
Kq = 8.4/(LAME*(1. - exp(-(D - l60))/LAME) + LAMH*(1. - exp(-l60/LAMH))) # just a guess



def createnew():
    einlist, dflist, eoutlist, xylist, mask = pickle.load(open("events60.pickle", "rb"))

    Ein = np.logspace(log10(Eth), 2.3, 1025)
    # tst = np.linspace(0., 1., 1024)
    # plt.plot(tst, 6.*tst*e*np.exp(-tst))
    # plt.axhline(5.)
    # plt.show()
    print("searching for solution")
    tth = broyden1(lambda t: Eth - Ein*t*e*np.exp(-t), np.ones(Ein.size)/2.)
    print("done")
    # plt.plot(Ein, tth)
    # plt.show()
    RCCR_transfer_function = interp1d(Ein, (tth + 1.)*np.exp(-tth)) # , fill_value = 1.)
    RCCR_den_transfer_function = interp1d(Ein, e*Ein/Eth*tth*(1. - tth))
    s1etmp = np.empty((einlist.size,  48), np.double)
    s1htmp = np.empty((einlist.size,  48), np.double)
    """
    RCCR transfer function define output energy from input signal
    E = Ein e t e^{-t} ; where t is time/\tau - characteristic time for signal e-folds
    Eth = Ein e \tau e^{-tau}
    Eout = Ein (\tau + 1) e^{-\tau}

    there is usual gaussian dispersion of signal taking place so dN/dEin = 1/(\sqrt(2 \pi) \sigma) e^{-(E_in - E_0)^2/2/\sigma^2}
    but we want to find dN/dEout = dEin/dEout dN/dEin
    dEout/dEin = (tau + 1) e^{-\tau} + Ein tau' e^{-tau} - Ein(\tau + 1)\tau' e^{-\tau} = (tau + 1)e^{-tau} - Ein\tau\tau'e^{-\tau} = (1 + 1/tau)Eth/eEin - Eth/e \tau'
    \tau e^{-\tau} = Eth/eEin
    \tau'e^{-\tau} - \tau \tau' e^{-\tau} = -Eth/eEin^2 => \tau' (1/tau - 1) Eth/eEin = -Eth/eEin^2 => \tau' = -tau/(Ein(1 - tau))

    dEout/dEin = (tau + 1)Eth/(e Ein tau) - Eth*tau/(e Ein (1 - tau)) = Eth ((tau + 1)/tau + tau/(1 - tau)) /(e Ein) = (Eth/(e Ein))/(tau(1 - tau))
    """

    data = np.loadtxt("CdTe_csection.txt", skiprows = 3)
    tabs = interp1d(np.log(data[:,0]), np.log(data[:,-1]*CdTeden))

    print("read events, create lag table")

    random = np.random.uniform(0., 1., dflist.shape)
    # dnorm = np.ones(random.size, np.double)
    # dnorm[dflist < 6e-5] = 5e-5
    # dnorm[dflist > 6e-5] = 1e-3
    enlist = np.empty(eoutlist.shape, np.double)
    enlist[:, 0] = einlist
    enlist[:,1:] = einlist[:, np.newaxis] - eoutlist[:, :-1]
    enlist[np.logical_not(mask)] = 1e-3
    enlist = np.maximum(enlist, 1e-3)

    print("create list of n-th order of events, assuming all of them are photones")

    y = (xylist%192)*(3./192.)
    x = (xylist//192)*(3./192.)

    print(mask.sum(), "total number of events", mask.shape)
    """
    plt.hist(x[mask], 64)
    plt.show()
    plt.hist2d(x[mask], y[mask], [np.linspace(0., 3., 49), np.linspace(0., 3., 49)])
    plt.show()
    """
    tga = np.ones(eoutlist.shape, np.double)
    """
    P(x) = 1 - e^{-x/lam} -> f(x) = dP/dx = 1/lam e^(-x/lam)
    dP(x + dy/dx dx) = dP(y + dy) = f(x) dx = f(y(x)) dy/dx dx -> f(x) = f(y(x)) dy/dx  f(x) = 1/lam e^{-x/lam) f(y) = const dy/dx = 1/lam e^{-x/lam} y = e^{-x/lam}
    limits we would like to put x in the range (0, d) -> P(x) = (1. - e^{-x/lam})/(1. - e^{-d/lam}) y(x = 0) = 0 y(x = d) = 1
    f(x) = e^{-x/lam}/(1. - e^{-d/lam}) = 1. dy/dx - > y|_{0}^{Y} = lam*(e^{-x/lam} - 1)/(1. - e^{-x/lam}) -> x =
    """
    mod = np.sqrt((x[:, 1:] - x[:, :-1])**2. + (y[:, 1:] - y[:, :-1])**2. + np.maximum(np.abs(dflist[:, 1:] - dflist[:, :-1]), 1e-5)**2.)/np.maximum(np.abs(dflist[:, 1:] - dflist[:, :-1]), 1e-5)
    print(np.max(np.abs(dflist[:,1:] - dflist[:,:-1])))
    print("arr shapes", tga.shape, eoutlist.shape, x.shape, y.shape, dflist.shape, mod.shape)
    tga[:, 1:] = mod
    # plt.hist(tga[mask], 128)
    # plt.show()
    # tga = tga[mask]
    #enlist[np.logical_not(mask)] = 10.
    # enlist = enlist[mask]

    d =  np.empty(enlist.shape, np.double)
    m2 = dflist < 6e-3
    print(m2.shape, d.shape, dflist.shape)
    d[m2] = 5e-5
    d[np.logical_not(m2)] = 1e-3


    lam = 1./np.exp(tabs(np.log(enlist))) # photon free path in the CdTe
    print(np.min(lam), np.max(lam))
    print(np.min(tga), np.max(tga))
    random = -lam*np.log(1. - (1. - np.exp(-d*tga/lam))*random)/tga #exponential distributed random number from 0 to d
    random[:,:] -= d/2. #shift it d/2 back since each bin cordinate is at its center
    random[:,1:][dflist[:,1:] < dflist[:, :-1]] *= -1. # revert values if photones moves to the upper edge of the detector

    print("check if obtained depth are in the desired range", np.all(np.abs(random) < d))
    print(dflist.dtype, random.dtype, dflist.shape, random.shape)
    dflist += random # add random location of the photones in side the layer

    ebins = np.linspace(0.5e-3, 66.e-3, 1025) #bins for output energy
    shifts = ebins[:257] - ebins[128] #set of the energetic bins to smooth with electronics energy resolution

    """
    Hecht rule
        ::math::
            n_0 (\lam_e/d (1 - exp(-(d - z)/lam_e)) + lam_h/d*(1 - exp(-z/lam_h))
            n(z) = n0 exp(-z/lam)   dn/dz = -1/lam exp(-z/lam)

            now I want to define coeficcient transfering E_in in to q
    """

    eoutlist[np.logical_not(mask)] = 0.
    s1e = Kq*eoutlist*(LAME*(1. - np.exp(-(D - dflist)/LAME)) + LAMH*(1. - np.exp(-dflist/LAMH)))
    print("obtained from electrones and holles charges")
    # note that both e and h strips have their own RCCR multipliers


    xlist = np.apply_along_axis(np.digitize, 1, x, np.linspace(0., 3., 49)) - 1
    ylist = np.apply_along_axis(np.digitize, 1, y, np.linspace(0., 3., 49)) - 1
    print("xlist.sahpe", xlist.shape, xlist.max(), xlist.min())
    print("s1etmp", s1etmp.shape)
    print("s1e", s1e.shape)

    s1etmp[:, :] = 0.
    s1htmp[:, :] = 0.

    tstart = time()
    for i, j in np.ndindex(s1e.shape):
        s1etmp[i, xlist[i, j]] += s1e[i, j]
        s1htmp[i, ylist[i, j]] += s1e[i, j]
        # print(s1etmp[i, :].shape, xlist[i, :].shape, s1e[i, :].shape)
        """
        np.add.at(s1etmp[i, :], xlist[i, :], s1e[i,:])
        np.add.at(s1htmp[i, :], ylist[i, :], s1h[i,:])
        """
        """
        s1etmp[i, j] += s1e[i, xlist[i, j]]
        s1htmp[i, j] += s1h[i, ylist[i, j]]
        """

    print(s1etmp)

    print("collecting charge loop size:", einlist.size, " time:" , time() - tstart)

    empos = np.argmax(s1etmp, axis=1)
    hmpos = np.argmax(s1htmp, axis=1)
    print(empos.size, hmpos.size, hmpos.max(), empos.max())
    m2 = np.logical_and(np.logical_and(empos != 0, empos != 47), np.logical_and(hmpos != 0, hmpos != 47))
    empos, hmpos = empos[m2], hmpos[m2]
    print(empos.size, hmpos.size, m2.shape)
    print("max min position", empos.shape, hmpos.shape, s1etmp.shape, empos.min(), empos.max(), hmpos.min(), hmpos.max())
    print(s1etmp.shape, m2.shape)
    s1etmp = s1etmp[m2]
    s1htmp = s1htmp[m2]
    print("max min position", empos.shape, hmpos.shape, s1etmp.shape, empos.min(), empos.max(), hmpos.min(), hmpos.max())
    #s1qe, s1qh = s1etmp[mask], s1htmp[mask]
    #s1qe = np.empty((empos.shape[0], 3), np.double)
    #s1qh = np.empty((empos.shape[0], 3), np.double)
    sqe = s1etmp[np.tile(np.arange(empos.size), (3, 1)), np.array([empos - 1, empos, empos + 1])]
    sqh = s1etmp[np.tile(np.arange(hmpos.size), (3, 1)), np.array([hmpos - 1, hmpos, hmpos + 1])]
    print(sqe)
    print(sqe.shape)
    print(s1etmp.shape)
    #plt.hist(sqe.sum(axis = 0), 64)
    #plt.show()
    print(sqe.shape, sqh.shape)

    # s1mod = RCCR_den_transfer_function(s1)
    # I want that after all this actions s1mod == Ein, so I need to adjust this
    # k coefficeint somehow
    sigmae = (0.4813e-3 + 0.172e-3*np.sqrt(sqe/1e-3))
    sigmah = (0.4813e-3 + 0.172e-3*np.sqrt(sqh/1e-3))

    de = shifts[1] - shifts[0]

    spec = np.zeros(ebins.size - 1, np.double)

    for j in range(257):
        mask = sqe[1] - shifts[j] > 5e-3
        sqe = sqe[:, mask]
        sqh = sqh[:, mask]
        sigmae = sigmae[:, mask]
        sigmah = sigmah[:, mask]


        twght = RCCR_den_transfer_function((sqe[1] - shifts[j])*1e3)
        print(twght)
        print(twght.max(), twght.min())
        tfunc = RCCR_transfer_function((sqe[1] - shifts[j])*1e3)
        #print(tfunc)

        #plt.hist(sqe[1], 64, histtype="step")
        #plt.hist(sqe.sum(axis = 0), 64, histtype="step")
        #plt.hist(sqe.sum(axis = 0), 64, histtype="step")

        weights = np.sum(np.exp(-shifts[j]**2./sigmae**2./2.)/sigmae*twght, axis=0) #*twght, axis = 0)

        spec += np.histogram((sqe.sum(axis=0) - shifts[j])*tfunc, ebins, weights = weights)[0]

        #weights = np.sum(np.exp(-shifts[j]**2./sigmah**2./2.)*de/sqrt(2.*pi)/sigmah, axis=0) #*twght, axis = 0)
        #spec += np.histogram(sqh.sum(axis=0)*tfunc, ebins, weights = weights)[0]


    plt.plot((ebins[1:] + ebins[:-1])/2., spec)
    plt.show()

    pickle.dump([ebins, spec], open("test60kevlinespec.pickle", "wb"))

if __name__ == "__main__":
    print("")
    createnew()
