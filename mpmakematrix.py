import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
from time import clock
import pyfits as pf
import re
import pickle

#...oooOOO000OOOooo......oooOOO000OOOooo......oooOOO000OOOooo...#
#Try to open specified file
path = '/srg/a1/work/hard/ART-XC_rmfsim/newgeom_60keV.eventdata' #<- name of file to open
#path = "tevents.dat"

"""
def hecth_smear(Einc, z):
    epsilon = 4.43*1e-3 #eV, energy per electron-hole pair
    N0      = Einc/epsilon #number of pairs from incoming particle
    F       = 0.15 #Fano factor
    variance = np.sqrt(F*N0)
    pairs    = np.random.normal(N0, variance)
#    print N0, pairs,'FWHM< keV ->', 2.35*variance*epsilon
    Efield = 100/0.1 #V/cm
    mete = 2.*1e-3#cm2/V
    mhth = 4.*1e-4#cm2/V
    d = 0.1#cm
    PH = pairs*Efield*(mete*(1-np.exp(-(d-z)/(mete*Efield))) + mhth*(1-np.exp(-z/(mhth*Efield))))
    if PH>0:
        retval =  np.random.normal(PH, (0.015/2.35)*PH)
        return retval
    else:
        return 0
"""

einlist = []
eoutlist = []
doutlist = []
xylist = []
msize = 0

for line in open(path,"r").readlines():
    res = re.findall("gamma (0.\d+) (-?\d.\d+e?[+-]?\d*) (-?\d.\d+e?[+-]?\d*) ([\w+_\d+_d+ \d.\d+e?\+?\-?\d*]+)", line)
    if not len(res):
        print("l and r", line, res)
        continue
    ein, x, y, res = res[0]
    ein = float(ein)
    x, y = float(x), float(y)
    einlist.append(ein)
    res = " " + res

    #print(res)

    #if "thCdTe" in res:
    #    print(res)
    #print(([], [], []) if (not "thCdTe" in res) else zip(*re.findall("thCdTe_(\d+)_(\d+) (\d.\d+e?[+-]?\d*)", res)))
    xy, depth, eout = ([], [], []) if (not "thCdTe" in res) else zip(*re.findall("thCdTe_(\d+)_(\d+) (\d.\d+e?[+-]?\d*)", res))
    #if " CdTe" in res:
    #    print(res)
    xy1, depth1, eout1 = ([], [], []) if (not " CdTe" in res) else zip(*re.findall(" CdTe_(\d+)_(\d+) (\d.\d+e?[+-]?\d*)", res))
    #print(depth, depth1)
    #print([(float(d) + 0.5)*5e-5 for d in depth])
    #print([(float(d) + 5.5)*1e-3 for d in depth1])
    lsize = len(eout) + len(eout1)

    doutlist.append(np.array([(float(d) + 0.5)*5e-5 for d in depth] + [(float(d) + 5.5)*1e-3 for d in depth1]))
    eoutlist.append(np.array([float(e) for e in (list(eout) + list(eout1))]))
    xylist.append(np.array([int(p) for p in (list(xy) + list(xy1))]))

    msize = msize if msize > lsize else lsize

einlist = np.array(einlist)
eflist = np.zeros((len(eoutlist), msize), np.double)
dflist = np.zeros((len(eoutlist), msize), np.double)
xylist1 = np.zeros((len(eoutlist), msize), np.int)
mask = np.zeros((len(eoutlist), msize), np.bool)

for i in range(len(eoutlist)):
    eflist[i, :eoutlist[i].size] = eoutlist[i][:]
    dflist[i, :eoutlist[i].size] = doutlist[i][:]
    xylist1[i, :eoutlist[i].size] = xylist[i][:]
    mask[i, :eoutlist[i].size] = True

pickle.dump([einlist, dflist, eflist, xylist1, mask], open("/srg/a1/work/andrey/srg_detector/BKI_data/events60.pickle", "wb"))
