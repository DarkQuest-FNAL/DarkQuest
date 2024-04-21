"""
script used to calculate the xsec (i.e., NAp)
given the production mode, coupling, and mass
"""

import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import os, sys
plt.style.use(hep.style.CMS)

INPUTDIR = "/seaquest/users/yfeng/DarkQuest/DarkQuest/lhe"

def loadBR(lep="muons"):
    """
    load the branching ratio
    (only depends on the final states and the mass.
    No dependence on the coupling)
    """
    masses = np.empty(0, dtype=np.float64)
    bfs = np.empty(0, dtype=np.float64)
    if lep == 'muons':
        iBFfile = f"{INPUTDIR}/data/BFtoMuons.txt"
    else:
        iBFfile = f"{INPUTDIR}/data/BFtoElectrons.txt"
    iBR = open(iBFfile)
    # Log[10,mAp/GeV]	Log[10,BF]
    for l in iBR:
        temp = l.split()
        if '#' in l or len(temp) == 0:
            continue
        mass = pow(10, float(temp[0]))
        bf = pow(10, float(temp[1]))
        masses = np.append(masses, mass)
        bfs = np.append(bfs, bf)
    return masses, bfs


def plotBR():
    """
    plot the branching ratio
    """
    colors = ['blue', 'red', 'green', 'purple']
    markers = ['o', 's', '^', 'v']
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for idx, lep in enumerate(['electrons', 'muons']):
        masses, bfs = loadBR(lep)
        ax.plot(masses, bfs, label=lep, color=colors[idx], marker=markers[idx])
    ax.set_xlim([1e-1, 1e1])
    ax.set_ylim([1e-3, 0.5e2])
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('$m_{(A^{\prime})}$ [GeV]')
    ax.set_ylabel('BR ($A^{\prime} \\rightarrow l^{+}l^{-}$)')
    ax.legend()
    plt.show()
    plt.savefig("BR.pdf")
    
    
def loadRates(epsilon, POT = 1.44e18, mech="Brem"):
    """
    load the Ap rates
    does not depend on the lepton type
    """
    masses = np.empty(0, dtype=np.float64)
    rates = np.empty(0, dtype=np.float64)
    
    assert mech in ["Brem", "Eta"], "Error: unknown production mechanism. Choose from Brem or Eta"
    
    ifileYield = f"{INPUTDIR}/data/{mech}Yield.txt"
    if not os.path.exists(ifileYield):
        print(f"Error: {ifileYield} does not exist")
        sys.exit(1)
        
    contents = open(ifileYield)
    for l in contents:
        temp = l.split()
        if '#' in l or len(temp) == 0:
            continue
        mass = pow(10, float(temp[0]))
        xsec = pow(10, float(temp[1]))
        NAp = xsec * (epsilon / 1e-6) ** 2 * (POT / 1.44e18)
        
        masses = np.append(masses, mass)
        rates = np.append(rates, NAp)
    return masses, rates


def getBR(masses_input, lep="muons"):
    """
    get the branching ratio for given masses
    """
    if isinstance(masses_input, float):
        masses_input = [masses_input]
    if isinstance(masses_input, list):
        masses_input = np.array(masses_input)
    return np.interp(masses_input, *loadBR(lep))


def getRates(masses_input, epsilon, POT = 1.44e18, mech="Brem"):
    """
    get the rates for given masses
    """
    if isinstance(masses_input, float):
        masses_input = [masses_input]
    if isinstance(masses_input, list):
        masses_input = np.array(masses_input)
    return np.interp(masses_input, *loadRates(epsilon, POT, mech))


def plotRates(POT = 1.44e18):
    """
    plot the rates for given epsilon
    """
    colors = ['blue', 'red', 'green', 'purple', 'orange', 'black']
    markers = ['o', 's', '^', 'v', 'x', '+']
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    idx = 0
    for mech in ['Brem', 'Eta']:
        for epsilon in [1.0e-6, 1.0e-7, 1.0e-8]:
            masses, rates = loadRates(epsilon, POT, mech)
            ax.plot(masses, rates, label=f"{mech}, $\epsilon={epsilon}$", color=colors[idx], marker=markers[idx])
            
            idx += 1
    ax.set_xlim([1e-2, 1e1])
    ax.set_ylim([1e0, 1e8])
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('$m_{(A^{\prime})}$ [GeV]')
    ax.set_ylabel('Rate ($A^{\prime}$ / 1e18 POT)')
    ax.legend()
    plt.show()
    plt.savefig(f"Rates.pdf")
    
class event:
    def __init__(self,values):
        self._mass = float(values[0])
        self._eps = float(values[1])
        self._ctau = float(values[2])
        self._minvz = float(values[3])
        self._maxvz = float(values[4])
        self._acpE = float(values[5])
    
    def getEps(self):
        return self._eps
    
    def getMass(self):
        return self._mass
    
    def getAccpt(self):
        return self._acpE
    
    def __repr__(self):
        return str([self._acpE,self._mass,self._eps,self._minvz,self._maxvz])

def isfloat(num):
    try:
        float(num)
        return True
    except ValueError:
        return False
    
def parseAccptFile(minVz,maxVz,lep="muons",mech="Brem"):
    fileName = f"{INPUTDIR}/{mech}_{lep}_{minVz}_{maxVz}.txt"
    if not os.path.exists(fileName):
        print(f"Error: {fileName} does not exist")
        sys.exit(1)
        
    file = open(fileName)
    eventsList = []
    for l in file.readlines():
        spStr_org = l.split()
        #print(spStr_org)
        spstr = [sp for sp in spStr_org if isfloat(sp)]
        #print(spstr)
        newEvent = event(spstr)
        #if newEvent.getMech()==mech:
        eventsList.append(newEvent)
    return eventsList
