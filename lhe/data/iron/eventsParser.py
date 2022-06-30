import numpy as np
import math
import numpy
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import random

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# this is how the event files index information:
# for mu- (1) N (2) -> mu- (3) N (4) S0 (5)
# each row is: run number, p1.z, p2.z, p3.x, p3.y, p3.z, p4.x, p4.y, p4.z, p5.x, p5.y, p5.z,
            #Q factor, alpha QCD
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

# read in and store data

# signal events
# 15 GeV Beam
sig15M3text = open("15beam0.3mass.txt", "r") # 0.3 S0 mass
strSig15M3 = sig15M3text.read().replace("|", "").strip()
sig15M3 = list([[ float(x) for x in row ] for row in list([ row.split() for row in strSig15M3.split("\n") ]) ])
sig15M5text = open("15beam0.5mass.txt", "r") # 0.5 S0 mass
strSig15M5 = sig15M5text.read().replace("|", "").strip()
sig15M5 = list([[ float(x) for x in row ] for row in list([ row.split() for row in strSig15M5.split("\n") ]) ])
sig15M7text = open("15beam0.7mass.txt", "r") # 0.7 S0 mass
strSig15M7 = sig15M7text.read().replace("|", "").strip()
sig15M7 = list([[ float(x) for x in row ] for row in list([ row.split() for row in strSig15M7.split("\n") ]) ])
sig15M9text = open("15beam0.9mass.txt", "r") # 0.9 S0 mass
strSig15M9 = sig15M9text.read().replace("|", "").strip()
sig15M9 = list([[ float(x) for x in row ] for row in list([ row.split() for row in strSig15M9.split("\n") ]) ])

# 20 GeV Beam
sig20M3text = open("20beam0.3mass.txt", "r") # 0.3 S0 mass
strSig20M3 = sig20M3text.read().replace("|", "").strip()
sig20M3 = list([[ float(x) for x in row ] for row in list([ row.split() for row in strSig20M3.split("\n") ]) ])
sig20M5text = open("20beam0.5mass.txt", "r") # 0.5 S0 mass
strSig20M5 = sig20M5text.read().replace("|", "").strip()
sig20M5 = list([[ float(x) for x in row ] for row in list([ row.split() for row in strSig20M5.split("\n") ]) ])
sig20M7text = open("20beam0.7mass.txt", "r") # 0.7 S0 mass
strSig20M7 = sig20M7text.read().replace("|", "").strip()
sig20M7 = list([[ float(x) for x in row ] for row in list([ row.split() for row in strSig20M7.split("\n") ]) ])
sig20M9text = open("20beam0.9mass.txt", "r") # 0.9 S0 mass
strSig20M9 = sig20M9text.read().replace("|", "").strip()
sig20M9 = list([[ float(x) for x in row ] for row in list([ row.split() for row in strSig20M9.split("\n") ]) ])
sig20M10text = open("20beam1.0mass.txt", "r") # 0.10 S0 mass
strSig20M10 = sig20M10text.read().replace("|", "").strip()
sig20M10 = list([[ float(x) for x in row ] for row in list([ row.split() for row in strSig20M10.split("\n") ]) ])

# 30 GeV Beam
sig30M3text = open("30beam0.3mass.txt", "r") # 0.3 S0 mass
strSig30M3 = sig30M3text.read().replace("|", "").strip()
sig30M3 = list([[ float(x) for x in row ] for row in list([ row.split() for row in strSig30M3.split("\n") ]) ])
sig30M5text = open("30beam0.5mass.txt", "r") # 0.5 S0 mass
strSig30M5 = sig30M5text.read().replace("|", "").strip()
sig30M5 = list([[ float(x) for x in row ] for row in list([ row.split() for row in strSig30M5.split("\n") ]) ])
sig30M7text = open("30beam0.7mass.txt", "r") # 0.7 S0 mass
strSig30M7 = sig30M7text.read().replace("|", "").strip()
sig30M7 = list([[ float(x) for x in row ] for row in list([ row.split() for row in strSig30M7.split("\n") ]) ])
sig30M9text = open("30beam0.9mass.txt", "r") # 0.9 S0 mass
strSig30M9 = sig30M9text.read().replace("|", "").strip()
sig30M9 = list([[ float(x) for x in row ] for row in list([ row.split() for row in strSig30M9.split("\n") ]) ])

# print("this is number of physical events for background, no cuts:", numPhysEventsBkgnd15)
# print("this is number of physical event for signal, no cuts:", numPhysEventsSig15)

# simulation parameters for tungsten nucleus
mMuon = 0.10566 # muon mass in GeV
beamFlux = 10**14 # muon beam flux
tungDen = 7.874 # target density (g/mL)
tungMM = 55.845 # tungsten molar mass (g/mol)
tungThck = 1.7 # target thickness (cm)
avoNum = 6.022 * 10**23 # avogardros num (mol^-1)
tungDen = tungDen * avoNum / tungMM # tungsten density (nuclei/cm^3)
numEvents = 10000
pbToCmSquared = 1.0 * 10**(-36) # 1 pb = 10^-36 cm^2

# calcHEP cross sections in cm^2
# background
hepCrossBkgnd10 = 8.3241 * 10**2 * pbToCmSquared # eBeam = 10 GeV
hepCrossBkgnd20 = 1.6617 * 10**3 * pbToCmSquared # eBeam = 20 GeV
hepCrossBkgnd30 = 2.3369 * 10**3 * pbToCmSquared # eBeam = 30 GeV

# signal event
hepCrossSig10M3 = 8.9045 * 10**0 * pbToCmSquared # eBeam = 10 GeV, S0 mass = 0.3 GeV signal process
hepCrossSig10M5 = 3.5987 * 10**0 * pbToCmSquared # eBeam = 10 GeV, S0 mass = 0.5 GeV signal process
hepCrossSig10M7 = 1.4843 * 10**0 * pbToCmSquared # eBeam = 10 GeV, S0 mass = 0.7 GeV signal process
hepCrossSig10M9 = 6.0602 * 10**-1 * pbToCmSquared # eBeam = 10 GeV, S0 mass = 0.9 GeV signal process
hepCrossSig20M3 = 1.5692 * 10**1 * pbToCmSquared # eBeam = 20 GeV, S0 mass = 0.3 GeV signal process
hepCrossSig20M5 = 8.1408 * 10**0 * pbToCmSquared # eBeam = 20 GeV, S0 mass = 0.5 GeV signal process
hepCrossSig20M7 = 4.3887 * 10**0 * pbToCmSquared # eBeam = 20 GeV, S0 mass = 0.7 GeV signal process
hepCrossSig20M9 = 2.3747 * 10**0 * pbToCmSquared # eBeam = 20 GeV, S0 mass = 0.9 GeV signal process
hepCrossSig20M10 = 1.6305 * 10**0 * pbToCmSquared
hepCrossSig30M3 = 2.0389 * 10**1 * pbToCmSquared # eBeam = 30 GeV, S0 mass = 0.3 GeV signal process
hepCrossSig30M5 = 1.1725 * 10**1 * pbToCmSquared # eBeam = 30 GeV, S0 mass = 0.5 GeV signal process
hepCrossSig30M7 = 7.0572 * 10**0 * pbToCmSquared # eBeam = 30 GeV, S0 mass = 0.7 GeV signal process
hepCrossSig30M9 = 4.2930 * 10**0 * pbToCmSquared # eBeam = 30 GeV, S0 mass = 0.9 GeV signal process
# numPhysEventsBkgnd10 = hepCrossBkgnd10 * beamFlux * tungDen * tungThck



# functions to boost to S0 rest frame, decay S0, and boost back to lab frame
def getGamma(v):
    return 1.0 / (1 - v**2)**0.5

def getS0Momentum(event, nEvent):
    px = event[nEvent][6]
    py = event[nEvent][7]
    pz = event[nEvent][8]
    vec = [px, py, pz]
    return np.array(vec)

# get magnitude of momentum for outgoing S0 four-vector. S0 is particle 4 in enumeration
def getMagMomentum(event, nEvent):
    vec = getS0Momentum(event, nEvent)
    return (vec[0]**2 + vec[1]**2 + vec[2]**2)**0.5

# getFourVector()[0] gives first component of four vector
def getFourVector(event, nEvent, S0mass):
    E = (getMagMomentum(event, nEvent)**2 + S0mass**2)**0.5
    vec = getS0Momentum(event, nEvent)
    lst = [E, vec[0], vec[1], vec[2]]
    return np.array(lst)

# get speed of S0 rest frame. v = p/E
def getSpeedS0Frame(event, nEvent, S0mass):
    E = (getMagMomentum(event, nEvent)**2 + S0mass**2)**0.5
    momentumVec = getS0Momentum(event, nEvent)
    px = momentumVec[0]
    py = momentumVec[1]
    pz = momentumVec[2]
    lst = [px/E, py/E, pz/E]
    return np.array(lst)

# boost transformation into S0 rest frame
def boostMat(event, nEvent, S0mass):
    S0FrameSpeed = getSpeedS0Frame(event, nEvent, S0mass)
    vx = S0FrameSpeed[0]
    vy = S0FrameSpeed[1]
    vz = S0FrameSpeed[2]
    vLen2 = vx*vx + vy*vy + vz*vz
    gammaV = getGamma(vLen2**0.5)
    boost = [[gammaV, -gammaV*vx, -gammaV*vy, -gammaV*vz],
             [-gammaV*vx, 1 + (gammaV-1)*vx**2/vLen2, (gammaV-1)*vx*vy/vLen2, (gammaV-1)*vx*vz/vLen2],
             [-gammaV*vy, (gammaV-1)*vy*vx/vLen2, 1 + (gammaV-1)*vy**2/vLen2, (gammaV-1)*vy*vz/vLen2],
             [-gammaV*vz, (gammaV-1)*vz*vx/vLen2, (gammaV-1)*vz*vy/vLen2, 1 + (gammaV-1)*vz**2/vLen2]]
    return np.array(boost)

# boost to go back into lab frame
def invBoostMat(event, nEvent, S0mass):
    S0FrameSpeed = getSpeedS0Frame(event, nEvent, S0mass)
    vx = -S0FrameSpeed[0]
    vy = -S0FrameSpeed[1]
    vz = -S0FrameSpeed[2]
    vLen2 = vx*vx + vy*vy + vz*vz
    gammaV = getGamma(vLen2**0.5)
    invBoost = [[gammaV, -gammaV*vx, -gammaV*vy, -gammaV*vz],
             [-gammaV*vx, 1 + (gammaV-1)*vx**2/vLen2, (gammaV-1)*vx*vy/vLen2, (gammaV-1)*vx*vz/vLen2],
             [-gammaV*vy, (gammaV-1)*vy*vx/vLen2, 1 + (gammaV-1)*vy**2/vLen2, (gammaV-1)*vy*vz/vLen2],
             [-gammaV*vz, (gammaV-1)*vz*vx/vLen2, (gammaV-1)*vz*vy/vLen2, 1 + (gammaV-1)*vz**2/vLen2]]
    return np.array(invBoost)

# P' = TP
# P = T^-1 P'

def boostToRestFrame(event, nEvent, S0mass):
     return boostMat(event, nEvent, S0mass).dot(getFourVector(event, nEvent, S0mass))

def boostToLabFrame(event, nEvent, S0mass, vec):
    return invBoostMat(event, nEvent, S0mass).dot(vec)

# get a random point on unit sphere
def rand():
  return np.random.uniform(0, 1)

def sphereSample():
    theta = 2*np.pi* rand()
    phi = 0.5*np.pi * rand()**0.5
    x = np.cos(theta)*np.sin(phi)
    y = np.sin(theta)*np.sin(phi)
    z = np.cos(phi)
    if (rand() < 0.5):
        z = -z
    return [x, y, z]

def getMuPmag(massS0):
    return ((massS0/2.0)**2 - mMuon**2)**0.5

# generate lists of muons from S0 decay in S0 rest frame then boost back to lab frame
def generateMuVectors(event, numEvents, massS0):
    muPlusPrime = [] # in S0 rest frame
    muMinusPrime = []
    muPlus = [] # in lab frame
    muMinus = []
    for x in range(numEvents):
        ranPVec = getMuPmag(massS0) * np.array(sphereSample())
        muPlusPrime.append(np.array([massS0/2.0, ranPVec[0], ranPVec[1], ranPVec[2]]))
        muMinusPrime.append(np.array([massS0/2.0, -1*ranPVec[0], -1*ranPVec[1], -1*ranPVec[2]]))
        muMinus.append(boostToLabFrame(event, x, massS0, muMinusPrime[x]))
        muPlus.append(boostToLabFrame(event, x, massS0, muPlusPrime[x]))
    return [muPlus, muMinus]

# GENERATE MU+ MU- DECAY PAIR IN SIGNAL EVENT - - - - - - - - - - - - - - - - - - - - - -
# do this in global scope so mu plus and mu minus arrays don't use different random values

# 10 GeV Beam
muVecs10Sig3 = generateMuVectors(sig10M3, numEvents, 0.3) # signal event, 10 GeV Beam 0.3 S0 Mass
muPlusVec10Sig3 = muVecs10Sig3[0]
muMinusVec10Sig3 = muVecs10Sig3[1]
muVecs10Sig5 = generateMuVectors(sig10M5, numEvents, 0.5) # signal event, 10 GeV Beam 0.5 S0 Mass
muPlusVec10Sig5 = muVecs10Sig5[0]
muMinusVec10Sig5 = muVecs10Sig5[1]
muVecs10Sig7 = generateMuVectors(sig10M7, numEvents, 0.7) # signal event, 10 GeV Beam 0.7 S0 Mass
muPlusVec10Sig7 = muVecs10Sig7[0]
muMinusVec10Sig7 = muVecs10Sig7[1]
muVecs10Sig9 = generateMuVectors(sig10M9, numEvents, 0.9) # signal event, 10 GeV Beam 0.9 S0 Mass
muPlusVec10Sig9 = muVecs10Sig9[0]
muMinusVec10Sig9 = muVecs10Sig9[1]

# 20 GeV beam
muVecs20Sig3 = generateMuVectors(sig20M3, numEvents, 0.3) # signal event, 20 GeV Beam 0.3 S0 Mass
muPlusVec20Sig3 = muVecs20Sig3[0]
muMinusVec20Sig3 = muVecs20Sig3[1]
muVecs20Sig5 = generateMuVectors(sig20M5, numEvents, 0.5) # signal event, 20 GeV Beam 0.5 S0 Mass
muPlusVec20Sig5 = muVecs20Sig5[0]
muMinusVec20Sig5 = muVecs20Sig5[1]
muVecs20Sig7 = generateMuVectors(sig20M7, numEvents, 0.7) # signal event, 20 GeV Beam 0.7 S0 Mass
muPlusVec20Sig7 = muVecs20Sig7[0]
muMinusVec20Sig7 = muVecs20Sig7[1]
muVecs20Sig9 = generateMuVectors(sig20M9, numEvents, 0.9) # signal event, 20 GeV Beam 0.9 S0 Mass
muPlusVec20Sig9 = muVecs20Sig9[0]
muMinusVec20Sig9 = muVecs20Sig9[1]
muVecs20Sig10 = generateMuVectors(sig20M10, numEvents, 1.0) # signal event, 20 GeV Beam 0.10 S0 Mass
muPlusVec20Sig10 = muVecs20Sig10[0]
muMinusVec20Sig10 = muVecs20Sig10[1]

# 30 GeV beam
muVecs30Sig3 = generateMuVectors(sig30M3, numEvents, 0.3) # signal event, 30 GeV Beam 0.3 S0 Mass
muPlusVec30Sig3 = muVecs30Sig3[0]
muMinusVec30Sig3 = muVecs30Sig3[1]
muVecs30Sig5 = generateMuVectors(sig30M5, numEvents, 0.5) # signal event, 30 GeV Beam 0.5 S0 Mass
muPlusVec30Sig5 = muVecs30Sig5[0]
muMinusVec30Sig5 = muVecs30Sig5[1]
muVecs30Sig7 = generateMuVectors(sig30M7, numEvents, 0.7) # signal event, 30 GeV Beam 0.7 S0 Mass
muPlusVec30Sig7 = muVecs30Sig7[0]
muMinusVec30Sig7 = muVecs30Sig7[1]
muVecs30Sig9 = generateMuVectors(sig30M9, numEvents, 0.9) # signal event, 30 GeV Beam 0.9 S0 Mass
muPlusVec30Sig9 = muVecs30Sig9[0]
muMinusVec30Sig9 = muVecs30Sig9[1]

# generate histograms of cos(theta) for lab frame decayed muons
# cos(theta) = pz/p
def muCosHist(muPlus, muMinus, len):
    muPlusCos = []
    muMinusCos = []
    for x in range(0, len):
        pxPlus = muPlus[x][1]
        pyPlus = muPlus[x][2]
        pzPlus = muPlus[x][3]
        magPPlus = (pxPlus**2 + pyPlus**2 + pzPlus**2)**0.5
        muPlusCos.append(pzPlus/magPPlus)
        pxMinus = muMinus[x][1]
        pyMinus = muMinus[x][2]
        pzMinus = muMinus[x][3]
        magPMinus = (pxMinus**2 + pyMinus**2 + pzMinus**2)**0.5
        muMinusCos.append(pzMinus/magPMinus)
    histPlusArray = numpy.array(muPlusCos)
    histMinusArray = numpy.array(muMinusCos)
    plt.hist(histPlusArray, edgecolor='black', bins=200, alpha = 0.8, label = 'mu plus')
    plt.hist(histMinusArray, edgecolor='black', bins=200, alpha = 0.6, label = 'mu minus')
    plt.xlabel('Cos[Theta]')
    matplotlib. pyplot. ticklabel_format(axis="both", style="", scilimits=None)
    plt.ylabel('Counts')
    plt.title('Histogram of Cos[Theta] for Muons in Decay')
    plt.xlim(0.98, 1)
    plt.legend(loc='upper right')
    plt.show()

def getMuMinusOutSig(event, nEvent):
    px = event[nEvent][3]
    py = event[nEvent][4]
    pz = event[nEvent][5]
    magP2 = px**2 + py**2 + pz**2
    E = (magP2 + mMuon**2)**0.5
    lst = [E, px, py, pz]
    return np.array(lst)

# get outgoing mu- for signal event

# 10 GeV Beam
muMinusOut10Sig3 = [] # 0.3 S0 mass
for x in range(0, numEvents):
    muMinusOut10Sig3.append(getMuMinusOutSig(sig10M3, x))
muMinusOut10Sig5 = [] # 0.5 S0 mass
for x in range(0, numEvents):
    muMinusOut10Sig5.append(getMuMinusOutSig(sig10M5, x))
muMinusOut10Sig7 = [] # 0.7 S0 mass
for x in range(0, numEvents):
    muMinusOut10Sig7.append(getMuMinusOutSig(sig10M7, x))
muMinusOut10Sig9 = [] # 0.9 S0 mass
for x in range(0, numEvents):
    muMinusOut10Sig9.append(getMuMinusOutSig(sig10M9, x))

# 20 GeV Beam
muMinusOut20Sig3 = [] # 0.3 S0 mass
for x in range(0, numEvents):
    muMinusOut20Sig3.append(getMuMinusOutSig(sig20M3, x))
muMinusOut20Sig5 = [] # 0.5 S0 mass
for x in range(0, numEvents):
    muMinusOut20Sig5.append(getMuMinusOutSig(sig20M5, x))
muMinusOut20Sig7 = [] # 0.7 S0 mass
for x in range(0, numEvents):
    muMinusOut20Sig7.append(getMuMinusOutSig(sig20M7, x))
muMinusOut20Sig9 = [] # 0.9 S0 mass
for x in range(0, numEvents):
    muMinusOut20Sig9.append(getMuMinusOutSig(sig20M9, x))
muMinusOut20Sig10 = [] # 1.0 S0 mass
for x in range(0, numEvents):
    muMinusOut20Sig10.append(getMuMinusOutSig(sig20M10, x))

# 30 GeV Beam
muMinusOut30Sig3 = [] # 0.3 S0 mass
for x in range(0, numEvents):
    muMinusOut30Sig3.append(getMuMinusOutSig(sig30M3, x))
muMinusOut30Sig5 = [] # 0.5 S0 mass
for x in range(0, numEvents):
    muMinusOut30Sig5.append(getMuMinusOutSig(sig30M5, x))
muMinusOut30Sig7 = [] # 0.7 S0 mass
for x in range(0, numEvents):
    muMinusOut30Sig7.append(getMuMinusOutSig(sig30M7, x))
muMinusOut30Sig9 = [] # 0.9 S0 mass
for x in range(0, numEvents):
    muMinusOut30Sig9.append(getMuMinusOutSig(sig30M9, x))


def getEnergy(px, py, pz, mass):
    magP2 = px**2 + py**2 + pz**2
    return (magP2 + mass**2)**0.5

# get muon four vectors from background event in enumeration order m n -> m m M n
def getMuOutVecsBkgnd(event):
    length = len(event)
    muMinus1Out = []
    muMinus2Out = []
    muPlusOut = []
    for x in range(length):
        px1 = event[x][3]
        py1 = event[x][4]
        pz1 = event[x][5]
        muMinus1Out.append(np.array([getEnergy(px1, py1, pz1, mMuon), px1, py1, pz1]))
        px2 = event[x][6]
        py2 = event[x][7]
        pz2 = event[x][8]
        muMinus2Out.append(np.array([getEnergy(px2, py2, pz2, mMuon), px2, py2, pz2]))
        px3 = event[x][9]
        py3 = event[x][10]
        pz3 = event[x][11]
        muPlusOut.append(np.array([getEnergy(px3, py3, pz3, mMuon), px3, py3, pz3]))
    return [muMinus1Out, muMinus2Out, muPlusOut]

# m^2 = (P1 + P2)^2 Lorentz Dot Productd
def getInvariantMassPair(mu1, mu2, nEvent):
    E = mu1[nEvent][0] + mu2[nEvent][0]
    px = mu1[nEvent][1] + mu2[nEvent][1]
    py = mu1[nEvent][2] + mu2[nEvent][2]
    pz = mu1[nEvent][3] + mu2[nEvent][3]
    return E**2 - px**2 - py**2 - pz**2

# sort background muon pairs by which pair has the higher energy mu-
# returns ordered pair of (higher energy pair, lower energy pair)
def sortMuons(event1, event2, nEvent):
    mu1Energy = event1[nEvent][0]
    mu2Energy = event2[nEvent][0]
    if mu1Energy >= mu2Energy:
        return [event1, event2]
    else:
        return [event2, event1]

def getInvariantMassArrays(muPlusDecaySig, muMinusDecaySig, muMinusOutSig, bkgndEvent):
    length = len(muPlusDecaySig)
    invMassHigherEnergyPairSignal = []
    invMassLowerEnergyPairSignal = []
    invMassHigherEnergyPairBkgnd = []
    invMassLowerEnergyPairBkgnd = []
    mu1Outbkgnd = getMuOutVecsBkgnd(bkgndEvent)[0]
    mu2Outbkgnd = getMuOutVecsBkgnd(bkgndEvent)[1]
    muPlusOutbkgnd = getMuOutVecsBkgnd(bkgndEvent)[2]
    for x in range(0, length):
        invMassHigherEnergyPairSignal.append(getInvariantMassPair(muPlusDecaySig, sortMuons(muMinusDecaySig, muMinusOutSig, x)[0], x))
        invMassLowerEnergyPairSignal.append(getInvariantMassPair(muPlusDecaySig, sortMuons(muMinusDecaySig, muMinusOutSig, x)[1], x))
        invMassHigherEnergyPairBkgnd.append(getInvariantMassPair(sortMuons(mu1Outbkgnd, mu2Outbkgnd, x)[0], muPlusOutbkgnd, x))
        invMassLowerEnergyPairBkgnd.append(getInvariantMassPair(sortMuons(mu1Outbkgnd, mu2Outbkgnd, x)[1], muPlusOutbkgnd, x))
    return [invMassHigherEnergyPairSignal, invMassLowerEnergyPairSignal, invMassHigherEnergyPairBkgnd, invMassLowerEnergyPairBkgnd]


def invariantMassHist(muPlusDecaySig, muMinusDecaySig, muMinusOutSig, bkgndEvent, massS0, eBeam):
    histPair1Array = getInvariantMassArrays(muPlusDecaySig, muMinusDecaySig, muMinusOutSig, bkgndEvent)[0]
    histPair2Array = getInvariantMassArrays(muPlusDecaySig, muMinusDecaySig, muMinusOutSig, bkgndEvent)[1]
    histPair1Bkgnd = getInvariantMassArrays(muPlusDecaySig, muMinusDecaySig, muMinusOutSig, bkgndEvent)[2]
    histPair2Bkgnd = getInvariantMassArrays(muPlusDecaySig, muMinusDecaySig, muMinusOutSig, bkgndEvent)[3]
    plt.hist(histPair1Array, edgecolor='black',bins=1000, alpha = 0.7, label = 'Higher Mu- Energy Pair Signal Process')
    plt.hist(histPair2Array, edgecolor='black', bins=1000, alpha = 0.7, label = 'Lower Mu- Energy Pair Signal Process')
    plt.hist(histPair1Bkgnd, edgecolor='black', bins=1000, alpha = 0.7, label = 'Higher Mu- Energy Pair Background Process')
    plt.hist(histPair2Bkgnd, edgecolor='black', bins=1000, alpha = 0.7, label = 'Lower Mu- Energy Pair Background Process')
    plt.xlim(massS0**2 - 0.05, massS0**2 + .05)
    plt.xlabel('Invariant Mass Squared [GeV^2]')
    plt.ylabel('Counts')
    matplotlib. pyplot. ticklabel_format(axis="both", style="", scilimits=None)
    plt.title('Histogram of Invariant Mass Squared for Muon Pairs in Signal and Background Processes for MS0 mass '
     + str(massS0) + ' [GeV] and beam energy ' + str(eBeam) + ' [GeV]' )
    # plt.title('Histogram of Invariant Mass Squared for Muon Pairs in Signal and Background Processes for MS0 mass of %i GeV, %j GeV Beam' % massS0 % eBeam)
    plt.legend(loc='upper right')
    plt.show()

# invariantMassHist(muPlusVec10Sig9, muMinusVec10Sig9, muMinusOut10Sig9, bkgnd10, 0.9, 10)

# function to find number of events that pass the invariant mass cut
# we want at least one muon pair to be inside the mass window
def getNumInvMassPassCutEvents(muPlusDecaySig, muMinusDecaySig, muMinusOutSig, bkgndEvent, S0Mass, radius):
    length = len(muPlusDecaySig)
    min = S0Mass**2 - radius
    max = S0Mass**2 + radius
    numSigInRange = 0.0
    numBkgndInRange = 0.0
    higherESignal = getInvariantMassArrays(muPlusDecaySig, muMinusDecaySig, muMinusOutSig, bkgndEvent)[0]
    lowerESignal = getInvariantMassArrays(muPlusDecaySig, muMinusDecaySig, muMinusOutSig, bkgndEvent)[1]
    higherEBackground = getInvariantMassArrays(muPlusDecaySig, muMinusDecaySig, muMinusOutSig, bkgndEvent)[2]
    lowerEBackground = getInvariantMassArrays(muPlusDecaySig, muMinusDecaySig, muMinusOutSig, bkgndEvent)[3]
    for x in range(0, length):
        # look at signal events first
        # this condition says if both muon pairs are outside the range
        if (higherESignal[x] < min or higherESignal[x] > max) and (lowerESignal[x] < min or lowerESignal[x] > max):
            numSigInRange = numSigInRange
        else:
            numSigInRange = numSigInRange + 1
        # look at background events next
        if (higherEBackground[x] < min or higherEBackground[x] > max) and (lowerEBackground[x] < min or lowerEBackground[x] > max):
            numBkgndInRange = numBkgndInRange
        else:
            numBkgndInRange = numBkgndInRange + 1
    return [numSigInRange, numBkgndInRange]

# getNumInvMassPassCutEvents(muPlusVecSig15, muMinusVecSig15, muMinusOutSig15, bkgnd15, 0.5, 0.1))

def plotNPassInvtMassCut(muPlusDecaySig, muMinusDecaySig, muMinusOutSig, bkgndEvent, S0Mass, eBeam, hepCrossSig, hepCrossBkgnd, nEvents):
    radius = []
    nPassCutSig = []
    nPassCutBkgnd = []
    for x in range(0, 10):
        scaledRadius = x*.01
        radius.append(scaledRadius/S0Mass**2)
        numSigInRange = getNumInvMassPassCutEvents(muPlusDecaySig, muMinusDecaySig, muMinusOutSig, bkgndEvent, S0Mass, scaledRadius)[0]
        numBackInRange = getNumInvMassPassCutEvents(muPlusDecaySig, muMinusDecaySig, muMinusOutSig, bkgndEvent, S0Mass, scaledRadius)[1]
        sigCrossAcc = hepCrossSig * (numSigInRange/nEvents)
        bkgndCrossAcc = hepCrossBkgnd * (numBackInRange/nEvents)
        nPassCutSig.append(sigCrossAcc * beamFlux * tungDen * tungThck)
        nPassCutBkgnd.append((bkgndCrossAcc * beamFlux * tungDen * tungThck)**0.5)
    plt.plot(radius, nPassCutSig)
    plt.plot(radius, nPassCutBkgnd)
    plt.xlabel('Radius of mass window / MS0^2')
    matplotlib. pyplot. ticklabel_format(axis="both", style="sci", scilimits=None)
    plt.ylabel('Number of physical events that pass the invariant mass cut')
    plt.title('Number of Physical Events that Pass Invariant Mass Cut as Function of Mass Window for MS0 ' + str(S0Mass)
    + ' [GeV] and beam energy ' + str(eBeam) + ' [GeV]')
    plt.legend(loc='upper right')
    matplotlib.pyplot.legend(["Num physical signal events", "Sqrt Num physical background events"])
    plt.show()

# plotNPassInvtMassCut(muPlusVec10Sig3, muMinusVec10Sig3, muMinusOut10Sig3, bkgnd10, 0.3, 10, hepCrossSig10M3, hepCrossBkgnd10, numEvents)
# plotNPassInvtMassCut(muPlusVec10Sig5, muMinusVec10Sig5, muMinusOut10Sig5, bkgnd10, 0.5, 10, hepCrossSig10M5, hepCrossBkgnd10, numEvents)
# plotNPassInvtMassCut(muPlusVec10Sig7, muMinusVec10Sig7, muMinusOut10Sig7, bkgnd10, 0.7, 10, hepCrossSig10M7, hepCrossBkgnd10, numEvents)
# plotNPassInvtMassCut(muPlusVec10Sig9, muMinusVec10Sig9, muMinusOut10Sig9, bkgnd10, 0.9, 10, hepCrossSig10M9, hepCrossBkgnd10, numEvents)
# plotNPassInvtMassCut(muPlusVec20Sig3, muMinusVec20Sig3, muMinusOut20Sig3, bkgnd20, 0.3, 20, hepCrossSig20M3, hepCrossBkgnd20, numEvents)
# plotNPassInvtMassCut(muPlusVec20Sig5, muMinusVec20Sig5, muMinusOut20Sig5, bkgnd20, 0.5, 20, hepCrossSig20M5, hepCrossBkgnd20, numEvents)
# plotNPassInvtMassCut(muPlusVec20Sig7, muMinusVec20Sig7, muMinusOut20Sig7, bkgnd20, 0.7, 20, hepCrossSig20M7, hepCrossBkgnd20, numEvents)
# plotNPassInvtMassCut(muPlusVec20Sig9, muMinusVec20Sig9, muMinusOut20Sig9, bkgnd20, 0.9, 20, hepCrossSig20M9, hepCrossBkgnd20, numEvents)
# plotNPassInvtMassCut(muPlusVec20Sig10, muMinusVec20Sig10, muMinusOut20Sig10, bkgnd20, 1.0, 20, hepCrossSig20M10, hepCrossBkgnd20, numEvents)
# plotNPassInvtMassCut(muPlusVec30Sig3, muMinusVec30Sig3, muMinusOut30Sig3, bkgnd30, 0.3, 30, hepCrossSig30M3, hepCrossBkgnd30, numEvents)
# plotNPassInvtMassCut(muPlusVec30Sig5, muMinusVec30Sig5, muMinusOut30Sig5, bkgnd30, 0.5, 30, hepCrossSig30M5, hepCrossBkgnd30, numEvents)
# plotNPassInvtMassCut(muPlusVec30Sig7, muMinusVec30Sig7, muMinusOut30Sig7, bkgnd30, 0.7, 30, hepCrossSig30M7, hepCrossBkgnd30, numEvents)
# plotNPassInvtMassCut(muPlusVec30Sig9, muMinusVec30Sig9, muMinusOut30Sig9, bkgnd30, 0.9, 30, hepCrossSig30M9, hepCrossBkgnd30, numEvents)
