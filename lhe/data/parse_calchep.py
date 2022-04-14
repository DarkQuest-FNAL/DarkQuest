import numpy as np
import glob

mMuon = 0.10566 # muon mass in GeV

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

def generateMuVectors(events, massS0):
    """
    Generate lists of muons from S0 decay in S0 rest frame, then boost back to lab frame
    """
    muPlusPrime = [] # in S0 rest frame
    muMinusPrime = []
    muPlus = [] # in lab frame
    muMinus = []
    for x in range(len(events)):
        ranPVec = getMuPmag(massS0) * np.array(sphereSample())
        muPlusPrime.append(np.array([massS0/2.0, ranPVec[0], ranPVec[1], ranPVec[2]]))
        muMinusPrime.append(np.array([massS0/2.0, -1*ranPVec[0], -1*ranPVec[1], -1*ranPVec[2]]))
        muMinus.append(boostToLabFrame(events, x, massS0, muMinusPrime[x]))
        muPlus.append(boostToLabFrame(events, x, massS0, muPlusPrime[x]))
    return [muPlus, muMinus]

def getMuMinusOutSig(events):
    muminus = []
    for nEvent in range(len(events)):
        px = events[nEvent][3]
        py = events[nEvent][4]
        pz = events[nEvent][5]
        magP2 = px**2 + py**2 + pz**2
        E = (magP2 + mMuon**2)**0.5
        muminus.append([E, px, py, pz])
    return np.array(muminus)

def getMuMinusInSig(events):
    return np.array([ev[1] for ev in events])

def parseHepmc(fname,muon_in_pz,muon_beam,mu_plus,mu_minus):
    import pyhepmc_ng as hep
    with hep.open(fname, "w") as f:
        for ev in range(len(muon_in_pz)):
            evt = hep.GenEvent(hep.Units.GEV, hep.Units.CM)
            evt.event_number = ev
            # 1 => 2 + 3 + 4
            # px      py       pz        e   pdgid status
            p1 = hep.GenParticle((0, 0, muon_in_pz[ev], muon_in_pz[ev]), 1001, 1)
            p2 = hep.GenParticle((muon_beam[ev][1], muon_beam[ev][2], muon_beam[ev][3], muon_beam[ev][0]), -13, 3)
            p3 = hep.GenParticle((mu_plus[ev][1], mu_plus[ev][2], mu_plus[ev][3], mu_plus[ev][0]), -13, 3)
            p4 = hep.GenParticle((mu_minus[ev][1], mu_minus[ev][2], mu_minus[ev][3], mu_minus[ev][0]), +13, 3)
            evt.add_particle(p1)
            evt.add_particle(p2)
            evt.add_particle(p3)
            evt.add_particle(p4)
            
            v1 = hep.GenVertex((1.0, 1.0, 1.0, 1.0));
            v1.add_particle_in (p1)
            v1.add_particle_out(p2)
            v1.add_particle_out(p3)
            v1.add_particle_out(p4)
            evt.add_vertex(v1)
            
            evt.weights = [1.0]
            f.write(evt)
            
    #print(oss)
    
def parseFile(fname):
    with open(fname) as f:
        input_rows = f.read().splitlines() 
    header = input_rows[:12]
    mass = [float(row.split()[-2]) for row in header if 'MASS' in row][0]
    if mass<=0.1:
        print('Cannot process visible decays for masses below muon mass')
        return
    signal_events = [[float(x) for x in row.replace('|','').split()] for row in input_rows[12:]]
    mu_plus,mu_minus = generateMuVectors(signal_events, mass)
    mu_minus_beam_out = getMuMinusOutSig(signal_events)
    mu_minus_beam_in_pz = getMuMinusInSig(signal_events)
    out_fname = fname.replace('.txt','.hepmc')
    parseHepmc(out_fname,mu_minus_beam_in_pz,mu_minus_beam_out,mu_plus,mu_minus)

files = glob.glob("iron/*beam*.txt")
for f in files:
    parseFile(f)
