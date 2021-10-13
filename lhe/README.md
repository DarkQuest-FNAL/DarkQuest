# LHE to HEPMC

## Setup 

To install locally HepMC2:
```
cd DarkQuest/lhe/
wget https://hepmc.web.cern.ch/hepmc/releases/hepmc2.06.10.tgz
tar -zxvf hepmc2.06.10.tgz
mkdir HepMC/
cd HepMC/
cp -r ../HepMC-2.06.10/cmake .
cmake ../HepMC-2.06.10/ -DCMAKE_INSTALL_PREFIX=$YOURPATH/DarkQuest/lhe/HepMC -Dmomentum:STRING=GEV -Dlength:STRING=CM
make
make install
```

Then edit Makefile in `src/` to include your HepMC3 bin/include path e.g. replace `-L/Users/cristina/darkquest/DarkQuest/lhe/HepMC/lib/ -lHepMC  -I/Users/cristina/darkquest/DarkQuest/lhe/HepMC/include/`, with your paths.

To setup ROOT (in the spinquest cluster):
```
source /e906/app/software/osg/software/e1039/this-e1039.sh
```

To compile:
```
mkdir bin
cd src/
make
```

Note, for MAC OS (e.g. if you are working locally), you can replace @rpath in `bin/displacedHepmc` with:
```
sudo install_name_tool -change @rpath/libHepMC.4.dylib  /Users/cristina/darkquest/DarkQuest/lhe/HepMC3/lib/libHepMC.4.dylib bin/displacedHepmc
```

## Input data

The location of A' files lhe files should be in the public `data/` directory in the spinquest cluster:
```
/seaquest/users/cmantill/DarkQuest/lhe/data/Aprime_Electrons/
/seaquest/users/cmantill/DarkQuest/lhe/data/Aprime_Muons/
```
and can also be downloaded from dropbox:
https://www.dropbox.com/sh/92xuab9e37gyujt/AAAogZsPCfrtvOTV1rEtAella?dl=0

## Creating HepMC from LHE

- The script `displacedHepmc` creates a HepMC file from the input LHE A' file. The input LHE file is sampled 2000 times 

To run:
```
# create output directories
mkdir output/displaced_Aprime_Muons/
mkdir output/displaced_Aprime_Electrons/

./bin/displacedHepmc inputdata ${mech}_${mass}_${eps}.root ${mech} ${lep} ${rseed} ${eps} ${mass} ${vtx1} ${vtx2}

# for example
./bin/displacedHepmc data/Aprime_Electrons/SeaQuestAprimeToElectronsLHE_Eta_mAp_0.54_GeV.txt Eta_0.54_-7.6.root Eta electron 0 -7.6 0.54 500 600
```
## Acceptance

The ratio of `number of accepted events / number of sampled events` can be used to calculate the acceptance. 