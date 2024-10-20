# LHE to HEPMC

### Setup  (outside the spinquest cluster)

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


Note, for MAC OS (e.g. if you are working locally), you can replace @rpath in `bin/displacedHepmc` with:
```
sudo install_name_tool -change @rpath/libHepMC.4.dylib  /Users/cristina/darkquest/DarkQuest/lhe/HepMC3/lib/libHepMC.4.dylib bin/displacedHepmc
```

### Setup (on the spinquest cluster)
make changes in `setup_e1039.sh`, namely `MY_E1039` and `HepMC_DIR` if needed. Then run:
```
source setup_e1039.sh
```

To compile:
```
mkdir bin
cd src/
make -j4
```

## Input data

```
# old samples
/seaquest/users/yfeng/LHEFiles/Old/Aprime_Muons
# new samples
/seaquest/users/yfeng/LHEFiles
```

## Creating HepMC from LHE

- The script `displacedHepmc` creates a HepMC file from the input LHE A' file. The input LHE file is sampled 2000 times 

To run:
```
# create output directories
mkdir -p output/displaced_Aprime_Muons/
mkdir -p output/displaced_Aprime_Electrons/

./bin/displacedHepmc inputdata ${mech}_${mass}_${eps}.root ${mech} ${lep} ${rseed} ${eps} ${mass} ${vtx1} ${vtx2} ${isNewFile}

# for example
./bin/displacedHepmc /seaquest/users/yfeng/LHEFiles/electrons/pion/m_0.01.lhe Pion_0.01_-5.6.root Pion electron 0 -5.6 0.01 500 600 1 
```

There is also the script `run_signal_newSamples.sh` available to run the script for all the A' files in the `eta/` channel.
```
bash run_signal_newSamples.sh  > run_eta_electrons.log
```

## Acceptance

The ratio of `number of accepted events / number of sampled events` can be used to calculate the acceptance. 