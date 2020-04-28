# A' lhe studies

## Convert LHE to HepMC

You will need HepMC to run the `displacedHepmc` script.

If you are in `spinquestgpvm01` you can use Dylan's HepMC version here: `LDLIBS=$(ROOTLIBS) -lgsl -lgslcblas -L/seaquest/users/dylrus/HepMC2/lib/ -lHepMC  -I/seaquest/users/dylrus/HepMC2/include/`. You will need to modify this in `src/Makefile`.

If not, you can install locally HepMC2:
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

To setup ROOT in spinquest dir:
```
source /e906/app/software/osg/software/e1039/this-e1039.sh
```

To compile:
```
mkdir bin
cd src/
make
```

Note, for OS, you can replace @rpath in `bin/displacedHepmc` with:
```
sudo install_name_tool -change @rpath/libHepMC.4.dylib  /Users/cristina/darkquest/DarkQuest/lhe/HepMC3/lib/libHepMC.4.dylib bin/displacedHepmc
```

The location of A' files lhe files should be in data:
```
data/Aprime_Electrons
data/Aprime_Muons
```
and are currently in dropbox
https://www.dropbox.com/sh/92xuab9e37gyujt/AAAogZsPCfrtvOTV1rEtAella?dl=0

To run `displacedHepmc`, which creates a HepMC displaced file sampled 2000 times, create the output directories
```
mkdir displaced_Aprime_Muons
mkdir displaced_Aprime_Electrons
```
and then run like
```
./bin/displacedHepmc inputdata ${mech}_${mass}_${eps}.root ${mech} ${lep} ${rseed} ${eps} ${mass} ${vtx1} ${vtx2}
```
customizing it with options such as
```
./bin/displacedHepmc data/Aprime_Electrons/SeaQuestAprimeToElectronsLHE_Eta_mAp_0.54_GeV.txt Eta_0.54_-7.6.root Eta electron 0 -7.6 0.54 500 600
```

There are short scripts that run all masses and couplings:
```
source run_dimuons.sh
source run_dielectrons.sh
```
These produce inputs to the SpinQuest simulation framework. Be warned that the outputs can be quite large!

The ratio of `number of accepted events / number of sampled events` can be used to calculate the acceptance. 
I currently save the output of `source run_dimuons.sh` e.g. `source run_dimuons.sh > data/output_muons_500_600.txt` and draw contours using the `DQ-acceptance.ipynb` notebook.
