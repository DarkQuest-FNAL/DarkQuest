# A' lhe studies

## Convert LHE to HepMC

You will need HepMC to run the `displacedHepmc` script.

If you are in `spinquestgpvm01` you can use Dylan's HepMC version here: `LDLIBS=$(ROOTLIBS) -lgsl -lgslcblas -L/seaquest/users/dylrus/HepMC2/lib/ -lHepMC  -I/seaquest/users/dylrus/HepMC2/include/`. You will need to modify this in `src/Makefile`. Also you should probably modify HepMC3 libraries -> HepMC.

If not, you can install locally HepMC3:
```
cd DarkQuest/lhe/
wget https://hepmc.web.cern.ch/hepmc/releases/HepMC3-3.2.1.tar.gz
tar -zxvf HepMC3-3.2.1.tar.gz 
mkdir HepMC3/
cp -r ../HepMC3-3.2.1/cmake HepMC3/
cd HepMC3/
cmake ../HepMC3-3.2.1/ -DCMAKE_INSTALL_PREFIX=$YOURPATH/DarkQuest/lhe/HepMC3
make
make install
```

Then edit Makefile in `src/` to include your HepMC3 bin/include path e.g. replace `-L/Users/cristina/darkquest/DarkQuest/lhe/HepMC3/lib/ -lHepMC  -I/Users/cristina/darkquest/DarkQuest/lhe/HepMC3/include/`, with your paths.

To compile:
```
cd src/
make
```

Note, for OS, you can replace @rpath in `bin/displacedHepmc` with:
```
sudo install_name_tool -change @rpath/libHepMC3.3.dylib /Users/cristina/darkquest/DarkQuest/lhe/HepMC3/lib/libHepMC3.dylib bin/displacedHepmc
```

The location of A' files lhe files should be in data:
```
data/Aprime_Electrons
data/Aprime_Muons
```
and are currently in dropbox
https://www.dropbox.com/sh/92xuab9e37gyujt/AAAogZsPCfrtvOTV1rEtAella?dl=0

To run `displacedHepmc`, this creates a HepMC displaced file sampled 2000 times, run:
```
./bin/displacedHepmc inputdata ${mech}_${mass}_${eps}.root ${mech} ${lep} ${rseed} ${eps} ${mass} ${vtx1} ${vtx2}
```
and customize it ith options e.g.
```
./bin/displacedHepmc data/Aprime_Electrons/SeaQuestAprimeToElectronsLHE_Eta_mAp_0.54_GeV.txt Eta_0.54_-7.6.root Eta electron 0 -7.6 0.54 500 600
```

A short script that runs for all masses and couplings:
```
mkdir displaced_Aprime_Muons
mkdir displaced_Aprime_Electrons
source run_dimuons.sh
source run_dielectrons.sh
```
produces the input to SpinQues simulation framework. It is quite heavy.

The ratio of `number of accepted events / number of sampled events` can be used to calculate the acceptance. 
I currently save the output of `source run_dimuons.sh` e.g. `source run_dimuons.sh > data/output_muons_500_600.txt` and draw contours using the `DQ-acceptance.ipynb` notebook.
