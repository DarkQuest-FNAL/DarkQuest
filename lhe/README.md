# A' lhe studies

## Convert LHE to HepMC

You will need HepMC to run the `displacedHepmc` script.
If you are in `spinquestgpvm01` you can use Dylan's version here: `LDLIBS=$(ROOTLIBS) -lgsl -lgslcblas -L/seaquest/users/dylrus/HepMC2/lib/ -lHepMC  -I/seaquest/users/dylrus/HepMC2/include/`. You will need to modify this in `src/Makefile`.

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

Then edit Makefile in `src/` to include your HepMC3 bin/include path e.g. replace `-L/Users/cristina/darkquest/DarkQuest/lhe/HepMC3/lib/ -lHepMC  -I/Users/cristina/darkquest/DarkQuest/lhe/HepMC3/include/`.
Note, for MAC, you can replace @rpath with:
```
sudo install_name_tool -change @rpath/libHepMC3.3.dylib /Users/cristina/darkquest/DarkQuest/lhe/HepMC3/lib/libHepMC3.dylib bin/displacedHepmc
```

The location of A' files lhe files is in data:
```
data/Aprime_Electrons
data/Aprime_Muons
```

To compile:
```
cd src/
make
```

To create HepMC run:
```
./bin/displacedHepmc inputdata ${mech}_${mass}_${eps}.root ${mech} ${lep} ${rseed} ${eps} ${mass} ${vtx1} ${vtx2}
```
e.g.
```
./bin/displacedHepmc data/Aprime_Electrons/SeaQuestAprimeToElectronsLHE_Eta_mAp_0.54_GeV.txt Eta_0.54_-7.6.root Eta electron 0 -7.6 0.54 500 600
```

For all masses and couplings:
```
mkdir displaced_Aprime_Muons
mkdir displaced_Aprime_Electrons
source run_dimuons.sh
source run_dielectrons.sh
```
these will be the input to SpinQues simulation framework.

To draw contours use `DQ-acceptance.ipynb` notebook.
