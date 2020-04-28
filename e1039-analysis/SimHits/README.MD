# SimHits
Creates simple ROOT TTree with basic g4SimHit information for input HepMC files.

To compile
```
source /e906/app/software/osg/software/e1039/this-e1039.sh
mkdir work
mkdir install
cmake ../src/ cmake -DCMAKE_INSTALL_PREFIX=../install
cd work
make
make install
export LD_LIBRARY_PATH=`pwd`/install/lib/:$LD_LIBRARY_PATH
```

Or to setup after code has compiled (every time) in spinquestgpvm01.fnal.gov:
```
source setup.sh
```

You can change the HepMC input in `Fun4Sim.C` (probably should be an option) and then run it as:
```
root -l Fun4Sim.C\(10000\)
```
It will look it in your `$DIR_TOP/../../lhe/displaced_*` files.
