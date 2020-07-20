# SimHits
Creates simple ROOT TTree with basic g4SimHit information for input HepMC files.

To get started:
```
git clone https://github.com/cmantill/DarkQuest.git
cd DarkQuest/e1039-analysis/SimHits/
```

To compile:
```
source /e906/app/software/osg/software/e1039/this-e1039.sh
mkdir work
mkdir install
cd work
cmake ../src/ -DCMAKE_INSTALL_PREFIX=../install
make
make install
cd ../
export DIR_TOP=`pwd`
export LD_LIBRARY_PATH=$DIR_TOP/install/lib/:$LD_LIBRARY_PATH
```

To recompile after changes:
```
cd work
make
make install
```

To re-setup after code has compiled (every time) in spinquestgpvm01.fnal.gov:
```
source setup.sh # (or setup_mye1039.sh - if using private version of code)
```

## Current private version 
The current version of the code uses a private version of e1039-core located in: https://github.com/liuk/e1039-core.git

If you want to install the private e1039-core version in your spinquest machine follow:
```
# make directory in your home dir
cd 
mkdir mye1039

# clone repo                                                                                                                                                                      
git clone https://github.com/E1039-Collaboration/e1039-core.git

# add kun's changes                                                                                                                                                              
git remote add kun https://github.com/liuk/e1039-core.git
git fetch kun EMCal
git merge kun/EMCal

# to setup installation
source /seaquest/users/cmantill/mye1039/core-inst/this-e1039.sh

# build                                                                                                                                                                           
./build.sh

# When using an new shell environment (i.e. text terminal) you should source this script (change location correspondingly)
source /seaquest/users/cmantill/mye1039/core-inst/this-e1039.sh
```

If you do not want to install your private version of e1039 you can call cmantill's version by sourcing the following scripts
```
source setup_mye1039.sh
```

The main scripts are located in `macro`.

## In screen
To run a process in `screen` (in SimHits directory): 
```
kinit
cd SimHits/
export DIR_TOP=`pwd`
source setup_mye1039.sh
export LD_LIBRARY_PATH=$DIR_PWD/install/lib/:$LD_LIBRARY_PATH
export DIR_CMANTILL=/seaquest/users/cmantill/DarkQuest/e1039-analysis/SimHits
cd macro/
mkdir -p simeval_electrons_emcal/
mkdir -p output_electrons_emcal/
```

To test
```
root -b -q Fun4SimHit.C 
```

to run over a number of files `source run_electrons.sh $MECH $NEVENTS` e.g.
```
source run_electrons.sh Brem 10000
```