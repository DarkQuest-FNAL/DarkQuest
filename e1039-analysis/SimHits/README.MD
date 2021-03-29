# SimHits
Collection of code to run e1039-simulation.

## Getting started

### Connecting to cluster

### Downloading code
```
git clone git@github.com:cmantill/DarkQuest.git
cd DarkQuest/e1039-analysis/SimHits/
```

You will want to install e1039-core privately (or call the already installed official simulation).

### Installing e1039-core 

<details><summary>Private installation in the spinquest machine</summary>
<p>
  
  ```bash
  # make directory in your home dir
  cd 
  mkdir mye1039
  
  # clone repo
  git clone git@github.com:cmantill/e1039-core.git
  
  # build and install
  cd /path/to/directory_where_you_download_e1039-core
  ./script/setup-install.sh auto
  source ../core-inst/this-e1039.sh
  ./build.sh
  ```
</p>
</details>

### Setting up e1039-core
When using an new shell environment (i.e. text terminal) you should source the `this-e1039.sh` script.
You can change the location correspondingly.
```
source /seaquest/users/cmantill/mye1039/core-inst/this-e1039.sh
```

If you do not want to install your private version of e1039 you can call cmantill's version by sourcing the following scripts
```
source setup_mye1039.sh
```

To call the official e1039-core version do:
```
source setup.sh
```

### Compiling SimAna
The main scripts are located in `macro`.

To compile (after `setup_mye1039.sh`):
```
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
### To run

The main script `Fun4Sim.C` depends on the location of signal files that depends on `DIR_CMANTILL` so make sure to setup that variable as indicated below.

To setup the environment variables (after compiling):
```
cd SimHits/
export DIR_TOP=`pwd`
source setup_mye1039.sh
export LD_LIBRARY_PATH=$DIR_TOP/install/lib/:$LD_LIBRARY_PATH
export DIR_CMANTILL=/seaquest/users/cmantill/DarkQuest/e1039-analysis/SimHits
cd macro/
```

To test
```
root -b -q Fun4Sim.C 
```

To run over a number of signal files, do `source run_electrons.sh $MECH $NEVENTS` e.g.
```
source run_electrons.sh Brem 10000
```
