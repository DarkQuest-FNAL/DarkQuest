# SimHits
Collection of code to run e1039-simulation, where 
- `src` is the directory that controls information to be saved in the output root file
- `macro` is the directory that contains some config files that specify the geometry setup and reconstruction steps

## Getting started

### Downloading code
```
git clone git@github.com:cmantill/DarkQuest.git
cd DarkQuest/e1039-analysis/SimHits/
```

SimHits depends on [e1039-core](https://github.com/E1039-Collaboration/e1039-core), and we have our customized version [here](https://github.com/cmantill/e1039-core). You will want to install e1039-core privately (or call the already installed official simulation).

### Installing e1039-core 

<details><summary>Private installation in the spinquest machine (needed for the displaced tracking)</summary>
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
If you have the private version of e1039 installed, change the `MY_E1039` in `setup_mye1039.sh` to your directory; if not, you can leave it unchaged - the default will be
cmantill's version. Then source the script:

```
source setup_mye1039.sh
```

Or to call the official e1039-core version do:
```
source setup.sh
```

The source has to be done everytime you log in.

### Compiling SimAna
To compile the analysis code under `SimAna` (after `setup_mye1039.sh`):
```
mkdir work
mkdir install
cd work
cmake ../src/ -DCMAKE_INSTALL_PREFIX=../install
make
make install
cd ../
```

Later to recompile if you have made changes:
```
cd work
make
make install
```

### To run

The main script `Fun4Sim.C` depends on the location of signal files that depends on `DIR_CMANTILL`, which has been setup in the source script.

To test
```
cd macro
root -b -q Fun4Sim.C 
```

The `Fun4Sim.C` runs on the Aprime to dielectron signal samples; `Fun4Sim_Ap2Muons.C` runs on the Aprime to dimuon signal samples.

To run over a number of signal files, do `source run_electrons.sh $MECH $NEVENTS` e.g.
```
source run_electrons.sh Brem 10000
```
