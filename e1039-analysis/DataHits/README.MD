# DataHits
Creates simple ROOT TTree with basic Hit information 

To compile
```
source /cvmfs/seaquest.opensciencegrid.org/seaquest/software/current/setup.sh
mkdir build
cd build
cmake ..
make
```

The main executable is `dataTuple.cxx` in the `exe` folder.
This script takes a 2017 data file, and saves basic information of the hits.

2017 data files are located in:
```
/pnfs/e906/persistent/users/liuk/darkp/digit/R008/02/*/digit*.root
```

To run over one file one can do:
```
./dataTuple /pnfs/e906/persistent/users/liuk/darkp/digit/R008/02/44/digit_024484_R008.root digit_024484_R008.root
```

An example can be found in:
```
ssh spinquestgpvm01.fnal.gov
ls -lrth /seaquest/users/cmantill/DarkQuest/e1039-analysis/DataHits/digit_024484_R008.root 
```