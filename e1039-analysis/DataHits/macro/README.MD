# DataToDST

To compile AnaModule. You can use eitehr private or shared version of e1039-core code.
```
cd macro/
source setup.sh # (sets up shared version of e1039-core)
source setup_mye1039.sh
mkdir work
mkdir install
cd work
cmake ../AnaModule/ -DCMAKE_INSTALL_PREFIX=../install
make
make install
cd ../
# then set up the library
export DIR_TOP=`pwd`
export LD_LIBRARY_PATH=$DIR_TOP/install/lib/:$LD_LIBRARY_PATH
# these last 2 commands are summarized in setup.sh or setup_mye1039.sh
```

To recompile after changes in the module:
```
cd work
make
make install
```

To re-setup after code has compiled (every time) in spinquestgpvm01.fnal.gov:
```
source setup.sh # (or setup_mye1039.sh - if using private version of code)
```

To test (from digit.root to reconstruction). This file should output both a .DST file and a flat .root file.
```
root -b -q RecoE906Data.C 
```

To flatten a .DST file with the AnaModule you can use
```
root -b -q DSTToTree.C
```