#!/bin/bash

nevents=10000
nmu=1
ifile=Brem_0.36_z500_600_eps_-5.6

if [ -z ${CONDOR_DIR_INPUT+x} ];
  then
    CONDOR_DIR_INPUT=./input;
    echo "CONDOR_DIR_INPUT is initiallized as $CONDOR_DIR_INPUT"
  else
    echo "CONDOR_DIR_INPUT is set to '$CONDOR_DIR_INPUT'";
fi

if [ -z ${CONDOR_DIR_OUTPUT+x} ];
  then
    CONDOR_DIR_OUTPUT=./out;
    mkdir -p $CONDOR_DIR_OUTPUT
    echo "CONDOR_DIR_OUTPUT is initiallized as $CONDOR_DIR_OUTPUT"
  else
    echo "CONDOR_DIR_OUTPUT is set to '$CONDOR_DIR_OUTPUT'";
fi

echo "hello, grid." | tee out.txt $CONDOR_DIR_OUTPUT/out.txt
echo "HOST = $HOSTNAME" | tee -a out.txt $CONDOR_DIR_OUTPUT/out.txt
pwd | tee -a out.txt $CONDOR_DIR_OUTPUT/out.txt

tar -xzvf $CONDOR_DIR_INPUT/input.tar.gz
ls -lh | tee -a out.txt $CONDOR_DIR_OUTPUT/out.txt

FN_SETUP=/e906/app/software/osg/software/e1039/this-e1039.sh
if [ ! -e $FN_SETUP ] ; then # On grid
    FN_SETUP=/cvmfs/seaquest.opensciencegrid.org/seaquest/${FN_SETUP#/e906/app/software/osg/}
fi
echo "SETUP = $FN_SETUP"
source $FN_SETUP pr.69
echo $PWD
#DIR_TOP=/seaquest/users/cmantill/DarkQuest/e1039-analysis/SimHits/
#DIR_TOP=/pnfs/e906/scratch/users/cmantill/SimHits/
#DIR_TOP=$PWD
#echo $DIR_TOP
#export LD_LIBRARY_PATH=$DIR_TOP/install/lib/:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=$DIR_TOP/:$LD_LIBRARY_PATH
#echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"

time root -b -q Fun4Sim.C\(${nevents},\"${ifile}\"\)
mv *.root $CONDOR_DIR_OUTPUT/

echo "gridrun.sh finished!"
