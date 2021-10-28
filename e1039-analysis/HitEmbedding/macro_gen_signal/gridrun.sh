#!/bin/bash

JOB_ID=$1
N_EVT=$2

if [ -z "$CONDOR_DIR_INPUT" -o -z "$CONDOR_DIR_OUTPUT" ] ; then
    echo "!ERROR!  CONDOR_DIR_INPUT/OUTPUT is undefined.  Abort."
    exit 1
fi
echo "INPUT  = $CONDOR_DIR_INPUT"
echo "OUTPUT = $CONDOR_DIR_OUTPUT"
echo "HOST   = $HOSTNAME"
echo "PWD    = $PWD"

tar xzf $CONDOR_DIR_INPUT/input.tar.gz

E1039_CORE_VERSION=pr.116
FN_SETUP=/e906/app/software/osg/software/e1039/this-e1039.sh
if [ ! -e $FN_SETUP ] ; then # On grid
    FN_SETUP=/cvmfs/seaquest.opensciencegrid.org/seaquest/${FN_SETUP#/e906/app/software/osg/}
fi
echo "SETUP = $FN_SETUP"
source $FN_SETUP

time root -b -q "Fun4Sim.C($JOB_ID, $N_EVT)"
RET=$?
if [ $RET -ne 0 ] ; then
    echo "Error in Fun4Sim.C: $RET"
    exit $RET
fi

mv  *.root *.tsv  $CONDOR_DIR_OUTPUT

echo "gridrun.sh finished!"
