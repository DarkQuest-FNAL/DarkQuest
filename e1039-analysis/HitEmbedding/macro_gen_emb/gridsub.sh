#!/bin/bash
DIR_MACRO=$(dirname $(readlink -f $BASH_SOURCE))

JOB_NAME=$1
DO_SUB=$2
N_JOB=$3
N_EVT=$4 # 10000 for 2 hours
echo "JOB_NAME = $JOB_NAME"
echo "DO_SUB   = $DO_SUB"
echo "N_JOB    = $N_JOB"
echo "N_EVT    = $N_EVT"
if [ $DO_SUB == 1 ]; then
    DIR_DATA=/pnfs/e1039/scratch/$USER/HitEmbedding/data_emb
    DIR_WORK=$DIR_DATA/$JOB_NAME
    ln -nfs $DIR_DATA data # for convenience
else
    DIR_WORK=$DIR_MACRO/scratch/$JOB_NAME
fi

rm -rf         $DIR_WORK
mkdir -p       $DIR_WORK
chmod -R 01755 $DIR_WORK

cd $DIR_MACRO
tar czf $DIR_WORK/input.tar.gz  *.C ../inst

for (( I_JOB = 1; I_JOB <= $N_JOB; I_JOB++ )) ; do
    mkdir -p       $DIR_WORK/$I_JOB/out
    chmod -R 01755 $DIR_WORK/$I_JOB
    cp -p $DIR_MACRO/gridrun.sh $DIR_WORK/$I_JOB
    
    if [ $DO_SUB == 1 ]; then
	CMD="/e906/app/software/script/jobsub_submit_spinquest.sh"
	CMD+=" --expected-lifetime='medium'" # medium=8h, short=3h, long=23h
	CMD+=" -L $DIR_WORK/$I_JOB/log_gridrun.txt"
	CMD+=" -f $DIR_WORK/input.tar.gz"
	CMD+=" -d OUTPUT $DIR_WORK/$I_JOB/out"
	CMD+=" file://$DIR_WORK/$I_JOB/gridrun.sh $N_EVT $I_JOB"
	$CMD |& tee $DIR_WORK/$I_JOB/log_jobsub_submit.txt
	RET_SUB=${PIPESTATUS[0]}
	test $RET_SUB -ne 0 && exit $RET_SUB
    else
	export  CONDOR_DIR_INPUT=$DIR_WORK/$I_JOB/in
	export CONDOR_DIR_OUTPUT=$DIR_WORK/$I_JOB/out
	mkdir -p $DIR_WORK/$I_JOB/in
	cp -p $DIR_WORK/input.tar.gz $DIR_WORK/$I_JOB/in
	cd $DIR_WORK/$I_JOB
	$DIR_WORK/$I_JOB/gridrun.sh $N_EVT $I_JOB |& tee $DIR_WORK/$I_JOB/log_gridrun.txt
    fi
done
