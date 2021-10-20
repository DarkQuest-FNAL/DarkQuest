#!/bin/bash
DIR_MACRO=$(dirname $(readlink -f $BASH_SOURCE))

FN_LIST_SIG=list_sig_file.txt
FN_LIST_EMB=list_emb_file.txt
test ! -f $FN_LIST_SIG && echo "Cannot find '$FN_LIST_SIG'.  Abort." && exit
test ! -f $FN_LIST_EMB && echo "Cannot find '$FN_LIST_EMB'.  Abort." && exit

DO_OVERWRITE=no
USE_GRID=no
JOB_B=
JOB_E=  # 0 = All available signal and/or embedding files
N_EVT=0 # 0 = All events in each signal+embedding file
OPTIND=1
while getopts ":ogj:e:" OPT ; do
    case $OPT in
	o ) DO_OVERWRITE=yes ;;
        g ) USE_GRID=yes ;;
        j ) JOB_E=$OPTARG ;;
        e ) N_EVT=$OPTARG ;;
    esac
done
shift $((OPTIND - 1))
JOB_NAME=$1
test -z $JOB_NAME  && echo "Need a job name.  Abort." && exit

if [ "${JOB_E%-*}" != "$JOB_E" ] ; then # Contain '-'
    JOB_B=${JOB_E%-*} # Before '-'
    JOB_E=${JOB_E#*-} # After '-'
fi

##
## Check and decide N of files/jobs to be processed
##
LIST_SIG=( $(cat $FN_LIST_SIG) )
LIST_EMB=( $(cat $FN_LIST_EMB) )
N_SIG=${#LIST_SIG[*]}
N_EMB=${#LIST_EMB[*]}
N_BOTH=$(( N_SIG < N_EMB ? N_SIG : N_EMB ))
test -z "$JOB_B" || test $JOB_B -lt 1       && JOB_B=1
test -z "$JOB_E" || test $JOB_E -gt $N_BOTH && JOB_E=$N_BOTH

echo "JOB_NAME     = $JOB_NAME"
echo "DO_OVERWRITE = $DO_OVERWRITE"
echo "USE_GRID     = $USE_GRID"
echo "N_EVT        = $N_EVT"
echo "JOB_B...E    = $JOB_B...$JOB_E"
echo "N_SIG        = $N_SIG"
echo "N_EMB        = $N_EMB"

##
## Prepare and execute the job submission
##
if [ $USE_GRID == yes ]; then
    DIR_DATA=/pnfs/e1039/scratch/$USER/HitEmbedding/data_embedded
    DIR_WORK=$DIR_DATA/$JOB_NAME
    ln -nfs $DIR_DATA data # for convenience
else
    DIR_WORK=$DIR_MACRO/scratch/$JOB_NAME
fi

cd $DIR_MACRO
mkdir -p $DIR_WORK
rm -f    $DIR_WORK/input.tar.gz
tar czf  $DIR_WORK/input.tar.gz  *.C  ../inst

for (( JOB_I = $JOB_B; JOB_I <= $JOB_E; JOB_I++ )) ; do
    FN_SIG=${LIST_SIG[$((JOB_I-1))]}
    FN_EMB=${LIST_EMB[$((JOB_I-1))]}
    DIR_WORK_JOB=$DIR_WORK/$(printf "%04d" $JOB_I)

    if [ -e $DIR_WORK_JOB ] ; then
	echo -n "  DIR_WORK_JOB already exists."
	if [ $DO_OVERWRITE = yes ] ; then
	    echo "  Clean up."
	    rm -rf $DIR_WORK_JOB
	else
	    echo "  Skip."
	    continue
	fi
    fi

    mkdir -p $DIR_WORK_JOB/out
    cp -p $DIR_MACRO/gridrun.sh $DIR_WORK_JOB

    if [ $USE_GRID == yes ]; then
	CMD="/e906/app/software/script/jobsub_submit_spinquest.sh"
	CMD+=" --expected-lifetime='medium'" # medium=8h, short=3h, long=23h
	CMD+=" -L $DIR_WORK_JOB/log_gridrun.txt"
	CMD+=" -f $DIR_WORK/input.tar.gz"
	CMD+=" -f $FN_SIG"
	CMD+=" -f $FN_EMB"
	CMD+=" -d OUTPUT $DIR_WORK_JOB/out"
	CMD+=" file://$DIR_WORK_JOB/gridrun.sh $N_EVT"
	$CMD |& tee $DIR_WORK_JOB/log_jobsub_submit.txt
	RET_SUB=${PIPESTATUS[0]}
	test $RET_SUB -ne 0 && exit $RET_SUB
    else
	export  CONDOR_DIR_INPUT=$DIR_WORK_JOB/in
	export CONDOR_DIR_OUTPUT=$DIR_WORK_JOB/out
	mkdir -p $DIR_WORK_JOB/in
	cp -p $DIR_WORK/input.tar.gz $DIR_WORK_JOB/in
	ln    -s  $FN_SIG            $DIR_WORK_JOB/in
	ln    -s  $FN_EMB            $DIR_WORK_JOB/in
	cd $DIR_WORK_JOB
	$DIR_WORK_JOB/gridrun.sh $N_EVT |& tee $DIR_WORK_JOB/log_gridrun.txt
    fi
done
