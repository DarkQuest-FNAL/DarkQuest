#!/bin/bash
DIR_MACRO=$(dirname $(readlink -f $BASH_SOURCE))

RUN_LIST=/data2/production/list/R009/R009_fy2017.list
RUN_N=$(cat $RUN_LIST | wc -l)
echo "RUN_LIST = $RUN_LIST"
echo "RUN_N    = $RUN_N"

N_EVT=0
RUN_B=85 # run 24279
RUN_E=90
DO_OVERWRITE=no
USE_GRID=no
OPTIND=1
while getopts ":n:e:og" OPT ; do
    case $OPT in
        n )    RUN_E=$OPTARG ;;
        e )    N_EVT=$OPTARG ;;
        o ) DO_OVERWRITE=yes ;;
        g ) USE_GRID=yes     ;;
        * ) PrintHelp ; exit ;;
    esac
done
shift $((OPTIND - 1))

if [ "${RUN_E%-*}" != "$RUN_E" ] ; then # Contain '-'
    RUN_B=${RUN_E%-*} # Before '-'
    RUN_E=${RUN_E#*-} # After '-'
fi
test -z "$RUN_B" || test $RUN_B -lt 1      && RUN_B=1
test -z "$RUN_E" || test $RUN_E -gt $RUN_N && RUN_E=$RUN_N
echo "N_EVT        = $N_EVT"
echo "RUN_B...E    = $RUN_B...$RUN_E"
echo "DO_OVERWRITE = $DO_OVERWRITE"
echo "USE_GRID     = $USE_GRID"

if [ $USE_GRID == yes ] ; then
    #DIR_WORK=/pnfs/e1039/persistent/users/$USER/data_emb_e906
    DIR_WORK=/pnfs/e1039/scratch/$USER/HitEmbedding/data_emb_e906
    ln -nfs $DIR_WORK data # for convenience
else
    DIR_WORK=$DIR_MACRO/scratch
fi
echo "DIR_WORK = $DIR_WORK"

cd $DIR_MACRO
mkdir -p $DIR_WORK
rm -f    $DIR_WORK/input.tar.gz
tar czf  $DIR_WORK/input.tar.gz  *.C *.txt ../inst

RUN_I=0
while read DB_SERVER DB_SCHEMA RUN ROADSET ; do
    (( RUN_I++ ))
    test $RUN_I -lt $RUN_B && continue
    test $RUN_I -gt $RUN_E && break
    test $RUN -lt 24279 && continue # No good spill before run 24279 at present
    echo "----------------------------------------------------------------"
    printf "Run %6d:  %5d / %5d\n" $RUN $RUN_I $RUN_N

    RUN6=$(printf '%06d' $RUN)
    DIR_ROOT=/pnfs/e906/production/digit/R009/${RUN6:0:2}/${RUN6:2:2}
    FN_ROOT=digit_${RUN6}_009.root
    
    DIR_WORK_RUN=$DIR_WORK/${RUN6:0:2}/${RUN6:2:2}/${RUN6:4:2}
    #echo "  DIR_WORK_RUN = $DIR_WORK_RUN"
    if [ -e $DIR_WORK_RUN ] ; then
	echo -n "  DIR_WORK_RUN already exists."
	if [ $DO_OVERWRITE = yes ] ; then
	    echo "  Clean up."
	    rm -rf $DIR_WORK_RUN
	else
	    echo "  Skip."
	    continue
	fi
    fi

    mkdir -p $DIR_WORK_RUN/out
    cp -p $DIR_MACRO/gridrun.sh $DIR_WORK_RUN
    
    if [ $USE_GRID == yes ] ; then
	CMD="/e906/app/software/script/jobsub_submit_spinquest.sh"
	CMD+=" --expected-lifetime='short'" # medium=8h, short=3h, long=23h
	CMD+=" -L $DIR_WORK_RUN/log_gridrun.txt"
	CMD+=" -f $DIR_WORK/input.tar.gz"
	CMD+=" -f $DIR_ROOT/$FN_ROOT"
	CMD+=" -d OUTPUT $DIR_WORK_RUN/out"
	CMD+=" file://$DIR_WORK_RUN/gridrun.sh $RUN $FN_ROOT $N_EVT"
	$CMD |& tee $DIR_WORK_RUN/log_jobsub_submit.txt
	RET_SUB=${PIPESTATUS[0]}
	test $RET_SUB -ne 0 && exit $RET_SUB
    else
	export  CONDOR_DIR_INPUT=$DIR_WORK_RUN/in
	export CONDOR_DIR_OUTPUT=$DIR_WORK_RUN/out
	mkdir -p $DIR_WORK_RUN/in
	cp -p $DIR_WORK/input.tar.gz $DIR_WORK_RUN/in
	ln -nfs $DIR_ROOT/$FN_ROOT $DIR_WORK_RUN/in
	cd $DIR_WORK_RUN
	$DIR_WORK_RUN/gridrun.sh $RUN $FN_ROOT $N_EVT |& tee $DIR_WORK_RUN/log_gridrun.txt
    fi
done <$RUN_LIST

