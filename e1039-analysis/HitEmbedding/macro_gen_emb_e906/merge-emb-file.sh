#!/bin/bash
## Script to merge the embedding files made from the E906 data.
## 1000 events * 3500 runs  -->  5000 events * 700 merged files
DIR_MACRO=$(dirname $(readlink -f $BASH_SOURCE))

N_EVT_FILE=5000 # N of events per merged file
RUN_LIST=/data2/production/list/R009/R009_fy2017.list
RUN_N=$(cat $RUN_LIST | wc -l)

DIR_IN=scratch
#DIR_OUT=/pnfs/e1039/scratch/$USER/data_emb_e906
DIR_OUT=/pnfs/e1039/persistent/users/$USER/data_emb_e906

##
## Function
##
function WriteOneFile {
    local -r II=$1
    local -r FN_LIST="$2"
    local -r DIR=$DIR_OUT/$(printf "%04d" $II)
    echo "Write to $DIR/embedding_data.root..."
    mkdir -p $DIR
    hadd -k $DIR/embedding_data.root $FN_LIST &>$DIR/log_hadd.txt
    echo "$FN_LIST" >$DIR/input.txt
}

##
## Main
##
if [ -e $DIR_OUT ] ; then
    echo "DIR_OUT exists."
    if [ "X$1" != 'Xgo' ] ; then
	echo "Give 'go' to clean it up and proceed.  Abort."
	exit
    fi
    echo "  Clean up..."
    rm -rf $DIR_OUT
fi

source /e906/app/software/osg/software/e1039/this-e1039.sh

OUT_I=1
LIST_ROOT=
N_EVT_TOT=0
N_EVT_NOW=0
while read DB_SERVER DB_SCHEMA RUN ROADSET ; do
    RUN6=$(printf '%06d' $RUN)
    DIR_RUN=$DIR_IN/${RUN6:0:2}/${RUN6:2:2}/${RUN6:4:2}/out
    FN_RUN=$DIR_RUN/embedding_data.root
    test -e $FN_RUN || continue
    
    N_EVT=$(tail -n 1 $DIR_RUN/stat.txt)
    (( N_EVT_TOT += N_EVT ))
    (( N_EVT_NOW += N_EVT ))
    printf "Run %6d: %6d %6d %6d\n" $RUN $N_EVT $N_EVT_NOW $N_EVT_TOT
    LIST_ROOT+=" $FN_RUN"

    if [ $N_EVT_NOW -ge $N_EVT_FILE ] ; then
	WriteOneFile $OUT_I "$LIST_ROOT"
	(( OUT_I++ ))
	LIST_ROOT=""
	N_EVT_NOW=0
    fi
done <$RUN_LIST

if [ $N_EVT_NOW -gt 0 ] ; then
    WriteOneFile $OUT_I "$LIST_ROOT"
fi

{
    echo "DIR_MACRO  = $DIR_MACRO"
    echo "N_EVT_FILE = $N_EVT_FILE"
    echo "DIR_IN     = $DIR_IN"
    echo "DIR_OUT    = $DIR_OUT"
    echo "N_EVT_TOT  = $N_EVT_TOT"
    date
} | tee $DIR_OUT/00info.txt
