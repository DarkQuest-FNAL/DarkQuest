#!/bin/bash

FN_LIST=list_dst.txt
FN_INFO=info_dst.txt

{
echo "Command: $0 $*"
echo "Time: $(date '+%F %H:%M:%S')"

NUM_TOT=0
SIZE_TOT=0
for DIR in $* ; do
    DIR_FULL=$(readlink -f $DIR)
    echo "Directory = $DIR_FULL"
    while read NAME SIZE ; do
	printf "  %20s  %6d MB \n" $NAME $(( $SIZE/1024/1024 ))
	echo $DIR_FULL/$NAME >&3
	NUM_TOT=$(( $NUM_TOT + 1 ))
	SIZE_TOT=$(( $SIZE_TOT + $SIZE / 1024 ))
    done < <(find $DIR_FULL -name DST.root -printf '%P %s\n')
done 3>$FN_LIST

echo
echo "Total: $NUM_TOT files, $(( $SIZE_TOT/1024 )) MB"

} | tee $FN_INFO
