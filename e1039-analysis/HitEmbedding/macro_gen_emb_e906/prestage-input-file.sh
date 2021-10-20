#!/bin/bash
## Script to prestage the E906 digit file.
DIR_BASE=/pnfs/e906/production/digit/R009
MANIP=/e906/app/software/script/manip-pnfs-file.sh

find $DIR_BASE -mindepth 2 -maxdepth 2 -type d | while read DIR_TGT ; do
    echo "----------------------------------------------------------------"
    echo "DIR_TGT = $DIR_TGT"
    $MANIP -p $DIR_TGT/digit_*.root
done

## A simpler way.
#RUN12=$1
#RUN34=$2
#DIR=$DIR_BASE/$RUN12/$RUN34
#if [ ! -e $DIR ] ; then
#    echo "Cannot find '$DIR'.  Abort."
#fi
#$MANIP -p $DIR/digit_*.root

