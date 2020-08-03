#!/bin/bash
FILES=/pnfs/e906/persistent/users/liuk/darkp/digit/R008/02/*/*.root
for f in $FILES
do
    x=$(basename $f)
    ./dataTuple $f $x
done
