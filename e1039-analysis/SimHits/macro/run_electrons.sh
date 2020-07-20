#!/bin/bash
mech=$1
input="${mech}_toSubmit.txt"
nevents=$2
while IFS= read -r ifile
do
    root -b -q Fun4SimHit.C\(${nevents},\"${ifile}\"\)
done < "$input"
