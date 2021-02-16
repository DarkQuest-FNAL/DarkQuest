#!/bin/bash
mech=$1
#input="${mech}_toSubmit.txt"
input="${mech}_-6.0.txt"
nevents=$2
while IFS= read -r ifile
do
    #root -b -q Fun4SimHit.C\(${nevents},\"${ifile}\"\)
    root -b -q Fun4Sim.C\(${nevents},\"${ifile}\"\)
    mv output.root simana_electrons/${ifile}.root
done < "$input"
