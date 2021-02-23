#!/bin/bash
mech=$1
input="${mech}_-6.0.txt"
nevents=$2
while IFS= read -r ifile
do
    root -b -q Fun4Sim.C\(${nevents},\"${ifile}\"\)
    mv output.root simana_electrons_Feb23/${ifile}.root
done < "$input"
