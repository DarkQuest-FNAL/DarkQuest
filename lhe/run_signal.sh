#!/bin/bash                                                                                                                                                                        
#FILES=data/Aprime_Electrons/*Brem*
#FILES=data/Aprime_Electrons/*Eta*
FILES=/seaquest/users/cmantill/DarkQuest/lhe/data/Aprime_Muons/*

for f in $FILES
do
    IFS='_' read -ra strarr <<< "$f" 
    mech=${strarr[2]}
    mass=${strarr[4]}
    for logeps in `seq -4.0 -0.2 -7.6`
    do
	    ./bin/displacedHepmc ${f} ${mech}_${mass}_${logeps}.root ${mech} muon 0 ${logeps} ${mass} 200 600
    done
done
