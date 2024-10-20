#!/bin/bash                                                                                                                                                                        
#FILES=data/Aprime_Electrons/*Brem*
#FILES=data/Aprime_Electrons/*Eta*
#FILES=/cms/data/seaquest/users/cmantill/DarkQuest/lhe/data/Aprime_Muons/*
FILES=/seaquest/users/yfeng/LHEFiles/electrons/eta/*

for f in $FILES
do
    IFS='_' read -ra strarr <<< "$f" 
    #mech=${strarr[2]}
    mech="Eta"
    #mass=${strarr[4]}
    mass=${strarr[1]}
    #for logeps in `seq -4.0 -0.2 -7.6`
    for logeps in `seq -4.0 -0.2 -4.0`		  
    do
	    ./bin/displacedHepmc ${f} ${mech}_${mass}_${logeps}.root ${mech} electron 0 ${logeps} ${mass} 500 600 1
    done
done
