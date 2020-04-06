#!/bin/bash                                                                                                                                                                        
FILES=data/Aprime_Electrons/*
for f in $FILES
do
    IFS='_' read -ra strarr <<< "$f" 
    mech=${strarr[2]}
    mass=${strarr[4]}
    for logeps in `seq -5.0 -0.2 -8.0`
    do
	./bin/displacedHepmc ${f} ${mech}_${mass}_${logeps}.root ${mech} electron 0 ${logeps} ${mass} 400 650
    done
done