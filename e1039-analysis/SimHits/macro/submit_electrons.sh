#!/bin/bash                                                                                                                                                                        
FILES=/seaquest/users/cmantill/DarkQuest/lhe/displaced_Aprime_Electrons/*
for f in $FILES
do
    nevents=10000
    basename $f .txt
    filename=$(basename "$f" ".txt")
    #echo $filename
    filejob="${filename}-${nevents}"
    #echo $filejob
    echo "source gridsub.sh ${filejob} 1 1 10000 ${filename} > logs/${filename}.log"
    #source gridsub.sh test3-100 1 1 10000 "Brem_0.06_z500_600_eps_-6.2" >test-grid-3.log
done
