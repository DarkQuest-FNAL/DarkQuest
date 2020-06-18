#!/bin/bash
#mkdir displaced_Aprime_Muons
while read mass yield
do
    for logeps in `seq -3.0 -0.2 -8.0`
    #for logeps in `seq -5.0 -0.2 -5.4`
    do
	#echo ./bin/displacedHepmc data/Aprime_Muons/SeaQuestAprimeToMuonsLHE_Brem_mAp_${mass}_GeV.txt brem_${mass}_${logeps}.root Brem muon 0 ${logeps} ${mass} 500 600
	./bin/displacedHepmc data/Aprime_Muons/SeaQuestAprimeToMuonsLHE_Brem_mAp_${mass}_GeV.txt brem_${mass}_${logeps}.root Brem muon 0 ${logeps} ${mass} 500 600
    done
done < data/brem_yields_muons.txt

while read mass yield
do
    for logeps in `seq -3.0 -0.2 -8.0`
    #for logeps in `seq -5.0 -0.2 -5.4`
    do
	#echo ./bin/displacedHepmc data/Aprime_Muons/SeaQuestAprimeToMuonsLHE_Eta_mAp_${mass}_GeV.txt eta_${mass}_${logeps}.root Eta muon 0 ${logeps} ${mass} 500 600
	./bin/displacedHepmc data/Aprime_Muons/SeaQuestAprimeToMuonsLHE_Eta_mAp_${mass}_GeV.txt eta_${mass}_${logeps}.root Eta muon 0 ${logeps} ${mass} 500 600
    done
done < data/eta_yields_muons.txt
