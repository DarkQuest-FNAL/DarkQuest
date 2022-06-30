mkdir output/iron/
for f in /seaquest/users/cmantill/dq/DarkQuest/lhe/data/iron/*_hepmc.txt
do
    filename=$(basename -- "$f")
    python run_sim.py --sim trimuon --inputfile $filename --inputpath /seaquest/users/cmantill/dq/DarkQuest/lhe/output/iron/ --zvertex 490 -n 5000 
    mv output.root output/iron/$filename.root
done
