mkdir output/iron/
for f in ../../../lhe/data/iron/*_hepmc.txt.hepmc
do
    filename=$(basename -- "$f")
    python run_sim.py --sim trimuon --inputfile $filename --inputpath ../../../lhe/data/iron/ --zvertex 490 --displaced -n 5000
    mv output.root output/iron/$filename.root
done
