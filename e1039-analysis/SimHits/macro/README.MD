# SimChainDev

`Fun4Sim.C` only takes the HepMC input and passes it through the simulation (no tracking and vertexing for now).

You can change the HepMC input in `Fun4Sim.C`:
```
root -b -q Fun4Sim.C\(100,\"Brem_0.06_z500_600_eps_-6.2\"\)
```

# Input

It will look for the HepMC file in the same directory so make sure you copy it there.

In `spinquestgpvm01.fnal.gov` the input HepMC files are located in:
```
/seaquest/users/cmantill/DarkQuest/lhe/displaced_Aprime_Electrons/
```
So:
```
ifile="Brem_0.06_z500_600_eps_-6.2"
cp /seaquest/users/cmantill/DarkQuest/lhe/displaced_Aprime_Electrons/$ifile.txt .
```
before running.

## Grid submission

To run a long job you can use the grid.
More instructions are here: https://github.com/E1039-Collaboration/e1039-wiki/wiki/Submit-jobs-to-the-grid

To setup:
```
source setup-jobsub-spinquest.sh 
```

To test:
```
source gridsub.sh 0 1 10000 Eta_0.54_z500_600_eps_-6.2 > logs/local_Eta_0.54_z500_600_eps_-6.2.log
```

To submit:
```
mkdir logs
source gridsub.sh Eta_0.54_z500_600_eps_-6.2-10000 1 1 10000 Eta_0.54_z500_600_eps_-6.2 > logs/Eta_0.54_z500_600_eps_-6.2.log
```

A submission for all files is summarized in `submit_electrons.sh`.

You can run `jobsub_q --user=$USER` (or its alias jobsub_q_mine) to check the status of your jobs.