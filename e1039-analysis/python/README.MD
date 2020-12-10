# Python scripts for plotting

Most of them are standalone and just take pre-processed simulation files.

Input files for different masses and production modes (for 5m displacement), w emcal info:
```
ssh -KXY cmantill@spinquestgpvm01.fnal.gov
/seaquest/users/cmantill/DarkQuest/e1039-analysis/SimHits/macro/output_electrons_emcal/*.root
```

For most of these notebooks, we only use one data file for Brem/Eta, i.e. one mass point.
You can find these in here:
```
scp @cmslpc.fnal.gov:~cmantill/nobackup/dq/forSergio/sim_eval_Brem_0.033428_z500_600_eps_-5.2.root
scp @cmslpc.fnal.gov:~cmantill/nobackup/dq/forSergio/sim_eval_Eta_0.019564_z500_600_eps_-5.2.root
```