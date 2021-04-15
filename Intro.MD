# Introduction to DarkQuest

DarkQuest is a proton fixed-target experiment especially sensitive to visible dark sector signatures in MeV-GeV mass range. 
It proposes a few low-cost upgrades to existing detectors from the SeaQuest/SpinQuest (E-1039) Fermilab fixed target beam-dump experiment.

# Theory and Experimental moitvation 

## Dark sector signatures

DarkQuest will directly search for the decays of unstable mediators and dark states to visible Standard Model particles.
It does so by shooting protons into a nuclear target. 
When protons collide with a nucleus, a new dark photon might be produced

DarkQuest distinguishes itself from other experimental initiatives in that it can cover a unique parameter space 
of dark matter scenarios as a compact beam dump spectrometer detector.

## DarkQuest experimental setup

DarkQuest would operate on the neutrino-muon beamline of the Fermilab Accelerator Complex (NM4) and 
receive a high-intensity beam of 120 GeV protons from the Main Injector. 

The next figure shows the location of the experiment in the Fermilab beam-line.

<p align="left">
  <img src="/images/location2.png" width="300" alt="NM4 Location">
</p>

The DarkQuest detector
The detector setup would be a version of the SeaQuest/SpinQuest spectrometer, but with an added Electromagnetic Calorimeter (Ecal).
The spectrometer is made of hodoscopes (scintillator planes) and drift chambers.

Here is the schematic of the DarkQuest experiment with the Ecal in red.

<p align="left">
  <img src="/images/detector.png" width="500" alt="Detector Setup">
</p>

Here:
* `DP` stands for `Dark Photon Hodoscope`. There are two 


# Computing 

## Getting started with computing accounts

If you do not have a Fermilab Kerberos account follow these instructions:

1. Fill out the [Fermilab Access Form for Users](https://get-connected.fnal.gov/users/access/). 

There are detailed instructions on how to get a computing 
account [here](hhttps://uscms.org/uscms_at_work/physics/computing/getstarted/getaccount_fermilab.shtml#ComputingAccountDetails). 
*Note that these have been written for the CMS experiment, so you need to change the affiliation to SpinQuest/E1039.
If your institute is in CMS then you can also get a CMS account first and later a SpinQuest account - see instructions below.*
For point-of-contact please put your supervisor who is associated/user w. Fermilab.

2. After the request has been accepted, you are given a default password for your Kerberos and Services account.

Contact the Service Desk (630-840-2345) to get your default passwords. Follow these links 
to change the [Kerberos password](https://fermi.servicenowservices.com/kb_view.do?sysparm_article=KB0010628) 
and the (Services password)[https://fermi.servicenowservices.com/kb_view.do?sysparm_article=KB0010542].

3. Setting up the Kerberos configuration and connecting to the cluster:

Find detailed instructions [here](https://uscms.org/uscms_at_work/physics/computing/getstarted/uaf.shtml#prerequisiteKerberosConfig) 
(note that again some of these instructions have been written for CMS accounts).
Once you have set up the kerberos configuration you can obtain a ticket by:
```
kinit username@FNAL.GOV
```    
4. Connecting to the spinQuest cluster
```
ssh -KXY username@spinquestgpvm01.fnal.gov
```    

If you already have a Fermilab Kerberos account and want to obtain access to the SpinQuest experiment:

1. Go to the Fermilab service desk to [request computing accounts](https://fermi.servicenowservices.com/wp).
2. Click in `Request Service`. Search for, or find in the list, `Experiment/Project/Collaboration Computing Account`. 
3. Fill in name, Kerberos Principal (username), desired affiliation/experiment, your institution name. By default, keep “No” for Computing/IT Professional Access Only section.
4. Wait for approval! It is recommended that we contact Kun Liu [liuk.pku@gmail.com] in advance so that he knows that the request is coming and can approve faster.
