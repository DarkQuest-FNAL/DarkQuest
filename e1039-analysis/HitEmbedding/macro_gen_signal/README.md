# Macro to generate signal events

This macro is to generate signal events to be used for the hit embedding.
The DST file (i.e. `DST.root`) outputted by this macro is used next.
The number of signal events per job (=file) should be matched with that of *embedding* events that you will use.
The number used below (5000) is suitable for the E906 NIM3 embedding events.

## Alternative

You can use another package for this step, such as `SimChainDev`.
In that case you are recommended to
* Match the number of events per job (=file) with that of *embedding* events that you will use, 
* Confirm that the geometry configuration is identical, and
* Save only necessary data nodes in DSTs to reduce the size of DSTs, via `Fun4AllDstOutputManager::AddNode()`.

## Default Setting of Event Generation

The default setting configured in `Fun4Sim.C` uses
* The J/psi production via `SQPrimaryParticleGen`,
* The standard geometry configuration, and
* The standard event reconstruction (i.e. tracking & vertexing).

You can modify them as you need.
The reconstruction here makes a result _not_ affected by the hit embedding.
It will be compared to another result reconstructed after the hit embedding.

## Usage

The usage is quite similar to `SimChainDev`.
You use two shell scripts to execute the macro on the grid.

Typically you first execute the following command to execute the macro on local for test,
which runs one job that generates 100 events;
```
./gridsub.sh jpsi_20211011
```

You then execute the following commands to submit grid jobs (`-g`),
which runs 100 jobs (`-j 100`) that generate 5000 events/job (`-e 5000`);
```
source /e906/app/software/script/setup-jobsub-spinquest.sh
./gridsub.sh -g -j 100 -e 5000 jpsi_20211011
```

You might submit more jobs when the statistics are not enough;
```
./gridsub.sh -g -j 101-200 -e 5000 jpsi_20211011
```
