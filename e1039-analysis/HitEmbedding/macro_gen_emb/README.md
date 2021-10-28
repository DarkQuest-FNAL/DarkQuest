# Macro to generate events to be embedded

This macro is to generate background events from which hits are taken and embedded.
A ROOT file that contains a compact TTree object (i.e. `embedding_data.root`) outputted by this macro is used in the hit embedding.

## Alternative

The NIM3 events in the E906 or E1039 data are better (i.e. more realistic) sources of background hits.

## Default Setting of Event Generation

The default setting configured in `Fun4Sim.C` uses
* `PHG4SimpleEventGenerator` to produce two mu+ per event, and
* The standard geometry configuration.

You can modify them as you need.

## Usage

The usage is quite similar to `macro_gen_signal` and thus `SimChainDev`.
You use two shell scripts to execute the macro on the grid.

Typically you first execute the following command to execute the macro on local for test;
```
./gridsub.sh emb_20210823 0 1 100
```

You then execute the following commands to submit grid jobs;
```
source /e906/app/software/script/setup-jobsub-spinquest.sh
./gridsub.sh emb_20210823 1 50 1000
```

