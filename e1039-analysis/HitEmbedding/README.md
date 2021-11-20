# A Package to Generate and Analyze Hit-Embedded MC Events

(**Lates Update**: all the essential code should already be moved to `SimHits`. This directory is obsolete)

The hit embedding is a procedure for embedding a set of hits into a simulated physics event (like J/psi).
It is to mimic a real event, which contains not only a physics process but also background hits.
It can be used in studying the effects of background hits, such fake track and efficiency loss.

## Terminology

For events, data or files:

| Term | Meaning |
| ---- | ------- |
| "Signal" events    | (Simulated) physics events                    |
| "Embedding" events | Events which (or whose hits) will be embedded |
| "Embedded" events  | Physics events where hits have been embedded  |

## Prerequisite

This program depends on the E1039 resource+share+core packages.
You are recommended to use "spinquestgpvm01.fnal.gov".
You are expected to be familiar with the procedure for generating and analyzing a set of simulated events, using `SimChainDev` and `AnaSimDst` in `e1039-analysis` for example.

You can check out the repository as usual;
```
cd /path/to/your_working_directory
git clone git@github.com:E1039-Collaboration/e1039-analysis.git
cd e1039-analysis/HitEmbedding
```

## Software structure

### Source codes

All source codes specific to this packages are stored in `src/`.
You have to compile them at least once by the following commands;
```
source setup.sh
# if you have DarkQuest's version of e1039-core:
source setup_mye1039.sh
```
and then:
```
cmake-this
make-this
```

`cmake-this` and `make-this` can be executed at any directory, since they don't depend on the current directory.
They have to be executed again after you modify the codes.

### Fun4All macros

Multiple Fun4All macros are available in this packages in order to carry out the hit embedding step-by-step.
Details of each macro are explained in `README.md` of each directory.
Probably you need not generate embedding events, since the E906 NIM3 events are already available and used by default.

1. [`macro_gen_signal/`](macro_gen_signal/):  Macro to generate signal events.
1. [`macro_gen_emb/`](macro_gen_emb/):  Macro to generate embedding events.
1. [`macro_gen_emb_e906/`](macro_gen_emb_e906/):  Macro to generate embedding events from E906 NIM3 events.
1. [`macro_embed/`](macro_embed/):  Macro to do the hit embedding.
1. [`macro_ana/`](macro_ana/):  Macro to analysis the hit-embedded events.

```
  [ Generate signal events ]                        ... `macro_gen_signal`
    |             |
    V             |
  [ Reconstruct ] |                                 ... `macro_gen_signal`
    |             |
    |             |  [ Generate embedding events ]  ... `macro_gen_emb` or `macro_gen_emb_e906`
    |             |    |
    |             V    V
    |           [ Embed hits ]                      ... `macro_embed`
    |             |
    |             V
    |           [ Reconstruct ]                     ... `macro_embed`
    |             |
    V             V
  [ Analyze ]   [ Analyze ]                         ... `macro_ana`
          |       |
          V       V
         [ Compare ]                                ... `macro_ana`
```

## To-do list

* Nothing planned.

## Contact

Feel free to send any questions/suggestions to Kenichi Nakano.

## Workflow for DQ:
- Generate signal events (e.g. particle gun, aprime signal). 
  - If possible generate them in batches of 5k events (easier to embed since e906 data is pre-processed in similar batches).
  - Execute `Fun4Sim.C` macro in `macro_gen_signal/` or execute `RecoE1039Sim.C` (with `run_sim.py` in `SimHits/`).

- Process events for embedding.
  - For normal e906 data we should not need to do this step again.
  - For a generated dataset (e.g. proton gun) we can use `macro_gen_emb/`.

- Embed events.
  - We can use the e906 2017 data that has been pre-processed.
  - Execute `Fun4Sim.C` macro in `macro_gen_embed/`.
