# Macro to analyze data before/after hit embedding

This macro is to analyze the data before/after the hit embedding.
You can apply the identical analysis program to both the embedded data and the non-embedded data.
The comparison between them shows you the effect of embedded hits.
Thus you will carry out the analysis twice in two sub-directories; `not_embedded` and `embedded`.


## Analysis of non-embedded data

You first make a list of the input data, which are `DST.root` that you generetd with `macro_gen_signal`.
You execute the following commands;
```
cd not_embedded
../make-dst-list.sh ../../macro_gen_signal/data/jpsi_20211011
```
You see `list_dst.txt` and `info_dst.txt` are created.

You then produce a compact ROOT file (`ana_tree.root`) from the input data, by executing the following command;
```
root -b ../Fun4All.C
```
The contents of `ana_tree.root` are defined by `src/AnaEmbeddedData.cc`.
You will most likely modify it for your own study.

You last draw a set of plots using `ana_tree.root`, by executing the following command;
```
root -b ../AnaOneData.C
```
The plots will appear in `result/`.


## Analysis of embedded data

You repeat the same procedure to analyze the embedded data as well.
The input data this time are `DST.root` that you generated with `macro_embed`.
```
cd ../embedded
../make-dst-list.sh ../../macro_embed/data/20211012
root -b ../Fun4All.C
root -b ../AnaOneData.C
```


## Simultaneous analysis of non-embedded and embedded data

You can/should analyze both the non-embedded and embedded data simultaneously, 
for example to evaluate the relative reconstruction efficiency versus rate.
Such analysis can be done via `AnaBothData.C`, where the contents of 
the analysis are implemented in `src/AnaCleanAndMessyData.cc`;
```
cd ..
root -b AnaBothData.C
```
The plots will appear in `result/`.
