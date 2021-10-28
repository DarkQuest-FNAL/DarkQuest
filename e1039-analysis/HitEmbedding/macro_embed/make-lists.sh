#!/bin/bash
## Script to make lists of signal DST files and embedding data files.
## Modify 'DIR_SIG' and 'DIR_EMB' and execute './make-lists.sh'.

## You can put multiple directories separated by space 
## and also use the shell wild card (*).
DIR_SIG="../macro_gen_signal/data/jpsi_20211011"
DIR_EMB="/pnfs/e1039/persistent/users/kenichi/data_emb_e906"

##
## Functions
##
function GetTreeEntries {
    local -r FILE=$1
    local -r TREE=$2
    local -r COMM="TFile* file = new TFile(\"$FILE\");
                   TTree* tree = (TTree*)file->Get(\"$TREE\");
                   cout << tree->GetEntries() << endl;"
    echo "$COMM" | root.exe -l -b 2>/dev/null
}

function MakeOneList {
    local -r FN_LIST=$1
    local -r DIR_DAT=$2
    local -r BN_DAT=$3
    local -r TREE_NAME=$4
    echo "List: $FN_LIST"
    echo "  Directory = $DIR_DAT"
    echo "  File name = $BN_DAT"
    echo "  Tree name = $TREE_NAME"

    local -r LIST_DAT=( $(find $DIR_DAT -name $BN_DAT) )
    local -r N_DAT=${#LIST_DAT[*]}
    echo "  N of files = $N_DAT"
    if [ $N_DAT -eq 0 ] ; then
    	echo "  Do nothing since no file is available."
	return
    fi

    local -r N_EVT=$(GetTreeEntries ${LIST_DAT[0]} $TREE_NAME)
    echo "  N of events in the 1st file = $N_EVT"
    echo "  N of total events ~= $(( N_DAT * N_EVT ))"
    
    echo "${LIST_DAT[*]}" | xargs -r readlink -m >$FN_LIST

    #local N_EVT_TOT=0
    #for FN_DAT in ${LIST_DAT[*]} ; do
    #	N_EVT=$(GetTreeEntries ${LIST_DAT[0]} $TREE_NAME)
    #	(( N_EVT_TOT += N_EVT ))
    #	echo "  $N_EVT $FN_DAT"
    #done
}

##
## Main
##
if ! which root.exe &>/dev/null ; then
    echo "Command 'root' not found.  Forget 'source setup.sh'?  Abort."
    exit
fi

{
echo "Time: $(date '+%F %H:%M:%S')"
echo
MakeOneList 'list_sig_file.txt' "$DIR_SIG" 'DST.root' 'T'
echo
MakeOneList 'list_emb_file.txt' "$DIR_EMB" 'embedding_data.root' 'tree'
} |& tee info_list.txt
