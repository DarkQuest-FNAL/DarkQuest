#export MY_E1039=/home/dsperka/DarkQuest-FNAL/core-inst/this-e1039.sh
export MY_E1039=/home/dsperka/wpmccormack/core-inst/this-e1039.sh
export DIR_TOP=$(dirname $(readlink -f $BASH_SOURCE))
echo $MY_E1039
echo $DIR_TOP
source $MY_E1039
export LD_LIBRARY_PATH=$DIR_TOP/install/lib/:$LD_LIBRARY_PATH
export DIR_CMANTILL=/cms/data/seaquest/users/cmantill/DarkQuest/e1039-analysis/SimHits
