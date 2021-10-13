export MY_E1039=/seaquest/users/$USER/core-inst/this-e1039.sh
export DIR_TOP=$(dirname $(readlink -f $BASH_SOURCE))
echo $MY_E1039
echo $DIR_TOP
source $MY_E1039
export LD_LIBRARY_PATH=$DIR_TOP/install/lib/:$LD_LIBRARY_PATH
export DIR_CMANTILL=/seaquest/users/cmantill/DarkQuest/e1039-analysis/SimHits
