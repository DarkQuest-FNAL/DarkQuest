export DIR_TOP=$(dirname $(readlink -f $BASH_SOURCE))
echo $DIR_TOP
source /seaquest/users/cmantill/core-inst/this-e1039.sh
export LD_LIBRARY_PATH=$DIR_TOP/install/lib/:$LD_LIBRARY_PATH
export DIR_CMANTILL=/seaquest/users/cmantill/DarkQuest/e1039-analysis/SimHits
