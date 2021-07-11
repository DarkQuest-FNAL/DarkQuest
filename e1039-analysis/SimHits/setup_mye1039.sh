export MY_E1039=/data/t3home000/dhoang/DQ/e1039-core/DarkQuest/e1039-analysis/SimHits/../../../../core-inst/this-e1039.sh
export DIR_TOP=$(dirname $(readlink -f $BASH_SOURCE))
echo $DIR_TOP
source $MY_E1039
export LD_LIBRARY_PATH=$DIR_TOP/install/lib/:$LD_LIBRARY_PATH
export DIR_CMANTILL=/data/t3home000/pharris/DQ/data/
