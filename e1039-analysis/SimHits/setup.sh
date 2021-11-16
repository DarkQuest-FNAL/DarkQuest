export DIR_TOP=$(dirname $(readlink -f $BASH_SOURCE))
echo $DIR_TOP
source /e906/app/software/osg/software/e1039/this-e1039.sh
export LD_LIBRARY_PATH=$DIR_TOP/install/lib/:$LD_LIBRARY_PATH
export ROOT_INCLUDE_PATH=$DIR_TOP/install/include:$ROOT_INCLUDE_PATH
export DIR_CMANTILL=/seaquest/users/cmantill/DarkQuest/e1039-analysis/SimHits
