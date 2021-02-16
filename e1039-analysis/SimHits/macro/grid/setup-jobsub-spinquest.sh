#!/bin/bash 
# A script to set up the use of the grid computer in E1039.
# It was extracted from /e906/app/users/yuhw/setup.sh.

GROUP=spinquest #legacy for jobsub_tools, who knows what else
test $(basename $BASH_SOURCE) = 'setup-jobsub-seaquest.sh' -o "X$1" = 'Xseaquest' && GROUP=seaquest
export GROUP

export JOBSUB_GROUP=$GROUP  #for jobsub_client
#export   KRB5CCNAME=/e906/app/users/$USER/.krbcc/my_cert

source /grid/fermiapp/products/common/etc/setups.sh
setup ifdhc
setup cpn
setup jobsub_client
#setup git

#test -e $KRB5CCNAME || mkdir -p $KRB5CCNAME

alias       jobsub_q_mine='jobsub_q       --user=$USER'
alias      jobsub_rm_mine='jobsub_rm      --user=$USER'
alias jobsub_history_mine='jobsub_history --user=$USER'

function jobsub_fetchlog_and_extract {
    local -r JOBID=$1
    mkdir -p $JOBID
    cd       $JOBID
    jobsub_fetchlog --job=$JOBID
    test -f $JOBID.tgz && tar xvzf $JOBID.tgz
}
