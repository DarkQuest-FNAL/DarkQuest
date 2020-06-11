#!/bin/bash

jobname=$1
do_sub=$2

njobs=$3
nevents=$4
ifile=$5

nmu=1
LIFE_TIME=long # short (3h), medium (8h) or long (23h)

if [ $do_sub == 1 ]; then
  echo "grid!"
else
  echo "local!"
fi

echo "njobs=$njobs"
echo "nevents=$nevents"

macros=$(dirname $(readlink -f $BASH_SOURCE))

sed "s/nevents=NAN/nevents=$nevents/"             $macros/gridrun.sh > $macros/gridrun_new.sh 
sed -i "s/nmu=NAN/nmu=$nmu/"                      $macros/gridrun_new.sh
sed -i "s/ifile=NAN/ifile=$ifile/"                $macros/gridrun_new.sh
chmod +x $macros/gridrun_new.sh

if [ $do_sub == 1 ]; then
work=/pnfs/e906/scratch/users/$USER/SimHits/$jobname
else
work=$macros/scratch/$jobname
fi

mkdir -p $work
chmod -R 01755 $work

cd $macros
cp /seaquest/users/cmantill/DarkQuest/lhe/displaced_Aprime_Electrons/$ifile.txt .
tar -czvf $work/input.tar.gz *.C *.opts $ifile.txt
rm $ifile.txt 
cd -

for (( id=1; id<=$njobs; id++ ))
do  
  mkdir -p $work/$id/log
  mkdir -p $work/$id/out
  chmod -R 01755 $work/$id

  rsync -av $macros/gridrun_new.sh $work/$id/gridrun_new.sh

  cmd="jobsub_submit"
  cmd="$cmd -g --OS=SL7 --use_gftp --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC,OFFSITE -e IFDHC_VERSION --expected-lifetime='$LIFE_TIME'"
  cmd="$cmd --mail_never"
  cmd="$cmd -L $work/$id/log/log.txt"
  cmd="$cmd -f $work/input.tar.gz"
  cmd="$cmd -d OUTPUT $work/$id/out"
  cmd="$cmd --append_condor_requirements='(TARGET.GLIDEIN_Site isnt \"UCSD\")'"
  cmd="$cmd file://`which $work/$id/gridrun_new.sh`"

  if [ $do_sub == 1 ]; then
    echo $cmd
    $cmd
  else
    mkdir -p $work/$id/input
    rsync -av $work/input.tar.gz $work/$id/input
    cd $work/$id/
    $work/$id/gridrun_new.sh | tee $work/$id/log/log.txt
    cd -
  fi
done 2>&1 | tee log_gridsub.txt

## When your job fails due to bad grid nodes,
## you can use the following option to exclude those nodes;
##   cmd="$cmd --append_condor_requirements='(TARGET.GLIDEIN_Site isnt \"UCSD\")'"
## Valid site names are listed here;
## https://cdcvs.fnal.gov/redmine/projects/fife/wiki/Information_about_job_submission_to_OSG_sites
## According to the Fermilab Service Desk, the "--blacklist" option has a known defect.
