############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
## SumbitCommand.sh is used for submitting DUNE Photon Library Generation ##
## using the FNAL Grid                                                    ##
## Author: Alex Himmel (ahimmel@fnal.gov)                                 ##
## Updated by Jason Stock (jason.stock@mines.sdsmt.edu) 2017-09-11        ##
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################

#!/bin/bash

tarfile=
debugVar=0
memory=2500MB
USER=${USER} #Set the user to the default USER from the environment unless over ridder

##This block handles flags given to the program.
# Allowed flags are:
#    -t | --tar    : Pass a tarfile of a larsoft installation to be setup on the cluster.
#                     User full path to file.
#    -u | --user   : Over ride the user directory to write to on dCache *NOT RECOMENDED
#    -d | --test   : run the test fcl file and a single job with short run time instead of building a new library
#    -m | --memory : the amount of memory to request from each node on the cluster. *NOT RECOMENDED 
#                     Allowed units are MB and GB
while :; do
  case $1 in
    --memory|-m)
      if [ "$2" ]; then
        memory=$2
      else
        printf 'ERROR: "--memory" requires an input value to pass on to the cluster.\n' >&2
        exit 10
      fi
      ;;
    --memory=?*)
      memory=${1#*=}
      ;;
    --tar|-t)
      if [ "$2" ]; then
        tarfile=$2
        shift
      else
        printf 'ERROR: "--tar" requires a path to a tar file of a larsoft installation.\n' >&2
        exit 10
      fi
      ;;
    --test|-d)
      debugVar=1
      ;;
    --tar=?*)
      tarfile=${1#*=}
      ;;
    --user|-u)
      if [ "$2" ]; then
        USER=$2
        shift
      else
        printf 'ERROR: "--user" requires a username to use for the dCache directory.\n' >&2
        exit 10
      fi
      ;;
    --user=?*)
      USER=${1#*=}
      ;;
    --)
      shift
      break
      ;;
    -?*)
      printf 'ERROR: Uknown option\n'
      exit 10
      ;;
    *)
      break
  esac
  shift
done


script=OpticalLibraryBuild_Grid_dune.sh
outdir=/pnfs/dune/scratch/users/${USER}/OpticalLibraries/OpticalLib_dune10kt_v2_1x2x6
fcl=$outdir/dune10kt_v2_1x2x6_buildopticallibrary_grid.fcl

if [ ! -d $outdir/root ]; then
    mkdir -p $outdir/root
    mkdir -p $outdir/fcl
    mkdir -p $outdir/log
fi

if [ ! -e $fcl ]; then
    cp `basename $fcl` $fcl
fi

# 57600 seconds = 16 hours, but it will not be a sharp cut-off
environmentVars="-e IFDH_CP_MAXRETRIES=5"
clientargs="--resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --OS=SL6 --group=dune -f $fcl --role=Analysis --memory=$memory "
if [ x$tarfile != x ]; then
  echo "Using tarball. Not setting LArSoft environment variables"
  larsoft=
  clientargs="${clientargs} --tar_file_name=dropbox://${tarfile} "
else
  larsoft="${environmentVars} -e mrb_top=$MRB_TOP -e mrb_project=dunetpc -e mrb_version=$MRB_PROJECT_VERSION -e mrb_quals=$MRB_QUALS "
fi

toolsargs="-q -g --opportunistic --OS=SL6 "
fileargs="-dROOT $outdir/root -dFCL $outdir/fcl -dLOG $outdir/log "

#Test job vs real job
if [ $debugVar -ne 0 ]; then
  echo "test Job"
  #Test job 1 - jobsub_client
  njobs=7200
  nphotons=10
  clientargs="$clientargs --expected-lifetime=600 "
  thisjob="-Q -N 1 file://$PWD/$script $njobs $nphotons"
else  #Debug var is set. Run the test job
  echo "Building Library"
  #Real job - jobsub_client
  njobs=6000
  nphotons=50000
  clientargs="$clientargs --expected-lifetime=10h "
  thisjob="-N $njobs file://$PWD/$script $njobs $nphotons"
fi

if [ x$tarfile != x ]; then
  echo "jobsub_submit $environmentVars $clientargs $fileargs $thisjob "
  jobsub_submit $environmentVars $clientargs $fileargs $thisjob 
  ret=$?
  printf "Exiting with status $ret\n"
  exit $ret
else
  echo "jobsub_submit $environmentVars $larsoft $clientargs $fileargs $thisjob"
  jobsub_submit $environmentVars $larsoft $clientargs $fileargs $thisjob 
  ret=$?
  printf "Exiting with status $ret\n"
  exit $ret
fi


