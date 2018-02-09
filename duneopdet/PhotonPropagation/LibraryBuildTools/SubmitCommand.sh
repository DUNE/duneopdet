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
checkVar=0
testVar=0
offsiteVar=0
memory=2500MB
expectedlifetime=8h
makeupJobs=0
scriptIn=OpticalLibraryBuild_Grid_dune.sh
outdir=/pnfs/dune/scratch/users/${USER}/OpticalLibraries/OpticalLib_dune10kt_v2_1x2x6
fclIn=dune10kt_v2_1x2x6_buildopticallibrary_grid.fcl
USER=${USER} #Set the user to the default USER from the environment unless over ridder
HELPFILE=SubmitCommand.hlp

##This block handles flags given to the program.
# Allowed flags are:
#    -t | --tar      : Pass a tarfile of a larsoft installation to be setup on the cluster.
#                        User full path to file.
#    -u | --user     : Over ride the user directory to write to on dCache *NOT RECOMENDED
#       | --test     : run the test fcl file and a single job with short run time instead of building a new library
#    -c | --check    : Preform a dry run, returning the jobsub command, but not actually running any grid jobs.
#    -l | --lifetime : The amount of time a job should be expected to run on the cluster. 
#    -s | --script   : The script to run on the grid for each job (By default OpticalLibraryBuild_Grid_dune.sh)
#    -f | --fcl      : The fcl file template to be given to $script. This will be used to build the individual fcl files for each grid job.
#    -o | --outdir   : The output directory for the results from the simulation. The default is the scratch space of the user who submitted the jobs.
#    -m | --memory   : the amount of memory to request from each node on the cluster. *NOT RECOMENDED 
#                       Allowed units are MB and GB
printf '////////////////////////////////////////////////////////////////////\n'
printf '////////       DUNE PhotonLibrary Build System /////////////////////\n'
while :; do
  case $1 in
    --debugSubmitCommand)
      printf "\nSetting script debugger. This is not normal run mode.\n"
      set -x
      set -v
      trap read debug
      ;;
    --help|-h)
      if [ -f $HELPFILE ]; then
        cat $HELPFILE
      else
        printf "Help File not found.\n" >&2
        exit 1
      fi
      exit 0
      ;;
    --script|-s)
      if [ "$2" ]; then
        scriptIn=$2
        printf "\nscriptIn set by user.\nscriptIn=$scriptIn\n"
        shift
      else
        printf 'ERROR: "--script" requires an input parameter to give the script to be executed on the grid nodes.\n' >&2
      fi
      ;;
    --script=?*)
      scriptIn=${1#*=}
      printf "\nscriptIn set by user.\nscriptIn=$scriptIn\n"
      ;;
    --fcl|-f)
      if [ "$2" ]; then
        fclIn=$2
        printf "\nInput fcl file set by user.\nfclIn=$fclIn\n"
        shift
      else
        printf 'ERROR: "--fcl" requires and input with the full path of the fcl file to be passed as the template for each grid job.\n' >&2
      fi
      ;;
    --fcl=?*)
      fclIn=${1#*=}
      printf "\nInput fcl file set by user.\nfclIn=$fclIn\n"
      ;;
    --outdir|-o)
      if [ "$2" ]; then
        outdir=$2
        outdir=$(sed 's/\/ *$//' <<<$outdir)
        printf "\noutput directory set by user.\noutput directory will be $outdir\n"
        shift
      else
        printf 'ERROR: "--outdir" requires an input with the full path of the directory where library generation results should be written.\n' >&2
      fi
      ;;
    --outdir=?*)
      outdir=${1#*=}
      outdir=$(sed 's/\/ *$//' <<<$outdir)
        printf "\noutput directory set by user.\noutput directory will be $outdir\n"
      ;;
    --memory|-m)
      if [ "$2" ]; then
        memory=$2
        printf "\nCluster memory requirement set by user.\nmemory request will be $memory\n"
        shift
      else
        printf 'ERROR: "--memory" requires an input value to pass on to the cluster.\n' >&2
        exit 10
      fi
      ;;
    --memory=?*)
      memory=${1#*=}
        printf "\nCluster memory requirement set by user.\nmemory request will be $memory\n"
      ;;
    --makeup|-n)
      if [ "$2" ]; then
        makeupJobs=$2
        printf "\nNumber Of Jobs required for Makeup Jobs set. Your OpticalLibraryBuild_Grid_dune.sh.\n If your OpticalLibraryBuild_Grid_dune.sh does not contain the correct makeup list, this step will not behave as expecte.\n"
        shift
      else
        printf 'ERROR: "--makeup" requires the number of makeup jobs to process.\n'
        exit 10
      fi
      ;;
    --makeup=?*)
        makeupJobs=${1#*=}
        printf "\nNumber Of Jobs required for Makeup Jobs set. Your OpticalLibraryBuild_Grid_dune.sh.\n If your OpticalLibraryBuild_Grid_dune.sh does not contain the correct makeup list, this step will not behave as expecte.\n"
      ;;
    --tar|-t)
      if [ "$2" ]; then
        tarfile=$2
        printf "\nInput tar file specified.\nLArSoft/dunetpc will be setup from $tarfile\n"
        shift
      else
        printf 'ERROR: "--tar" requires a path to a tar file of a larsoft installation.\n' >&2
        exit 10
      fi
      ;;
    --tar=?*)
      tarfile=${1#*=}
        printf "\nInput tar file specified.\nLArSoft/dunetpc will be setup from $tarfile\n"
      ;;
    --check|-c)
      checkVar=1
      printf "\nSetting check mode ON.\n"
      ;;
    --test)
      testVar=1
      printf "\nSetting test mode ON.\n"
      ;;
    --offsite)
      offsiteVar=1
      printf "\nAllow jobs to go offsite (e.g. the OSG).\n"
      ;;
    --lifetime|-l)
      if [ "$2" ]; then
        expectedlifetime=$2
        printf "\nCluster upper limit on runtime set by user\nThe requested runtime will be $expectedlifetime\n"
        shift
      else
        printf 'ERROR: "--lifetime" requires a parameter telling the cluster what the upper bound of each jobs runtime should be.\n' >&2
      fi
      ;;
    --lifetime=?*)
      expectedlifetime=${1#*=}
        printf "\nCluster upper limit on runtime set by user\nThe requested runtime will be $expectedlifetime\n"
      ;;
    --user|-u)
      if [ "$2" ]; then
        USER=$2
        printf "\nUSER for outputs over-ridden by user. CAUTION: This will likely not work. Only do this if you know exaclty why you are doing so.\nUser is set to $USER\n"
        shift
      else
        printf 'ERROR: "--user" requires a username to use for the dCache directory.\n' >&2
        exit 10
      fi
      ;;
    --user=?*)
      USER=${1#*=}
        printf "\nUSER for outputs over-ridden by user. CAUTION: This will likely not work. Only do this if you know exaclty why you are doing so.\nUser is set to $USER\n"
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



if [ ! -d $outdir/root ]; then
  mkdir -p $outdir/root
  mkdir -p $outdir/fcl
  mkdir -p $outdir/log
fi

printf "\nfclIn=$fclIn\n"
printf "basename fclIn=$(basename $fclIn)\n"

fcl="$outdir/$(basename $fclIn)"

printf "fcl set.\nfcl=$fcl\n"

printf "\nscriptIn=$scriptIn\n"
printf "basename scriptIn=$(basename $scriptIn)\n"

script="$outdir/`basename $scriptIn`"

printf "script set.\nscript=$script\n"

if [ -e $fcl ]; then
  printf "\n$fcl already exists. Removing old file and replacing with new.\n"
  rm -f $fcl
fi
printf "\nPreparing fcl for transfer to the grid.\ncp $fclIn $fcl\n"
if [ -e $fclIn ]; then
  cp $fclIn $fcl
else
  printf "\nExiting with error. Source file for fcl not found. \nPlease make sure the fcl \n$fclIn \nexists.\n"
  exit 10
fi

if [ -e $script ]; then
  printf "\n$script already exists. Removing old file and replacing with new.\n"
  rm -f $script
fi
printf "\nPreparing script for transfer to the grid.\ncp $scriptIn $script\n"
if [ -e $scriptIn ]; then
  cp $scriptIn $script
else
  printf "\nExiting with error. Source file for Script not found. \nPlease make sure the script \n$scriptIn \nexists.\n"
  exit 10
fi

environmentVars="-e IFDH_CP_MAXRETRIES=5"
usage="DEDICATED,OPPORTUNISTIC"
if [ $offsiteVar -ne 0 ]; then
    usage="DEDICATED,OPPORTUNISTIC,OFFSITE"
fi
clientargs="--resource-provides=usage_model=$usage --OS=SL6 --group=dune -f $fcl --role=Analysis --memory=$memory "
if [ x$tarfile != x ]; then
  printf "\nUsing tarball. Not setting LArSoft environment variables\n"
  larsoft=
  clientargs="${clientargs} --tar_file_name=dropbox://${tarfile} "
else
  larsoft="${environmentVars} -e mrb_top=$MRB_TOP -e mrb_project=dunetpc -e mrb_version=$MRB_PROJECT_VERSION -e mrb_quals=$MRB_QUALS "
fi

toolsargs="-q -g --opportunistic --OS=SL6 "
fileargs="-dROOT $outdir/root -dFCL $outdir/fcl -dLOG $outdir/log "

#Test job vs real job
if [ $testVar -ne 0 ]; then #TEST VAR IS SET. Run the test job
  printf "\n!!TEST JOB SET!!\nOnly a single job is being sent to the grid. This is for testing and debugging purposes, and will not build a complete library.\n"
  #Test job 1 - jobsub_client
  njobs=300000 #This is picked to select 10 voxels for 100x100x300 bins with 10 photons each. 
  nphotons=10
  clientargs="$clientargs --expected-lifetime=$expectedlifetime "
  thisjob="-Q -N 1 file://$script $njobs $nphotons $(basename $fcl)"
else  
  printf "Building Library\n"
  #Real job - jobsub_client
  njobs=6000
  nphotons=50000
  clientargs="$clientargs --expected-lifetime=$expectedlifetime "
  #  thisjob="-N $njobs file://$script $njobs $nphotons"
  if [ 0 -ne $makeupJobs ]; then
    echo "thisjob=\"-N $makeupJobs file://$script $njobs $nphotons $(basename $fcl) true\""
    thisjob="-N $makeupJobs file://$script $njobs $nphotons $(basename $fcl) true"
  else
    thisjob="-N $njobs file://$script $njobs $nphotons $(basename $fcl)"
  fi
fi

if [ x$tarfile != x ]; then
  printf "\n\njobsub_submit $environmentVars $clientargs $fileargs $thisjob \n\n\n"
  if [ $checkVar -ne 0 ]; then
    printf "CHECK Mode is set. The jobsub command will be printed, but will not be executed. Please check the command and run again without check mode. If you are trying to submit test jobs instead, the correct flag is -s or --test.\n"
  else
    jobsub_submit $environmentVars $clientargs $fileargs $thisjob 
  fi
  ret=$?
  printf "\nExiting with status $ret\n"
  exit $ret
else
  printf "jobsub_submit $environmentVars $larsoft $clientargs $fileargs $thisjob\n"
  if [ $checkVar -ne 0 ]; then
    printf "\n\nCHECK Mode is set. The jobsub command will be printed, but will not be executed. Please check the command and run again without check mode. If you are trying to submit test jobs instead, the correct flag is -s or --test.\n\n\n"
  else
    jobsub_submit $environmentVars $larsoft $clientargs $fileargs $thisjob 
  fi
  ret=$?
  printf "\nExiting with status $ret\n"
  exit $ret
fi


