#!/bin/bash
#
# A script to run the optical library building job
#
#
# To run this job:
#
# jobsub -N [NoOfJobs] -M -dROOT /out/dir/root -dFCL /out/dir/fcl -dLOG /out/dir/log OpticalLibraryBuild_Grid_dune.sh
#
# You will get outputs in the area specified by the "outstage" variable 
# which is specified below.
#
# The form of the output is one file for each few voxels. These then need 
# stitching together, which is done after all jobs are done, with a
# dedicated stitching script.
#

#
#
##
#if ! echo $missing | grep -w "$PROCESS" > /dev/null; then
#  echo "Exiting since $PROCESS is already complete."
#  exit 0
#fi


#
# Set up our environment
#

njobs=$1
label=${CLUSTER}_$(printf '%04d' $PROCESS)
## This sets all the needed FW and SRT and LD_LIBRARY_PATH envt variables. 
## Then we cd back to our TMP area. EC, 23-Nov-2010.

# In each voxel, run this many photons:
NPhotonsPerVoxel=$2
fclIn=$3


umask 0002

export GROUP=dune
export HOME=$CONDOR_DIR_ROOT
export CONDOR_SCRATCH_DIR=$TEMP

#LOG=./pd_library_gen_${label}.log
#FCL=./pd_library_gen_${label}.fcl
LOG=${CONDOR_DIR_LOG}/pd_library_gen_${label}.log
FCL=${CONDOR_DIR_FCL}/pd_library_gen_${label}.fcl
RTF=${CONDOR_DIR_ROOT}/pd_library_gen_${label}.root

#
# Prepare to run
#
cd $CONDOR_DIR_ROOT
touch ${LOG}

echo "Writing log to ${LOG}"        1>>${LOG} 2>&1 
echo "OpticalLibraryBuild_Grid_dune.sh command is:"     1>>${LOG} 2>&1 
echo "OpticalLibraryBuild_Grid_dune.sh $1 $2 $3"       1>>${LOG} 2>&1 
echo "njobs=$njobs"          1>>${LOG} 2>&1 
echo "label=$label"         1>>${LOG} 2>&1 
echo "NPhotonsPerVoxel=$NPhotonsPerVoxel"       1>>${LOG} 2>&1    
echo "log=$LOG"          1>>${LOG} 2>&1 
echo "FLC=$FCL"          1>>${LOG} 2>&1 
echo "RTF=$RTF"         1>>${LOG} 2>&1 


date 1>> ${LOG} 2>&1



#
# Library building parameters
#

# Copy fcl file and configure for this PROCESS 
echo "Create this job's fhicl file" 1>>${LOG} 2>&1
echo "checking for FHICL location." 1>>${LOG} 2>&1

ls ${CONDOR_DIR_INPUT}/*.fcl  1>> ${LOG} 2>&1
ls ${CONDOR_DIR_FCL}/$fclIn  1>> ${LOG} 2>&1
echo "Listing condor directories and contents for debugging." 1>>${LOG} 2>&1
echo "CONDOR_DIR_LOG=${CONDOR_DIR_LOG}" 2>>${LOG} 2>&1
ls ${CONDOR_DIR_LOG} 1>>${LOG} 2>&1
echo "CONDOR_DIR_FCL=${CONDOR_DIR_FCL}" 2>>${LOG} 2>&1
ls ${CONDOR_DIR_FCL} 1>>${LOG} 2>&1
echo "CONDOR_DIR_ROOT=${CONDOR_DIR_ROOT}" 2>>${LOG} 2>&1
ls ${CONDOR_DIR_ROOT} 1>>${LOG} 2>&1
echo "CONDOR_DIR_INPUT=${CONDOR_DIR_INPUT}" 2>>${LOG} 2>&1
ls ${CONDOR_DIR_INPUT} 1>>${LOG} 2>&1


cp -v ${CONDOR_DIR_INPUT}/$fclIn $FCL 1>> ${LOG} 2>&1
#cp -v ${CONDOR_DIR_FCL}/$fclIn $FCL 1>> ${LOG} 2>&1

NX=`awk '/NX/{ print $2 }' $FCL`
NY=`awk '/NY/{ print $2 }' $FCL`
NZ=`awk '/NZ/{ print $2 }' $FCL`
echo "NX=$NX"
echo "NY=$NY"
echo "NZ=$NZ"

# Total number of voxels
#NTopVoxel=216000
NTopVoxel=`echo "$NX*$NY*$NZ" | bc`

echo "Voxels: $NTopVoxel = $NX * $NY * $NZ" 1>> ${LOG} 2>&1





# In each grid job, do this many voxels:
NVoxelsPerJob=`echo "$NTopVoxel/$njobs" | bc`
echo "NVoxelsPerJob=$NVoxelsPerJob" 1>>${LOG} 2>&1

# This works out which voxels this job should focus on: 
FirstVoxel=`echo "($NVoxelsPerJob * $PROCESS ) % $NTopVoxel" | bc` 
LastVoxel=`echo "(($NVoxelsPerJob * $PROCESS ) + $NVoxelsPerJob - 1 ) % $NTopVoxel" | bc`

# Construct the run-time configuration file for this job;
# Make the random number seeds a function of the PROCESS number.
generatorSeed=$(( $PROCESS *23 + 31))
g4Seed=$(( $PROCESS *41 + 37))


echo "physics.producers.generator.FirstVoxel: $FirstVoxel" >> $FCL
echo "physics.producers.generator.LastVoxel: $LastVoxel"   >> $FCL
echo "physics.producers.generator.N: $NPhotonsPerVoxel"    >> $FCL





# And then tell the user about it:
echo "This job will run from voxel $FirstVoxel to $LastVoxel, generating $NPhotonsPerVoxel in each" 1>> ${LOG} 2>&1
echo "CLUSTER:    " $CLUSTER 1>> ${LOG} 2>&1
echo "PROCESS:    " $PROCESS 1>> ${LOG} 2>&1
echo "PWD:        " $PWD     1>> ${LOG} 2>&1


# No need to set random seeds - generated from machine state
#echo "physics.producers.generator.RandomSeed: $generatorSeed">> $FCL
#echo "physics.producers.largeant.RandomSeed: $g4Seed">> $FCL


echo "> Setup $GROUP environment"                                                              1>> ${LOG} 2>&1
echo "> source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh"                    1>> ${LOG} 2>&1
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh                             1>> ${LOG} 2>&1

echo "> INPUT_TAR_FILE=$INPUT_TAR_FILE"                                                        1>> ${LOG} 2>&1
echo ">x\$INPUT_TAR_FILE={x$INPUT_TAR_FILE}"                                                      1>> ${LOG} 2>&1
echo "Bool {x$INPUT_TAR_FILE != x}"                                                              1>>${LOG} 2>&1

if [ x$INPUT_TAR_FILE != x ]; then #This is how we handle users trying to pass their local installation as a tarball.
  echo "CONDOR_SCRATCH_DIR=$CONDOR_SCRATCH_DIR"   1>>${LOG} 2>&1
  echo "TMP=$TMP"  1>>${LOG} 2>&1
  echo "TEMP=$TEMP"  1>>${LOG} 2>&1
  mkdir $CONDOR_SCRATCH_DIR/local 1>>${LOG} 2>&1
  cd $CONDOR_SCRATCH_DIR/local  1>>${LOG} 2>&1
  pwd                         1>>${LOG} 2>&1
  echo "Extracting TAR"      1>>${LOG} 2>&1
  cd ${TMP}/local
  tar -xf $INPUT_TAR_FILE   1>>${LOG} 2>&1

  echo "Extracted TAR"     1>>${LOG} 2>&1
  echo "Files extracted are:"     1>>${LOG} 2>&1
  ls                       1>>${LOG} 2>&1
  echo "Initializing localProducts from tarball ${INPUT_TAR_FILE}."     1>>${LOG} 2>&1
  echo 'sed "s@setenv MRB_INSTALL.*@setenv MRB_INSTALL ${TMP}/local@" $TMP/local/setup | \'
  echo 'sed "s@setenv MRB_TOP.*@setenv MRB_TOP ${TMP}@" > $TMP/local/setup.local | \'
  sed "s@setenv MRB_INSTALL.*@setenv MRB_INSTALL ${TMP}/local@" $TMP/local/setup | \
    sed "s@setenv MRB_TOP.*@setenv MRB_TOP ${TMP}@" > $TMP/local/setup.local 

  echo "Setting up products"     1>>${LOG} 2>&1
  #. ${TMP}/local/setup             1>>${LOG} 2>&1
  . ${TMP}/local/setup.local     1>>${LOG} 2>&1
  echo "MRB_PROJECT=$MRB_PROJECT"                  1>>${LOG} 2>&1
  echo "MRB_PROJECT_VERSION=$MRB_PROJECT_VERSION"  1>>${LOG} 2>&1
  echo "MRB_QUALS=$MRB_QUALS"                      1>>${LOG} 2>&1
  echo "MRB_TOP=$MRB_TOP"                          1>>${LOG} 2>&1
  echo "MRB_SOURCE=$MRB_SOURCE"                    1>>${LOG} 2>&1
  echo "MRB_BUILDDIR=$MRB_BUILDDIR"                1>>${LOG} 2>&1
  echo "MRB_INSTALL=$MRB_INSTALL"                  1>>${LOG} 2>&1

  echo "mrbslp next"                               1>>${LOG} 2>&1
  mrbslp                                           1>>${LOG} 2>&1
  echo "Setup tarbal done"                         1>>${LOG} 2>&1
else  #If we are not using a local install, setup becomes much easier.
  echo "> ups list -aK+ $mrb_project"                                                            1>> ${LOG} 2>&1
  ups list -aK+ $mrb_project                                                                     1>> ${LOG} 2>&1
  echo "> setup $mrb_project $mrb_version -q$mrb_quals"                                          1>> ${LOG} 2>&1
  setup $mrb_project $mrb_version -q $mrb_quals                                                  1>> ${LOG} 2>&1
fi

echo "Setting up environment log" 1>>${LOG} 2>&1
ENVLOG=$CONDOR_DIR_LOG/environment_${label}.log
touch $ENVLOG
uname -a 1>> $LOG     2>&1
uname -a 1>> $ENVLOG  2>&1
cat /etc/redhat-release 1>> $LOG     2>&1
cat /etc/redhat-release 1>> $ENVLOG  2>&1
env >> $ENVLOG


# Run the job
echo ""                                                                 1>> ${LOG} 2>&1
echo "***** Starting job"                                               1>> ${LOG} 2>&1
echo "lar -c $FCL -n $NVoxelsPerJob -T ${RTF}" 1>> ${LOG} 2>&1
lar -c $FCL -n $NVoxelsPerJob -T ${RTF}        1>> ${LOG} 2>&1
ret=$?
echo   "***** Job completed ($ret)"                                     1>> ${LOG} 2>&1
echo                                                                    1>> ${LOG} 2>&1
date                                                                    1>> ${LOG} 2>&1

exit $ret
