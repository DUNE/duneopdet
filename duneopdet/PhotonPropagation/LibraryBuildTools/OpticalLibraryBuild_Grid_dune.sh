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
# Set up our environment
#

njobs=$1
label=${CLUSTER}_$(printf '%04d' $PROCESS)

umask 0002

export GROUP=dune
export HOME=$CONDOR_DIR_ROOT

LOG=${CONDOR_DIR_LOG}/pd_library_gen_${label}.log
FCL=${CONDOR_DIR_FCL}/pd_library_gen_${label}.fcl


#
# Library building parameters
#


# Copy fcl file and configure for this PROCESS 
echo "Create this job's fhicl file" 1>> ${LOG} 2>&1
mv -v ${CONDOR_DIR_INPUT}/*.fcl $FCL 1>> ${LOG} 2>&1

NX=`awk '/NX/{ print $2 }' $FCL`
NY=`awk '/NY/{ print $2 }' $FCL`
NZ=`awk '/NZ/{ print $2 }' $FCL`

# Total number of voxels
#NTopVoxel=216000
NTopVoxel=`echo "$NX*$NY*$NZ" | bc`

echo "Voxels: $NTopVoxel = $NX * $NY * $NZ" 1>> ${LOG} 2>&1


# In each voxel, run this many photons:
NPhotonsPerVoxel=$2



# In each grid job, do this many voxels:
NVoxelsPerJob=`echo "$NTopVoxel/$njobs" | bc`

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




#
# Prepare to run
#
cd $CONDOR_DIR_ROOT
echo "Writing log to ${LOG}"
touch ${LOG}
date 1>> ${LOG} 2>&1

# And then tell the user about it:
echo "This job will run from voxel $FirstVoxel to $LastVoxel, generating $NPhotonsPerVoxel in each" 1>> ${LOG} 2>&1
echo "CLUSTER:    " $CLUSTER 1>> ${LOG} 2>&1
echo "PROCESS:    " $PROCESS 1>> ${LOG} 2>&1
echo "PWD:        " $PWD     1>> ${LOG} 2>&1


# No need to set random seeds - generated from machine state
#echo "physics.producers.generator.RandomSeed: $generatorSeed">> $FCL
#echo "physics.producers.largeant.RandomSeed: $g4Seed">> $FCL



## This sets all the needed FW and SRT and LD_LIBRARY_PATH envt variables. 
## Then we cd back to our TMP area. EC, 23-Nov-2010.

mrb_top=$3
mrb_project=$4
mrb_version=$5
mrb_quals=$6

echo "> Setup $GROUP environment"                                                              1>> ${LOG} 2>&1
echo "> source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh"                    1>> ${LOG} 2>&1
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh                             1>> ${LOG} 2>&1

#echo "> source ${mrb_top}/localProducts_${mrb_project}_${mrb_version}_${mrb_quals/:/_}/setup"  1>> ${LOG} 2>&1
#source ${mrb_top}/localProducts_${mrb_project}_${mrb_version}_${mrb_quals/:/_}/setup           1>> ${LOG} 2>&1
#echo "> mrbslp"                                                                                1>> ${LOG} 2>&1
#source $MRB_DIR/bin/setup_local_products                                                       1>> ${LOG} 2>&1
echo "> ups list -aK+ $mrb_project"                                                            1>> ${LOG} 2>&1
ups list -aK+ $mrb_project                                                                     1>> ${LOG} 2>&1
echo "> setup $mrb_project $mrb_version -q$mrb_quals"                                          1>> ${LOG} 2>&1
setup $mrb_project $mrb_version -q$mrb_quals                                                   1>> ${LOG} 2>&1

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
echo "lar -c $FCL -n $NVoxelsPerJob -T PhotonLibraryFile_${label}.root" 1>> ${LOG} 2>&1
lar -c $FCL -n $NVoxelsPerJob -T PhotonLibraryFile_${label}.root        1>> ${LOG} 2>&1
ret=$?
echo   "***** Job completed ($ret)"                                     1>> ${LOG} 2>&1
echo                                                                    1>> ${LOG} 2>&1
date                                                                    1>> ${LOG} 2>&1

exit $ret

