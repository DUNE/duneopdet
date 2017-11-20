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


#Define functions

SetProcessNumber(){
  #All of these outputs don't make much sense on the node, since they can't be written to log.
  local  ourLog=$1
  local myProcess=$2
  printf "Beginning with PROCESS# $PROCESS\n" 1>>${ourLog} 2>&1
  printf "Correcting process number.\n" 1>>${ourLog} 2>&1
  #using an indexed bash array to make up jobs. 1>>${ourLog} 2>&1
  local missingJobs=(1082 1086 1087 1088 1530 1532 1535 1542 1543 1545 1547 1549 1552 1554 1555 1562 1564 1572 1594 1596 1601 1604 1606 1623 1624 1625 1626 1627 1628 1635 1636 1643 1644 1645 1646 1651 1656 1657 1663 1665 1666 1668 1670 1671 1672 1674 1675 1676 1677 1683 1684 1703 1704 1705 1717 1722 1741 1742 1743 1755 1757 1761 1767 1773 1774 1775 1776 1781 1787 1793 1796 1803 1808 1822 1823 1828 1829 1832 1861 1864 1865 1869 1882 1886 1894 1895 1901 1906 1915 1916 1918 1922 1925 1932 1934 1936 1940 1944 1947 1948 1950 1955 1960 1964 1965 1972 1978 1980 1983 1985 1987 1990 1994 1996 1998 2000 2002 2003 2006 2011 2012 2013 2014 2017 2021 2025 2026 2027 2030 2038 2040 2041 2044 2050 2051 2054 2056 2067 2074 2076 2077 2078 2081 2083 2085 2086 2089 2090 2091 2092 2093 2095 2096 2097 2098 2100 2101 2102 2103 2106 2112 2114 2115 2120 2121 2122 2123 2124 2125 2126 2127 2129 2130 2135 2137 2138 2140 2141 2143 2154 2157 2160 2161 2166 2175 2178 2180 2182 2186 2192 2196 2198 2200 2201 2202 2204 2210 2212 2213 2214 2218 2219 2229 2230 2236 2238 2239 2240 2246 2255 2265 2278 2279 2283 2296 2302 2304 2305 2309 2312 2317 2318 2319 2320 2325 2329 2335 2337 2338 2340 2341 2342 2345 2347 2349 2355 2357 2358 2361 2362 2380 2392 2393 2394 2396 2399 2410 2418 2420 2438 2440 2454 2456 2458 2460 2461 2473 2474 2475 2476 2477 2478 2479 2480 2481 2482 2483 2484 2485 2486 2494 2496 2500 2501 2518 2538 2540 2575 2578 2580 2587 2617 2618 2620 2621 2623 2638 2639 2642 2658 2660 2671 2676 2677 2678 2680 2681 2682 2689 2690 2693 2696 2697 2698 2699 2700 2701 2702 2703 2718 2719 2720 2731 2737 2738 2739 2756 2758 2759 2778 2779 2780 2799 2818 2819 2820 2839 2858 2859 2879 2899) 1>>${ourLog} 2>&1
  printf "Missing Job Numbers identified as \n" 1>>${ourLog} 2>&1
  for jobNum in missingJobs; do
    printf "$jobNum\n" 1>>${ourLog} 2>&1
  done
  printf "Old Process Number: $myProcess\n" 1>>${ourLog} 2>&1
  local myProcess=${missingJobs[$PROCESS]} 
  printf "New Process Number: $myProcess (for making the right voxels)\n"  1>>${ourLog} 2>&1
  echo "$myProcess"
}


#
# Set up our environment
#

njobs=$1
NPhotonsPerVoxel=$2
fclIn=$3
if [ "$4" ]; then
  if [ "$4" == "true" || "$4" == "True" || "$4" == "TRUE" ]; then
    makeupControl=true
  fi
else
  makeupControl==false
fi

label=${CLUSTER}_$(printf '%04d' $PROCESS)
## This sets all the needed FW and SRT and LD_LIBRARY_PATH envt variables. 
## Then we cd back to our TMP area. EC, 23-Nov-2010.

# In each voxel, run this many photons:


umask 0002

export GROUP=dune
export HOME=$CONDOR_DIR_ROOT
export CONDOR_SCRATCH_DIR=$TEMP

#we make a log here, just in case we need information from the process reassignment. This will make two log files for each makeup run. I can live with that.
#Don't forget to modify the makeup function to have the right array of missing job numbers.
if [ "$makeupControl" == "true" ]; then
  echo "FIRST"
  echo "$PROCESS"
  LOG=${CONDOR_DIR_LOG}/pd_library_gen_${label}.log
  #LOG=./pd_library_gen_${label}.log
  touch ${LOG}
  PROCESS="$(SetProcessNumber ${LOG} $PROCESS)"
  echo "SECOND"
  echo "$PROCESS"
fi

#need to update the label after we change the process
label=${CLUSTER}_$(printf '%04d' $PROCESS)
#LOG=./pd_library_gen_${label}.log
#FCL=./pd_library_gen_${label}.fcl
LOG=${CONDOR_DIR_LOG}/pd_library_gen_${label}.log
FCL=${CONDOR_DIR_FCL}/pd_library_gen_${label}.fcl
RTF=${CONDOR_DIR_ROOT}/pd_library_gen_${label}.root

echo "$LOG"
echo "$FCL"
echo "$RTF"


#
# Prepare to run
#
cd $CONDOR_DIR_ROOT
touch ${LOG}


echo "Writing log to ${LOG}"        1>>${LOG} 2>&1 
echo "Root file is ${RTF}"          1>>${LOG} 2>&1
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
