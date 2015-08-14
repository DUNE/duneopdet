#!/bin/bash

script=OpticalLibraryBuild_Grid_dune.sh
outdir=/pnfs/dune/scratch/users/ahimmel/OpticalLibraries/OpticalLib_dune10kt_v1_workspace/
fcl=$PWD/dune10kt_workspace_buildopticallibrary_grid.fcl

clientargs="--resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --OS=SL6 --group=dune -f $fcl --role=Analysis"
toolsargs="-q -g --opportunistic --OS=SL6 "
fileargs="-dROOT $outdir/root -dFCL $outdir/fcl -dLOG $outdir/log "

larsoft="$MRB_TOP $MRB_PROJECT $MRB_PROJECT_VERSION $MRB_QUALS"

#Test job - jobsub_tools
#njobs=216000
#thisjob="-T $PWD/$script $njobs"
#echo "jobsub $toolsargs $fileargs $thisjob"
#jobsub $fileargs $thisjob

#real job - jobsub_tools
#njobs=500
#thisjob="-N $njobs $PWD/$script $njobs"
#echo "jobsub $toolsargs $fileargs $thisjob"
#jobsub $fileargs $thisjob

#Test job - jobsub_client
#njobs=216000
#thisjob="-T file://$PWD/$script $njobs"
#echo "jobsub_submit $clientargs $fileargs $thisjob"
#jobsub_submit $clientargs $fileargs $thisjob

#Test job 1 - jobsub_client
#njobs=7200
#nphotons=10
#thisjob="-M -N 1 file://$PWD/$script $njobs $nphotons"

#Test job 2 - jobsub_client
#njobs=500
#nphotons=1000
#thisjob="-M -N 1 file://$PWD/$script $njobs $nphotons"

#Real job - jobsub_client
njobs=500
nphotons=50000
thisjob="-N $njobs file://$PWD/$script $njobs $nphotons"

echo "jobsub_submit $clientargs $fileargs $thisjob $larsoft"
jobsub_submit $clientargs $fileargs $thisjob $larsoft
