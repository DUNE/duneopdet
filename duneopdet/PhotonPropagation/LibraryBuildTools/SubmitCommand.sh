#!/bin/bash

script=OpticalLibraryBuild_Grid_dune.sh
outdir=/pnfs/lbne/scratch/users/ahimmel/OpticalLibraries/OpticalLib_protoDune_v1_workspace/
fcl=$outdir/protoDune_buildopticallibrary_grid.fcl

clientargs="--resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --OS=SL6 --group=lbne -f $fcl --role=Analysis"
toolsargs="-q -g --opportunistic --OS=SL6 "
fileargs="-dROOT $outdir/root -dFCL $outdir/fcl -dLOG $outdir/log "

larsoft="$MRB_TOP dunetpc $MRB_PROJECT_VERSION $MRB_QUALS"

#Test job 1 - jobsub_client
#njobs=7200
#nphotons=10
#thisjob="-M -N 1 file://$PWD/$script $njobs $nphotons"

#Real job - jobsub_client
njobs=500
nphotons=50000
thisjob="-N $njobs file://$PWD/$script $njobs $nphotons"

echo "jobsub_submit $clientargs $fileargs $thisjob $larsoft"
jobsub_submit $clientargs $fileargs $thisjob $larsoft
