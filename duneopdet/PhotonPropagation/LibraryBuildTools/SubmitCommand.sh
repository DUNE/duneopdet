#!/bin/bash

script=OpticalLibraryBuild_Grid_dune.sh
outdir=/pnfs/lbne/scratch/users/ahimmel/OpticalLibraries/OpticalLib_dune10kt_v1_1x2x6
fcl=$outdir/dune10kt_v1_1x2x6_buildopticallibrary_grid.fcl

# 57600 seconds = 16 hours, but it will not be a sharp cut-off
clientargs="--resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --OS=SL6 --group=dune -f $fcl --role=Analysis --memory=2500MB --expected-lifetime=0s --append_condor_requirements='((TARGET.GLIDEIN_ToDie-CurrentTime)>57600)'"
toolsargs="-q -g --opportunistic --OS=SL6 "
fileargs="-dROOT $outdir/root -dFCL $outdir/fcl -dLOG $outdir/log "

larsoft="$MRB_TOP dunetpc $MRB_PROJECT_VERSION $MRB_QUALS"

#Test job 1 - jobsub_client
#njobs=7200
#nphotons=10
#clientargs="$clientargs --expected-lifetime=600"
#thisjob="-M -N 1 file://$PWD/$script $njobs $nphotons"

#Real job - jobsub_client
njobs=1500
nphotons=50000
thisjob="-N $njobs file://$PWD/$script $njobs $nphotons"

echo "jobsub_submit $clientargs $fileargs $thisjob $larsoft"
jobsub_submit $clientargs $fileargs $thisjob $larsoft
