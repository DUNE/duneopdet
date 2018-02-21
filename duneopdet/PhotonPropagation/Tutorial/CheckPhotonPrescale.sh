#!/bin/bash

fhicl=$1
artfile=$2

artscale=`config_dumper -S $artfile | grep "ScintPreScale" | awk '{ print $2 }'`
fhiclscale=`ART_DEBUG_CONFIG=1 lar -c $fhicl 2>&1 | grep "ScintPreScale" | awk '{ print $2 }'`
fhiclqe=`ART_DEBUG_CONFIG=1 lar -c $fhicl  2>&1 | grep -m 1 "QuantumEfficiency" | awk '{ print $2 }'`

if [[ $artscale != $fhiclscale ]]; then
    echo "You need to set the prescale in your fhicl to match this file:"
    echo "services.LArProperties.ScintPreScale: $artscale"
    exit 1
fi

scalem=${fhiclscale/[eE]/*10^}
qem=${fhiclqe/[eE]/*10^}

if (( $(echo "$scalem < $qem" | bc -l) )); then
    echo "Your quantum efficiency, $fhiclqe, is too high. It needs to"
    echo "be below the Scintillation pre-scale $fhiclscale"
    exit 1
fi

echo "Good to go!"

