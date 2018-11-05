#!/usr/bin/env python

print """
#include "services_dune.fcl"
#include "singles_dune.fcl"
#include "largeantmodules_dune.fcl"
#include "detsimmodules_dune.fcl"
#include "mccheatermodules.fcl"
#include "photpropservices_dune.fcl"
#include "opticaldetectormodules_dune.fcl"
#include "opticaldetectorservices_dune.fcl"
#include "FlashMatchAna.fcl"


process_name: OpticalResim

services:
{
  # Load the service that manages root files for histograms.
  TFileService: { fileName: "dune1x2x6_optical_tutorial_resimulate_hist.root" }
  TimeTracker:       {}
  RandomNumberGenerator: {} #ART native random number generator
  message:      @local::standard_info
  @table::dunefd_simulation_services
}

# DUNE FD 1x2x6 workspace geometry
services.Geometry:                @local::dune10kt_1x2x6_geo
services.PhotonVisibilityService: @local::dune10kt_1x2x6_photonvisibilityservice
services.OpDigiProperties:        @local::dunefd_opdigiproperties


source:
{
  module_type: RootInput
  maxEvents:  -1        # Run over all events
  #specify from command line with -s or --source
   
}


outputs:
{
   out1:
   {
      module_type: RootOutput
      fileName:    "dune1x2x6_optical_tutorial_resimulate_gen.root"
      #default file name, can override from command line with -o or --output
   }
}

"""




preareas  = [ 15, 30, 60 ]
prenoises = [ 10, 300, 1000 ]
#prenoises = [ 100 ]
presnrs   = [ 3, 4, 7 ]
prerefs  = [ "Opt", "Pes" ]
signal = 18.18

defareas  = [ 45 ]
defnoises = [ 100 ]
defsnrs   = [ 5 ]
defrefs   = [ "Non" ]

areas = {}
noises = {}
snrs = {}
reflected = {}


sets = [ ("DEF", defareas, defnoises, defsnrs, defrefs),
         ("EFF", preareas, defnoises, defsnrs, defrefs),
         ("NSE", defareas, prenoises, defsnrs, defrefs),
         ("SNR", defareas, defnoises, presnrs, defrefs),
         ("REF", defareas, defnoises, defsnrs, prerefs) ]
    


for ( name, iareas, inoises, isnrs, irefs ) in sets:
    for area in iareas:
        for noise in inoises:
            for snr in isnrs:
                for ref in irefs:
                    tag = "{0}{1:02d}cm{2:04d}Hz{3:1d}snr{4}Refl".format(name, area, noise, snr, ref)
                    if ref == "Opt":
                        reflected[tag] = 1.44
                    elif ref == "Pes":
                        reflected[tag] = 0.84
                        area /= 2
                    areas[tag] = area
                    noises[tag] = noise
                    snrs[tag] = snr

tags = sorted(areas.keys())



print """
physics:
{

   # Run both detector simulation and reconstruction
   producers:
   {
"""
for tag in tags:
    print "      opdigi{0}:    @local::dunefd_opdigi_threegang".format(tag)
for tag in tags:
    print "      ophit{0}:     @local::dunefd_ophit".format(tag)
for tag in tags:
    print "      opflash{0}:   @local::dunefd_opflash".format(tag)

print """
   }
   
   analyzers: {
"""
for tag in tags:
    print "      flashmatch{0}:  @local::marley_flashmatchana".format(tag)
    #print "      flashmatch{0}:  @local::standard_flashmatchana".format(tag)

print """
   }

"""
simpaths = []
for tag in tags:
    print "   simPath{0}: [ opdigi{0}, ophit{0}, opflash{0}]".format(tag)
    simpaths.append("simPath{0}".format(tag))
endpaths = []
for tag in tags:
    print "   anaPath{0}: [ flashmatch{0} ]".format(tag)
    endpaths.append("anaPath{0}".format(tag))

print "   stream1:  [ out1 ]"
print ""
print "   trigger_paths: [" + ", ".join(simpaths) + "]"
print "   end_paths: [" + ", ".join(endpaths) + ", stream1 ]"
print "}"



for tag in tags:
    QE = areas[tag] * 0.00287 / 4.05 # Convert effective area to QE
    print "## Configs for {0}".format(tag)

    
    #print "physics.producers.opdigi{0}.SSP_LED_DigiTree:       true".format(tag)
    print "physics.producers.opdigi{0}.QEOverride:             {1:.6f}".format(tag, QE)

    if tag in reflected:
        refQE = QE*reflected[tag]
        print "physics.producers.opdigi{0}.QERefOverride:          {1:.6f}".format(tag, refQE)

    print "physics.producers.opdigi{0}.DarkNoiseRate:          {1:d} #Hz".format(tag, noises[tag])
    print "physics.producers.opdigi{0}.LineNoiseRMS:           {1:.3f}".format(tag, signal/snrs[tag])
    print "physics.producers.ophit{0}.InputModule:             opdigi{0}".format(tag)
    print "physics.producers.opflash{0}.InputModule:           ophit{0}".format(tag)
    print "physics.analyzers.flashmatch{0}.OpDetWaveformLabel: opdigi{0}".format(tag)
    print "physics.analyzers.flashmatch{0}.OpHitModuleLabel:   ophit{0}".format(tag)
    print "physics.analyzers.flashmatch{0}.OpFlashModuleLabel: opflash{0}".format(tag)
    print
    
    
