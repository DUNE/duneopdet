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




preareas  = [ 15, 30, 45, 60 ]
prenoises = [ 10, 100, 1000 ]
presnrs   = [ 4, 5, 7 ]
signal = 18.18

areas = {}
noises = {}
snrs = {}

for area in preareas:
    for noise in prenoises:
        for snr in presnrs:
            tag = "{0:02d}cm{1:04d}Hz{2:1d}snr".format(area, noise, snr)
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
    print "      ophit{0}:     @local::dunefd_ophit".format(tag)
    print "      opflash{0}:   @local::dunefd_opflash".format(tag)

print """
   }
   
   analyzers: {
"""
for tag in tags:
    print "      flashmatch{0}:  @local::marley_flashmatchana".format(tag)

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
print "   end_paths: [" + ", ".join(endpaths) + "]"
print "}"



for tag in tags:
    QE = areas[tag] * 0.00287 / 4.05 # Convert effective area to QE
    print "## Configs for {0}".format(tag)
    print "physics.producers.opdigi{0}.QEOverride:             {1:.6f}".format(tag, QE)
    print "physics.producers.opdigi{0}.DarkNoiseRate:          {1:d} #Hz".format(tag, noises[tag])
    print "physics.producers.opdigi{0}.LineNoiseRMS:           {1:.3f}".format(tag, signal/snrs[tag])
    print "physics.producers.ophit{0}.InputModule:             opdigi{0}".format(tag)
    print "physics.producers.opflash{0}.InputModule:           ophit{0}".format(tag)
    print "physics.analyzers.flashmatch{0}.OpHitModuleLabel:   ophit{0}".format(tag)
    print "physics.analyzers.flashmatch{0}.OpFlashModuleLabel: opflash{0}".format(tag)
    print
    
    
