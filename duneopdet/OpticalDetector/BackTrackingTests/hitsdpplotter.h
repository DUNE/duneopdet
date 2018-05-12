BEGIN_PROLOG

standard_hitsdpplotter:
{
  module_type:	HitSdpPlotter

}

END_PROLOG

#include "services_dune.fcl"

process_name: hitsdpplotter

services:
{
  @table::dunefd_services
  TFileService:          { fileName: "hitsdpplotter.root" }
  TimeTracker:           {}
  MemoryTracker:         {}
  RandomNumberGenerator: {}
}

physics:
{
  analyzers:
  {
    hitsdpplotter:     @local::standard_hitsdpplotter
  }

  ana: [ hitsdpplotter ]
  end_paths: [ ana ]

}

source:
{
  module_type: RootInput
  maxEvents:  -1      # Number of events to create
}

