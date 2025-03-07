add_subdirectory(fcl)

cet_make_library(
  LIBRARY_NAME
    SolarNuUtils
  SOURCE
    SolarAuxUtils.cc 
    AdjOpHitsUtils.cc
  LIBRARIES
    ${ROOT_CORE}
    ${ROOT_EVE}
    larcore::Geometry_Geometry_service
    larcorealg::Geometry
    lardata::Utilities
    lardataobj::RawData
    lardataobj::RecoBase
    larsim::MCCheater_BackTrackerService_service
    larsim::MCCheater_PhotonBackTrackerService_service
    larsim::MCCheater_ParticleInventoryService_service
    art::Persistency_Provenance
    art::Persistency_Common
    art::Utilities 
    art_root_io::tfile_support
    messagefacility::MF_MessageLogger
    cetlib::cetlib
  )

install_headers()
  
