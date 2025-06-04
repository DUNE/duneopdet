art_make(BASENAME_ONLY ALLOW_UNDERSCORES
  LIBRARY_NAME LowEPDSUtils
  LIB_LIBRARIES
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
    dunecore::ProducerUtils
  
  SERVICE_LIBRARIES
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
    dunecore::ProducerUtils
  
  MODULE_LIBRARIES LowEPDSUtils
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
    dunecore::ProducerUtils
  )

add_subdirectory(fcl)

install_headers()
install_fhicl()
install_scripts()
install_source()
  
