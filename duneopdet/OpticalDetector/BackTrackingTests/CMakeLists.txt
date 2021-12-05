
art_make( BASENAME_ONLY MODULE_LIBRARIES
  larcore_Geometry_Geometry_service
  larcorealg_Geometry
  lardataobj_RawData
  lardataobj_AnalysisBase
  dune_DuneObj
  larsim_MCCheater_BackTrackerService_service
  larsim_MCCheater_ParticleInventoryService_service
  larsim_MCCheater_PhotonBackTrackerService_service
  larsim_Simulation nug4::ParticleNavigation lardataobj_Simulation
  nusimdata::SimulationBase	
  ART_FRAMEWORK_CORE
 ART_FRAMEWORK_PRINCIPAL
  ART_FRAMEWORK_SERVICES_REGISTRY
  ART_ROOT_IO_TFILE_SUPPORT
  ROOT::Core
  ART_ROOT_IO_TFILESERVICE_SERVICE
  art::Persistency_Common canvas
  art::Persistency_Provenance canvas
  art::Utilities canvas
  messagefacility::MF_MessageLogger
  cetlib::cetlib cetlib_except
  ROOT_BASIC_LIB_LIST
  )

install_fhicl()
install_headers()
install_source()
