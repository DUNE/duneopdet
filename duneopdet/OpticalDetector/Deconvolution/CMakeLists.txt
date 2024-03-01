
art_make( BASENAME_ONLY ALLOW_UNDERSCORES
                LIBRARY_NAME duneopdet_OpticalDetector
                LIB_LIBRARIES 
                                lardataalg::DetectorInfo
                                lardataobj::RecoBase
                                larcorealg::Geometry
                                larcore::Geometry_Geometry_service
                                larana::OpticalDetector_OpHitFinder
                                nusimdata::SimulationBase
                                fhiclcpp::fhiclcpp
                                messagefacility::MF_MessageLogger
                                art::Framework_Core
                                art::Framework_Principal
                                art::Framework_Services_Registry
                                art_root_io::tfile_support
                                ROOT::Core
                                art_root_io::TFileService_service
                                art::Framework_Services_Optional_RandomNumberGenerator_service
                                art::Persistency_Common
                                art::Utilities 
				canvas::canvas
                                cetlib::cetlib 
				cetlib_except::cetlib_except 
                                Boost::filesystem
  
                MODULE_LIBRARIES
                                dunecore::DuneObj
                                lardataalg::DetectorInfo
                                lardataobj::RecoBase
                                lardata::headers
                                larana::OpticalDetector
                                larana::OpticalDetector_OpHitFinder
                                larcorealg::Geometry
                                larcore::Geometry_Geometry_service
                                duneopdet::OpticalDetector
                                larsim::MCCheater_ParticleInventoryService_service
                                larsim::MCCheater_PhotonBackTrackerService_service 
                                nusimdata::SimulationBase
                                nurandom::RandomUtils_NuRandomService_service
                                fhiclcpp::fhiclcpp
                                messagefacility::MF_MessageLogger
                                art::Framework_Core
                                art::Framework_Principal
                                art::Framework_Services_Registry
                                art_root_io::tfile_support
                                ROOT::Core
                                art_root_io::TFileService_service
                                art::Framework_Services_Optional_RandomNumberGenerator_service
                                art::Persistency_Provenance
                                art::Persistency_Common
                                art::Utilities 
				canvas::canvas
                                cetlib::cetlib 
				cetlib_except::cetlib_except
				CLHEP::CLHEP
                                Boost::filesystem

                SERVICE_LIBRARIES
                                larcorealg::Geometry
                                larcore::Geometry_Geometry_service
                                duneopdet::OpticalDetector
                                fhiclcpp::fhiclcpp
                                messagefacility::MF_MessageLogger
                                art::Framework_Core
                                art::Framework_Principal
                                art::Framework_Services_Registry
                                art_root_io::tfile_support
                                ROOT::Core
                                art_root_io::TFileService_service
                                art::Persistency_Common
                                art::Utilities 
				canvas::canvas
                                cetlib::cetlib cetlib_except::cetlib_except
				CLHEP::CLHEP
                                Boost::filesystem

)

install_headers()
install_fhicl()
install_source()
