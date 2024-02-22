
art_make(
          SERVICE_LIBRARIES larsim_PhotonPropagation
                        larsim::Simulation nug4::ParticleNavigation 
			lardataobj::Simulation
                        larevt::Filters
	                lardataobj::RawData
	                larcorealg::Geometry
	                larcore::Geometry_Geometry_service
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support
                        ROOT::Core
                        art_root_io::TFileService_service
                        art::Framework_Services_Optional_RandomNumberGenerator_service
			art::Persistency_Common
			art::Persistency_Provenance
			art::Utilities
                        canvas::canvas
			messagefacility::MF_MessageLogger
			fhiclcpp::fhiclcpp
                        CLHEP
			cetlib::cetlib cetlib_except
                        ROOT::EG
          MODULE_LIBRARIES
                        larsim::LegacyLArG4
                        lardataobj::Simulation
                        duneopdet::PhotonPropagation_PhotonVisibilityServiceS2_service
                        larsim::Simulation 
                        nug4::ParticleNavigation lardataobj_Simulation
                        larcorealg::Geometry
      	                larcore::Geometry_Geometry_service
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support
                        ROOT::Core
                        art_root_io::TFileService_service
                        art::Framework_Services_Optional_RandomNumberGenerator_service
			nurandom::RandomUtils_NuRandomService_service
			art::Persistency_Common
			art::Persistency_Provenance
			art::Utilities
                        canvas::canvas
			messagefacility::MF_MessageLogger
			fhiclcpp::fhiclcpp
			cetlib::cetlib cetlib_except
                        CLHEP
	 ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_fhicl()
install_source()

add_subdirectory(LibraryBuildTools)
add_subdirectory(ScintTimeTools)
