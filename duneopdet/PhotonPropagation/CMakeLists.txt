
art_make(
          SERVICE_LIBRARIES larsim_PhotonPropagation
                        larsim_Simulation nug4::ParticleNavigation lardataobj_Simulation
                        larevt_Filters
	                lardataobj_RawData
	                larcorealg_Geometry
	                larcore_Geometry_Geometry_service
                        ART_FRAMEWORK_CORE
			ART_FRAMEWORK_PRINCIPAL
			ART_FRAMEWORK_SERVICES_REGISTRY
			ART_ROOT_IO_TFILE_SUPPORT
			ROOT::Core
			ART_ROOT_IO_TFILESERVICE_SERVICE
                        ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE
			art::Persistency_Common
			art::Persistency_Provenance
			art::Utilities
                        canvas::canvas
			messagefacility::MF_MessageLogger
			fhiclcpp::fhiclcpp
                        CLHEP
			cetlib::cetlib cetlib_except
			ROOT_BASIC_LIB_LIST
                        ROOT_EG
          MODULE_LIBRARIES
                        ROOT_BASIC_LIB_LIST
                        larsim_LegacyLArG4
                        lardataobj_Simulation
                        dune_PhotonPropagation_PhotonVisibilityServiceS2_service
                        larsim_Simulation nug4::ParticleNavigation lardataobj_Simulation
                        larcorealg_Geometry
      	                larcore_Geometry_Geometry_service
                        ART_FRAMEWORK_CORE
			ART_FRAMEWORK_PRINCIPAL
			ART_FRAMEWORK_SERVICES_REGISTRY
			ART_ROOT_IO_TFILE_SUPPORT
			ROOT::Core
			ART_ROOT_IO_TFILESERVICE_SERVICE
                        ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE
			nurandom::RandomUtils_NuRandomService_service
			art::Persistency_Common
			art::Persistency_Provenance
			art::Utilities
                        canvas::canvas
			messagefacility::MF_MessageLogger
			fhiclcpp::fhiclcpp
			cetlib::cetlib cetlib_except
                        CLHEP
	 )

install_headers()
install_fhicl()
install_source()

add_subdirectory(LibraryBuildTools)
add_subdirectory(ScintTimeTools)
