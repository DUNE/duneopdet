
art_make(
          SERVICE_LIBRARIES larsim_PhotonPropagation
                        larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
                        larevt_Filters
	                lardataobj_RawData
	                larcorealg_Geometry
	                larcore_Geometry_Geometry_service
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			art_Persistency_Common
			art_Persistency_Provenance
			art_Utilities
                        canvas
			${MF_MESSAGELOGGER}
			${FHICLCPP}
                        ${CLHEP}
			cetlib cetlib_except
			${ROOT_BASIC_LIB_LIST}
                        ${ROOT_EG}
          MODULE_LIBRARIES
                        ${ROOT_BASIC_LIB_LIST}
                        larsim_LArG4
                        lardataobj_Simulation
                        dune_PhotonPropagation_PhotonVisibilityServiceS2_service
                        larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
                        larcorealg_Geometry
      	                larcore_Geometry_Geometry_service
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			nutools_RandomUtils_NuRandomService_service
			art_Persistency_Common
			art_Persistency_Provenance
			art_Utilities
                        canvas
			${MF_MESSAGELOGGER}
			${FHICLCPP}
			cetlib cetlib_except
                        ${CLHEP}
	 )

install_headers()
install_fhicl()
install_source()

add_subdirectory(LibraryBuildTools)

