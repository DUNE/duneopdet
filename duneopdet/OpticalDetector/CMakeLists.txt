art_make( BASENAME_ONLY
          LIBRARY_NAME OpticalDetector
          LIB_LIBRARIES 
                        RecoBase
                        Geometry
                        Geometry_service
                        ${SIMULATIONBASE}
                        ${ART_FRAMEWORK_CORE}
			                  ${ART_FRAMEWORK_PRINCIPAL}
			                  ${ART_FRAMEWORK_SERVICES_REGISTRY}
			                  ${ART_FRAMEWORK_SERVICES_BASIC}
			                  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			                  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			                  ${ART_PERSISTENCY_COMMON}
			                  ${ART_UTILITIES}
			                  ${MF_MESSAGELOGGER}
			                  ${MF_UTILITIES}
	                      ${FHICLCPP}
			                  ${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
          MODULE_LIBRARIES 
	                OpticalDetector
	  		OpHitFinder
                        Geometry
                        Geometry_service
                        LArProperties_service
                        ${SIMULATIONBASE}
                        ${ART_FRAMEWORK_CORE}
			                  ${ART_FRAMEWORK_PRINCIPAL}
			                  ${ART_FRAMEWORK_SERVICES_REGISTRY}
			                  ${ART_FRAMEWORK_SERVICES_BASIC}
			                  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			                  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			                  ${ART_PERSISTENCY_COMMON}
			                  ${ART_UTILITIES}
			                  ${MF_MESSAGELOGGER}
			                  ${MF_UTILITIES}
	                      ${FHICLCPP}
			                  ${CETLIB}
                        ${CLHEP}
                        ${ROOT_BASIC_LIB_LIST}
          SERVICE_LIBRARIES
                        Geometry
                        Geometry_service
                        ${ART_FRAMEWORK_CORE}
			                  ${ART_FRAMEWORK_PRINCIPAL}
                  			${ART_FRAMEWORK_SERVICES_REGISTRY}
			                  ${ART_FRAMEWORK_SERVICES_BASIC}
			                  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			                  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			                  ${ART_PERSISTENCY_COMMON}
			                  ${ART_UTILITIES}
			                  ${MF_MESSAGELOGGER}
			                  ${MF_UTILITIES}
	                      ${FHICLCPP}
			                  ${CETLIB}
                        ${CLHEP}
                        ${ROOT_BASIC_LIB_LIST}
         )

install_headers()
install_fhicl()
install_source()
install_fw( LIST sensl_1.dat
            SUBDIRNAME dune/OpticalDetector )
            
