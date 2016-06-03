art_make( BASENAME_ONLY
          LIBRARY_NAME larana_OpticalDetector
          LIB_LIBRARIES 
                        lardata_RecoBase
                        larcore_Geometry
                        larcore_Geometry_Geometry_service
                        nusimdata_SimulationBase
                        ${ART_FRAMEWORK_CORE}
			                  ${ART_FRAMEWORK_PRINCIPAL}
			                  ${ART_FRAMEWORK_SERVICES_REGISTRY}
			                  ${ART_FRAMEWORK_SERVICES_BASIC}
			                  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			                  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			                  art_Persistency_Common canvas_Persistency_Common
			                  art_Utilities canvas_Utilities
			                  ${MF_MESSAGELOGGER}
			                  ${MF_UTILITIES}
	                      ${FHICLCPP}
			                  ${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
          MODULE_LIBRARIES 
	                larana_OpticalDetector
	  		larana_OpticalDetector_OpHitFinder
                        larcore_Geometry
                        larcore_Geometry_Geometry_service
                        nusimdata_SimulationBase
                        ${ART_FRAMEWORK_CORE}
			                  ${ART_FRAMEWORK_PRINCIPAL}
			                  ${ART_FRAMEWORK_SERVICES_REGISTRY}
			                  ${ART_FRAMEWORK_SERVICES_BASIC}
			                  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			                  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			                  art_Persistency_Common canvas_Persistency_Common
			                  art_Utilities canvas_Utilities
			                  ${MF_MESSAGELOGGER}
			                  ${MF_UTILITIES}
	                      ${FHICLCPP}
			                  ${CETLIB}
                        ${CLHEP}
                        ${ROOT_BASIC_LIB_LIST}
          SERVICE_LIBRARIES
                        larcore_Geometry
                        larcore_Geometry_Geometry_service
                        ${ART_FRAMEWORK_CORE}
			                  ${ART_FRAMEWORK_PRINCIPAL}
                  			${ART_FRAMEWORK_SERVICES_REGISTRY}
			                  ${ART_FRAMEWORK_SERVICES_BASIC}
			                  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			                  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			                  art_Persistency_Common canvas_Persistency_Common
			                  art_Utilities canvas_Utilities
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
            
