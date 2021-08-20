art_make(NO_PLUGINS
  EXCLUDE
    ScintTimeXeDoping_tool.cc
  LIB_LIBRARIES
    ROOT::Core
    ART_FRAMEWORK_SERVICES_REGISTRY
    CANVAS
    cetlib
    cetlib_except
    fhiclcpp
    CLHEP
    ART_UTILITIES
  )


cet_build_plugin(ScintTimeXeDoping art::tool
                larsim_PhotonPropagation_ScintTimeTools
                messagefacility::MF_MessageLogger
                FHICLCPP
                CLHEP
              )


install_headers()
install_fhicl()
install_source()
