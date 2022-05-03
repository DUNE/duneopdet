art_make(NO_PLUGINS
  EXCLUDE
    ScintTimeXeDoping_tool.cc
  LIB_LIBRARIES
    ROOT::Core
    art::Framework_Services_Registry
    canvas::canvas
    cetlib::cetlib
    cetlib_except::cetlib_except
    fhiclcpp::fhiclcpp
    CLHEP
    art::Utilities
  )


cet_build_plugin(ScintTimeXeDoping art::tool
                messagefacility::MF_MessageLogger
                fhiclcpp::fhiclcpp
                CLHEP
              )


install_headers()
install_fhicl()
install_source()
