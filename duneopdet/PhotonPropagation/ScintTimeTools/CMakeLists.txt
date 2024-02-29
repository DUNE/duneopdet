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
    CLHEP::CLHEP
    art::Utilitiesa
  ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types


cet_build_plugin(ScintTimeXeDoping art::tool
                messagefacility::MF_MessageLogger
                fhiclcpp::fhiclcpp
                CLHEP::CLHEP
                larsim::PhotonPropagation
              )


install_headers()
install_fhicl()
install_source()
