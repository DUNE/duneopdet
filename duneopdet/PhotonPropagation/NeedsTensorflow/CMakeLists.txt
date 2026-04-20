add_subdirectory( VisibilityModelTools )

art_make(
  MODULE_LIBRARIES
    larcorealg::Geometry
    larcore::Geometry_Geometry_service
    art::Framework_Core
    art::Framework_Principal
    art::Framework_Services_Registry
    art_root_io::tfile_support
    ROOT::Core
    ROOT::Tree
    art_root_io::TFileService_service
    art::Framework_Services_Optional_RandomNumberGenerator_service
    nurandom::RandomUtils_NuRandomService_service
    art::Persistency_Common
    art::Persistency_Provenance
    art::Utilities
    canvas::canvas
    messagefacility::MF_MessageLogger
    fhiclcpp::fhiclcpp
    cetlib::cetlib cetlib_except::cetlib_except
    CLHEP::CLHEP

    VisibilityModel
    larsim::OpticalPath
    larsim::PhotonPropagation
    larsimdnn::TFLoaderTool
    TensorFlow::cc
    TensorFlow::framework
    larsim::PhotonPropagation
    larsim::PhotonPropagation_PhotonVisibilityService_service
    duneopdet::PhotonPropagation_PhotonVisibilityServiceS2_service
)

install_headers()
install_fhicl()
install_source()
