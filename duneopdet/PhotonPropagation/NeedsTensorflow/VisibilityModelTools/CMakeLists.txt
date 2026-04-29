cet_make_library(
  LIBRARY_NAME VisibilityModel 
  INTERFACE 
  SOURCE VisibilityModel.hh
  LIBRARIES
  larcorealg::CoreUtils
  larcorealg::Geometry
  larcoreobj::geo_vectors
  fhiclcpp::fhiclcpp
)

cet_write_plugin_builder(
  phot::vismodel::IVisibilityModel art::tool Modules
  INSTALL_BUILDER)

message(STATUS "CMAKE_MODULE_PATH is ${CMAKE_MODULE_PATH}")

include (phot::vismodel::IVisibilityModel)

cet_build_plugin( SemiAnalyticalVisModel phot::vismodel::IVisibilityModel
  LIBRARIES PRIVATE
  VisibilityModel
  larsim::OpticalPath
  larsim::PhotonPropagation_PhotonVisibilityService_service
  duneopdet::PhotonPropagation_PhotonVisibilityServiceS2_service
  larcore::Geometry_Geometry_service
)

cet_build_plugin( ComputableGraphVisModel phot::vismodel::IVisibilityModel
  LIBRARIES PRIVATE
  VisibilityModel
  larsimdnn::TFLoaderTool
  TensorFlow::cc
  TensorFlow::framework
  larcore::Geometry_Geometry_service
)

cet_build_plugin( PhotonLibraryVisModel phot::vismodel::IVisibilityModel
  LIBRARIES PRIVATE
  VisibilityModel
  art::Persistency_Common
  art::Framework_Core
  art::Framework_Principal
  art::Framework_Services_Registry
  art::Utilities
  larcore::Geometry_Geometry_service
  larsim::PhotonPropagation
  larsim::PhotonPropagation_PhotonVisibilityService_service
  duneopdet::PhotonPropagation_PhotonVisibilityServiceS2_service
)

install_headers()
install_fhicl()
install_source()
