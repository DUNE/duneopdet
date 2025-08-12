include (phot::OpticalPath)

cet_build_plugin(DUNEHD10ktOpticalPath phot::OpticalPath
  LIBRARIES PRIVATE
  larcoreobj::geo_vectors
)

cet_build_plugin(ProtoDUNEHDOpticalPath phot::OpticalPath
  LIBRARIES PRIVATE
  larcoreobj::geo_vectors
)

install_headers()
install_fhicl()
install_source()