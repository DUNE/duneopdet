# use modern cetmodules build commands

# PhotonLibraryPropagationS2_module.cc depends on LegacyLArG4 and cannot be built with geant 4.11

cet_build_plugin(PhotonVisibilityServiceS2 art::service
  LIBRARIES PUBLIC
                        larsim::PhotonPropagation
                        larsim::Simulation nug4::ParticleNavigation 
			lardataobj::Simulation
                        lardata::headers
                        larevt::Filters
	                lardataobj::RawData
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
			art::Persistency_Common
			art::Persistency_Provenance
			art::Utilities
                        canvas::canvas
			messagefacility::MF_MessageLogger
			fhiclcpp::fhiclcpp
                        CLHEP::CLHEP
			cetlib::cetlib cetlib_except::cetlib_except
                        ROOT::EG
)


install_headers()
install_fhicl()
install_source()

# NeedsTensorflow depends on LegacyLArG4 and cannot be built with geant 4.11
##if(DEFINED ENV{TENSORFLOW_DIR})
##  add_subdirectory(NeedsTensorflow)
##endif()
add_subdirectory(LibraryBuildTools)
add_subdirectory(ScintTimeTools)
add_subdirectory(OpticalPathTools)
