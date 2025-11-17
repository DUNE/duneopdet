art_make(LIB_LIBRARIES 
        lardataobj::RecoBase
        larcorealg::Geometry
        larcore::Geometry_Geometry_service
        larana::OpticalDetector_OpHitFinder
        larreco::PhotonCalibrator
        nusimdata::SimulationBase
        fhiclcpp::fhiclcpp
        messagefacility::MF_MessageLogger
        art::Framework_Core
        art::Framework_Principal
        art::Framework_Services_Registry
        art_root_io::tfile_support
        ROOT::Core
        art_root_io::TFileService_service
        art::Framework_Services_Optional_RandomNumberGenerator_service
        art::Persistency_Common
        art::Utilities 
	canvas::canvas
        cetlib::cetlib 
	cetlib_except::cetlib_except 
        Boost::filesystem
         MODULE_LIBRARIES
 
        duneopdet::OpticalDetector
	nusimdata::SimulationBase
#	dunereco::FDSensOpt_FDSensOptData
	art::Persistency_Provenance
	ROOT::Geom
	ROOT::XMLIO
	ROOT::Gdml
	ROOT::Core ROOT::Hist ROOT::Tree

  # lardataobj_Simulation
        SERVICE_LIBRARIES
        art_root_io::tfile_support

)

install_headers()
install_fhicl()
install_source()
