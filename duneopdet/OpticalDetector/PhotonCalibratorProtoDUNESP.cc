// dunetpc includes
#include "dune/OpticalDetector/PhotonCalibratorProtoDUNESP.h"

// LArSoft Includes
#include "larreco/Calibrator/IPhotonCalibrator.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()

// ART includes
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
   
#include <vector>


namespace calib {

  PhotonCalibratorProtoDUNESP::PhotonCalibratorProtoDUNESP(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg)
  {
    // Get the geometry service for getting max number of channels
    auto const& geometry(*lar::providerFrom< geo::Geometry >());

    // Initialize the SPE vectors with default values
    for (unsigned int channel = 0; channel < geometry.MaxOpChannel(); channel++) {
      fSPESizes.push_back(1.);
      fSPEShifts.push_back(0.);
    }

  }



  double PhotonCalibratorProtoDUNESP::PE(double adcs, int opchannel) const
  {
    return adcs/fSPESizes[opchannel] + fSPEShifts[opchannel];
  }

}

//DEFINE_ART_SERVICE_INTERFACE_IMPL(calib::PhotonCalibratorProtoDUNESP, calib::IPhotonCalibrator)
