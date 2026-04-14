// dunetpc includes
#include "duneopdet/OpticalDetector/PhotonCalibratorProtoDUNESP.h"

// LArSoft Includes
#include "larreco/Calibrator/IPhotonCalibrator.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()

// ART includes
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <vector>
#include <algorithm>

namespace calib {

  PhotonCalibratorProtoDUNESP::PhotonCalibratorProtoDUNESP(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg)
  {
    // Get the geometry service for getting max number of channels
    //auto const& geometry(*lar::providerFrom< geo::Geometry >());

    // Initialize the SPE vectors with default values
    //for (unsigned int channel = 0; channel < geometry.MaxOpChannel(); channel++) {
    // fSPESizes.push_back(1.);
    // fSPEShifts.push_back(0.);
    //}


    fBadChannels = pset.get<std::vector<int> >("BadChannels");
    auto chanNums = pset.get<std::vector<int> >("ChannelNumbers"); // to create a map
    unsigned int idx = 0;
    for (auto &ch: chanNums) {
        fChannelMap[ch] = idx;
        ++idx;
    }

    auto speSizes = pset.get<std::vector<float> >("SPESizes");
    auto speShifts = pset.get<std::vector<float> >("SPEShifts");
    if ( !speSizes.size() ) {
        // FIXME assert error. Don't use the calibrator if no calibrations provided.
    }
    fSPESizes = speSizes;
    if (speShifts.size() )
        fSPEShifts = speShifts;
    else
        // If no shifts provided, default to 0.
        fSPEShifts = std::vector<float>(speSizes.size(),0.);
  }
  double PhotonCalibratorProtoDUNESP::PE(double adcs, int opchannel) const
  {
    if (std::find(fBadChannels.begin(), fBadChannels.end(), opchannel) != fBadChannels.end()) {
      mf::LogDebug("PhotonCalibratorProtoDUNESP") << "Skipping bad channel " << opchannel;
      return 0;
    }

    unsigned int idx = opchannel;
    if (fChannelMap.size()) {
        idx = fChannelMap.at(opchannel);
    }

    return adcs/fSPESizes[idx] + fSPEShifts[idx];
  }

}

//DEFINE_ART_SERVICE_INTERFACE_IMPL(calib::PhotonCalibratorProtoDUNESP, calib::IPhotonCalibrator)
