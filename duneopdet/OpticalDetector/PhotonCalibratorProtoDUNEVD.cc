// dunetpc includes
#include "duneopdet/OpticalDetector/PhotonCalibratorProtoDUNEVD.h"

// LArSoft Includes
#include "larreco/Calibrator/IPhotonCalibrator.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "larcore/Geometry/WireReadout.h"

// ART includes
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <vector>
#include <algorithm>
#include <map>

namespace calib {

  PhotonCalibratorProtoDUNEVD::PhotonCalibratorProtoDUNEVD(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg)
  {

    fBadChannels = pset.get<std::vector<int> >("BadChannels");
    fHdwChannels = pset.get<std::vector<int> >("HdwChannels");
    fSPEAreas    = pset.get<std::vector<float> >("SPEAreas");
    fSPEShifts   = pset.get<std::vector<float> >("SPEShifts");
   
    //Map for associating each hardware channel with its correspondent SPE area 
    for(size_t i=0; i<fHdwChannels.size();i++){
      fAreaMap.emplace(fHdwChannels[i], fSPEAreas[i]);
      fShiftMap.emplace(fHdwChannels[i], fSPEShifts[i]);
    }
  
  }

  double PhotonCalibratorProtoDUNEVD::PE(double adcs, int opchannel) const
  {
    if (std::find(fBadChannels.begin(), fBadChannels.end(), opchannel) != fBadChannels.end()) {
      mf::LogDebug("PhotonCalibratorProtoDUNEVD") << "Skipping bad channel " << opchannel;
      return 0;
    }
    auto entry = fAreaMap.find(opchannel);
    float area = entry->second;
    auto entry2 = fShiftMap.find(opchannel);
    float shift = entry2->second;
    return adcs/area + shift;
  }

}

//DEFINE_ART_SERVICE_INTERFACE_IMPL(calib::PhotonCalibratorProtoDUNEVD, calib::IPhotonCalibrator)
