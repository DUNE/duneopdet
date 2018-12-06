////////////////////////////////////////////////////////////////////////
// \file PhotonCalibratorServiceProtoDUNESP.h
//
// \brief Framework interface to PhotonCalibratorProtoDUNESP
//
// \author ahimmel@fnal.gov
//
////////////////////////////////////////////////////////////////////////


#ifndef PHOTONCALIBRATORSERVICEPROTODUNESP
#define PHOTONCALIBRATORSERVICEPROTODUNESP

// LArSoft Includes
#include "dune/OpticalDetector/PhotonCalibratorProtoDUNESP.h"
#include "larreco/Calibrator/IPhotonCalibratorService.h"


// ART Includes
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Services/Registry/ServiceTable.h" 
#include "art/Framework/Principal/Run.h"


#include <vector>

namespace calib {

  
  class PhotonCalibratorServiceProtoDUNESP : public IPhotonCalibratorService
  {
  public:
    using provider_type = PhotonCalibratorProtoDUNESP;
    
    struct ServiceConfiguration_t
    {
      //fhicl::Atom<float> SPESize  { fhicl::Name("SPESize")  };
      //fhicl::Atom<float> SPEShift { fhicl::Name("SPEShift") };
      //fhicl::Atom<float> UseArea  { fhicl::Name("UseArea")  };
    };
    
    using Parameters = art::ServiceTable<ServiceConfiguration_t>;

  public:
    PhotonCalibratorServiceProtoDUNESP(Parameters const & config,
                                       art::ActivityRegistry& aReg)
      : fProvider( new PhotonCalibratorProtoDUNESP(config.get_PSet(), aReg) )
      { }
    
    provider_type const* provider() const override { return fProvider.get(); }

  private:
    std::unique_ptr<PhotonCalibratorProtoDUNESP> fProvider;
  };

}

DECLARE_ART_SERVICE_INTERFACE_IMPL(calib::PhotonCalibratorServiceProtoDUNESP, 
                                   calib::IPhotonCalibratorService, 
                                   LEGACY)

#endif // PHOTONCALIBRATORSERVICEPROTODUNESP
