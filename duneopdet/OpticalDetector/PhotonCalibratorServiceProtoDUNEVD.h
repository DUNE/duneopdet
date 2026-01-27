////////////////////////////////////////////////////////////////////////
// \file PhotonCalibratorServiceProtoDUNEVD.h
//
// \brief Framework interface to PhotonCalibratorProtoDUNEVD
//
// \author lpaulucc@fnal.gov based on PhotonCalibratorServiceProtuDUNESP by ahimmel@fnal.gov
//
////////////////////////////////////////////////////////////////////////


#ifndef PHOTONCALIBRATORSERVICEPROTODUNEVD
#define PHOTONCALIBRATORSERVICEPROTODUNEVD

// LArSoft Includes
#include "duneopdet/OpticalDetector/PhotonCalibratorProtoDUNEVD.h"
#include "larreco/Calibrator/IPhotonCalibratorService.h"


// ART Includes
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Services/Registry/ServiceTable.h" 
#include "art/Framework/Principal/Run.h"


#include <vector>

namespace calib {

  
  class PhotonCalibratorServiceProtoDUNEVD : public IPhotonCalibratorService
  {
  public:
    using provider_type = PhotonCalibratorProtoDUNEVD;
    
    struct ServiceConfiguration_t
    {
      fhicl::Sequence<int> HdwChannels  { fhicl::Name("HdwChannels")  };
      fhicl::Sequence<float> SPEAreas  { fhicl::Name("SPEAreas")  };
      fhicl::Sequence<float> SPEShifts { fhicl::Name("SPEShifts") };
      fhicl::Sequence<int> BadChannels { fhicl::Name("BadChannels"), fhicl::Comment("Channels to remove from reconstruction")  };
    };
    
    using Parameters = art::ServiceTable<ServiceConfiguration_t>;

  public:
    PhotonCalibratorServiceProtoDUNEVD(Parameters const & config,
                                       art::ActivityRegistry& aReg)
      : fProvider( new PhotonCalibratorProtoDUNEVD(config.get_PSet(), aReg) )
      { }
    
    provider_type const* provider() const override { return fProvider.get(); }

  private:
    std::unique_ptr<PhotonCalibratorProtoDUNEVD> fProvider;
  };

}

DECLARE_ART_SERVICE_INTERFACE_IMPL(calib::PhotonCalibratorServiceProtoDUNEVD, 
                                   calib::IPhotonCalibratorService, 
                                   LEGACY)

#endif // PHOTONCALIBRATORSERVICEPROTODUNEVD
