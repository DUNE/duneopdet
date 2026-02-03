////////////////////////////////////////////////////////////////////////
// \file PhotonCalibratorServiceDUNE.h
//
// \brief Framework interface to PhotonCalibratorDUNE
//
// \author lpaulucc@fnal.gov based on PhotonCalibratorServiceProtuDUNESP by ahimmel@fnal.gov
//
////////////////////////////////////////////////////////////////////////


#ifndef PHOTONCALIBRATORSERVICEDUNE
#define PHOTONCALIBRATORSERVICEDUNE

// LArSoft Includes
#include "duneopdet/OpticalDetector/PhotonCalibratorDUNE.h"
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

  
  class PhotonCalibratorServiceDUNE : public IPhotonCalibratorService
  {
  public:
    using provider_type = PhotonCalibratorDUNE;
    
    struct ServiceConfiguration_t
    {
      fhicl::Sequence<int> HdwChannels  { fhicl::Name("HdwChannels")  };
      fhicl::Sequence<float> SPEAreas  { fhicl::Name("SPEAreas")  };
      fhicl::Sequence<float> SPEShifts { fhicl::Name("SPEShifts") };
      fhicl::Sequence<int> BadChannels { fhicl::Name("BadChannels"), fhicl::Comment("Channels to remove from reconstruction")  };
    };
    
    using Parameters = art::ServiceTable<ServiceConfiguration_t>;

  public:
    PhotonCalibratorServiceDUNE(Parameters const & config,
                                       art::ActivityRegistry& aReg)
      : fProvider( new PhotonCalibratorDUNE(config.get_PSet(), aReg) )
      { }
    
    provider_type const* provider() const override { return fProvider.get(); }

  private:
    std::unique_ptr<PhotonCalibratorDUNE> fProvider;
  };

}

DECLARE_ART_SERVICE_INTERFACE_IMPL(calib::PhotonCalibratorServiceDUNE, 
                                   calib::IPhotonCalibratorService, 
                                   LEGACY)

#endif // PHOTONCALIBRATORSERVICEDUNE
