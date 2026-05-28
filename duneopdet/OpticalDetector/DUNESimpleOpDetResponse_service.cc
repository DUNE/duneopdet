////////////////////////////////////////////////////////////////////////
//
//  \file DUNESimpleOpDetResponse_service.cc
//
////////////////////////////////////////////////////////////////////////

#include "larana/OpticalDetector/OpDetResponseInterface.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"

#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// CLHEP includes
#include "CLHEP/Random/RandFlat.h"

namespace opdet {
  class DUNESimpleOpDetResponse : public OpDetResponseInterface {
  public:
    DUNESimpleOpDetResponse(fhicl::ParameterSet const& pset);

  private:
    void doReconfigure(fhicl::ParameterSet const& p) override;
    bool doDetected(int OpChannel, const sim::OnePhoton& Phot, int& newOpChannel) const override;
    bool doDetectedLite(int OpChannel, int& newOpChannel) const override;

    double fScintPreScale = 1.0;
    double fQE = 1.0;
    double fAcceptFraction = 1.0;
    
  }; // class DUNESimpleOpDetResponse

}

DECLARE_ART_SERVICE_INTERFACE_IMPL(opdet::DUNESimpleOpDetResponse,
                                   opdet::OpDetResponseInterface,
                                   LEGACY)

namespace opdet {
  //--------------------------------------------------------------------
  DUNESimpleOpDetResponse::DUNESimpleOpDetResponse(fhicl::ParameterSet const& pset)
  {
    this->doReconfigure(pset);
  }

  //--------------------------------------------------------------------
  void DUNESimpleOpDetResponse::doReconfigure(fhicl::ParameterSet const& pset)
  {
    auto const* LarProp = lar::providerFrom<detinfo::LArPropertiesService>();
    fScintPreScale = LarProp->ScintPreScale();
    fQE = pset.get<double>("QuantumEfficiency");
    if (fScintPreScale <= 0) {
      mf::LogError("DUNESimpleOpDetResponse_service")
        << "A prescale of " << LarProp->ScintPreScale()
        << " has been applied during optical MC production, but it is less than or equal to zero.";
      throw cet::exception("DUNESimpleOpDetResponse_service") << "Invalid scintillation prescale value";
    }
    fAcceptFraction = fQE / fScintPreScale;
    if (fAcceptFraction > 1.0001 || fAcceptFraction < 0)
      {
        mf::LogError("DUNESimpleOpDetResponse_service")
          << "A prescale of " << LarProp->ScintPreScale()
          << " has been applied and the Quantum efficiency is " << fQE << " QE/Prescale needs to be between 0 and 1";
        throw cet::exception("DUNESimpleOpDetResponse_service") << "Invalid QE and prescale combination";       
      }
  }

  //--------------------------------------------------------------------
  bool DUNESimpleOpDetResponse::doDetected(int OpChannel,
                                           const sim::OnePhoton& /*Phot*/,
                                           int& newOpChannel) const
  {
    return doDetectedLite(OpChannel, newOpChannel);
  }

  //--------------------------------------------------------------------
  bool DUNESimpleOpDetResponse::doDetectedLite(int OpChannel, int& newOpChannel) const
  {
    newOpChannel = OpChannel;
    if (CLHEP::RandFlat::shoot(1.0) < fAcceptFraction) return true;
    return false;
  }

} // namespace

DEFINE_ART_SERVICE_INTERFACE_IMPL(opdet::DUNESimpleOpDetResponse, opdet::OpDetResponseInterface)
