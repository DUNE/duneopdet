////////////////////////////////////////////////////////////////////////
// \file PhotonCalibratorProtoDUNESP.h
//
// \brief ProtoDUNESP service provider applying a flat scale factor to all optical hits.
//
// \author ahimmel@fnal.gov
//
////////////////////////////////////////////////////////////////////////


#ifndef PHOTONCALIBRATORPROTODUNESP_H
#define PHOTONCALIBRATORPROTODUNESP_H

#include "larreco/Calibrator/IPhotonCalibrator.h"

// ART includes
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"


namespace calib {

  class PhotonCalibratorProtoDUNESP : public IPhotonCalibrator
  {
  public:
    PhotonCalibratorProtoDUNESP(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);

    // Override base class functions
    double PE(double adcs, int opchannel) const override;
    bool   UseArea() const override { return true; } // ProtoDUNE always uses area

    /// Need a 3D position because result depends on position along length of
    /// bar. This is going to be pretty imprecise even so.
    // virtual double GeV(double PE, int opchannel, TVector3 pos) override;

  private:
    std::vector<float>  fSPESizes;
    std::vector<float>  fSPEShifts;

    
  }; // class PhotonCalibratorProtoDUNESP
}

//DECLARE_ART_SERVICE_INTERFACE_IMPL(calib::PhotonCalibratorProtoDUNESP, calib::IPhotonCalibrator, LEGACY)


#endif
