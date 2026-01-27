////////////////////////////////////////////////////////////////////////
// \file PhotonCalibratorProtoDUNEVD.h
//
// \brief ProtoDUNEVD service provider applying a channel based scale factor to optical hits.
//
// \author lpaulucc@fnal.gov based on PhotonCalibratorProtoduneSP by ahimmel
//
////////////////////////////////////////////////////////////////////////


#ifndef PHOTONCALIBRATORPROTODUNEVD_H
#define PHOTONCALIBRATORPROTODUNEVD_H

#include "larreco/Calibrator/IPhotonCalibrator.h"

// ART includes
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"


namespace calib {

  class PhotonCalibratorProtoDUNEVD : public IPhotonCalibrator
  {
  public:
    PhotonCalibratorProtoDUNEVD(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);

    // Override base class functions
    double PE(double adcs, int opchannel) const override;
    bool   UseArea() const override { return true; } // ProtoDUNE always uses area

    /// Need a 3D position because result depends on position along length of
    /// bar. This is going to be pretty imprecise even so.
    // virtual double GeV(double PE, int opchannel, TVector3 pos) override;

  private:
    std::vector<int>        fHdwChannels;
    std::vector<float>      fSPEAreas;
    std::vector<float>      fSPEShifts;
    std::vector<int>        fBadChannels;
    std::map<int, float>    fAreaMap;
    std::map<int, float>    fShiftMap;
  }; // class PhotonCalibratorProtoDUNEVD
}

//DECLARE_ART_SERVICE_INTERFACE_IMPL(calib::PhotonCalibratorProtoDUNEVD, calib::IPhotonCalibrator, LEGACY)


#endif
