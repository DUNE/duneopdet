////////////////////////////////////////////////////////////////////////
// \file ProtoDUNEOpDetResponse.h
//
// \brief service containing information about the response of optical detectors in ProtoDUNE
//
// \author L. Paulucci (adapted from DUNEOpDetResponse module)
//
////////////////////////////////////////////////////////////////////////

#ifndef PROTODUNE_OPDET_RESPONSE_H
#define PROTODUNE_OPDET_RESPONSE_H

// LArSoft includes
#include "lardataobj/Simulation/SimPhotons.h"
#include "larana/OpticalDetector/OpDetResponseInterface.h"
#include "dune/OpticalDetector/DUNEOpDetResponseInterface.h"

#include "TMath.h"
#include "TFile.h"
#include "TH1D.h"
#include "TVector3.h"
#include <vector>
#include <fstream>

namespace opdet
{
  class ProtoDUNEOpDetResponse : public opdet::OpDetResponseInterface {
    public:

      ProtoDUNEOpDetResponse(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
      ~ProtoDUNEOpDetResponse() throw();

      // virtual bool detectedLite(int OpDet, int &newOpChannel, int& hardwareChannel) const;
      //bool ProtoDUNEOpDetResponse::detectedLite(int OpDet, int &newOpChannel, int& hardwareChannel) const

      bool detectedLiteWithChannel(int OpDet, int &newOpChannel, int& hardwareChannel) const
      {
        return doDetectedLiteWithChannel( OpDet, newOpChannel, hardwareChannel);
      }



    private:

      virtual void doReconfigure(fhicl::ParameterSet const& p);

      virtual int  doNOpChannels() const;
      virtual bool doDetected(int OpDet, const sim::OnePhoton& Phot, int &newOpChannel) const;
      virtual bool doDetectedLite(int OpDet, int &newOpChannel) const;
    //bool ProtoDUNEOpDetResponse::doDetectedLite(int OpDet, int &newOpChannel, int &hardwareChannel) const
      bool doDetectedLiteWithChannel(int OpDet, int &newOpChannel, int& hardwareChannel) const;

      float fQE;                     // Quantum efficiency of paddle
      float fQEArapucaBeam;          // Quantum efficiency of Arapuca bar on beam-side (mean eff for all windows)
      float fQEArapucaNonBeam;       // Quantum efficiency of Arapuca bar on non-beam side (mean eff for all windows)

      float fWavelengthCutLow;       // Sensitive wavelength range
      float fWavelengthCutHigh;      //

      bool fLightGuideAttenuation;   // Flag to turn on position-dependent sensitivity
      double lambdaShort;
      double lambdaLong;
      double fracShort;
      double fracLong;

      //TY: Comment out the 5 unused variables
      //bool fWireAttenuation;             // Flag to turn on wire/mesh transmission probability
      //TH1D* hm;                          // histogram to store transmission probability
      //TVector3* Dist;
      //TVector3* VersorX;
      //TVector3* VersorY;
      std::string fWireTransmissionFile;

      std::string fChannelConversion;
      bool fFullSimChannelConvert;   // Flag to conver detector->electronics channels in full optical sim
      bool fFastSimChannelConvert;   // Flag to conver detector->electronics channels in fast optical sim

      int fLongAxis;                 // 0 = x, 1 = y, 2 = z

  }; // class ProtoDUNEOpDetResponse


} //namespace opdet


DECLARE_ART_SERVICE_INTERFACE_IMPL(opdet::ProtoDUNEOpDetResponse, opdet::OpDetResponseInterface, LEGACY)

#endif //OPDET_RESPONSE_PROTODUNE_H
