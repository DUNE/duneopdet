////////////////////////////////////////////////////////////////////////
// \file DUNE35tonOpDetResponse.h
//
// \brief service containing information about the response of optical detectors in DUNE35ton
//
// \author ahimmel@phy.duke.edu
//
////////////////////////////////////////////////////////////////////////

#ifndef DUNE35ton_OPDET_RESPONSE_H
#define DUNE35ton_OPDET_RESPONSE_H

// LArSoft includes
#include "lardataobj/Simulation/SimPhotons.h"
#include "larana/OpticalDetector/OpDetResponseInterface.h"



namespace opdet
{
    class DUNE35tonOpDetResponse : public opdet::OpDetResponseInterface {
    public:

        DUNE35tonOpDetResponse(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
        ~DUNE35tonOpDetResponse() throw();



    private:

        virtual void doReconfigure(fhicl::ParameterSet const& p);

        virtual int  doNOpChannels() const;
        virtual bool doDetected(int OpDet, const sim::OnePhoton& Phot, int &newOpChannel) const;
        virtual bool doDetectedLite(int OpDet, int &newOpChannel) const;

        float fQE;                     // Quantum efficiency of tube
        
        float fWavelengthCutLow;       // Sensitive wavelength range 
        float fWavelengthCutHigh;      // 
        
        bool fLightGuideAttenuation;   // Flag to turn on position-dependent sensitivity

        std::string fChannelConversion;
        bool fFullSimChannelConvert;   // Flag to conver detector->electronics channels in full optical sim
        bool fFastSimChannelConvert;   // Flag to conver detector->electronics channels in fast optical sim

        int fLongAxis;                 // 0 = x, 1 = y, 2 = z

    }; // class DUNE35tonOpDetResponse

    
} //namespace opdet


DECLARE_ART_SERVICE_INTERFACE_IMPL(opdet::DUNE35tonOpDetResponse, opdet::OpDetResponseInterface, LEGACY)

#endif //OPDET_RESPONSE_DUNE35ton_H
