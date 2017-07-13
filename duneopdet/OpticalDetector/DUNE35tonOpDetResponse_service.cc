// -*- mode: c++; c-basic-offset: 4; -*-
////////////////////////////////////////////////////////////////////////
//
//  \file DUNE35tonOpDetResponse_service.cc
//
////////////////////////////////////////////////////////////////////////


#include "dune/OpticalDetector/DUNE35tonOpDetResponse.h"
#include "TGeoNode.h"
#include "TGeoBBox.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "CLHEP/Random/RandFlat.h"


namespace opdet{


    //--------------------------------------------------------------------
    DUNE35tonOpDetResponse::DUNE35tonOpDetResponse(fhicl::ParameterSet const& pset, 
                                         art::ActivityRegistry &/*reg*/)
    {
        this->doReconfigure(pset);
    }
    
    //--------------------------------------------------------------------
    DUNE35tonOpDetResponse::~DUNE35tonOpDetResponse() throw()
    { }


    //--------------------------------------------------------------------
    void DUNE35tonOpDetResponse::doReconfigure(fhicl::ParameterSet const& pset)
    {
        double tempfQE =         pset.get<double>("QuantumEfficiency");
        fWavelengthCutLow =      pset.get<double>("WavelengthCutLow");
        fWavelengthCutHigh =     pset.get<double>("WavelengthCutHigh");
        fLightGuideAttenuation = pset.get<bool>("LightGuideAttenuation");
        fChannelConversion =     pset.get<std::string>("ChannelConversion");
        std::string tmpAxis =    pset.get<std::string>("LongAxis"); 

        boost::algorithm::to_lower(tmpAxis);

        if (tmpAxis == "x") fLongAxis = 0;
        if (tmpAxis == "y") fLongAxis = 1;
        if (tmpAxis == "z") fLongAxis = 2;
        
        // Only allow channel conversion once - so it must be set to happen
        // either during full simulation (library generation) or during
        // fast simulation (library use).
        
        boost::algorithm::to_lower(fChannelConversion);
        
        fFullSimChannelConvert = false;
        fFastSimChannelConvert = false;
        
        if (fChannelConversion == "full") fFullSimChannelConvert = true;
        if (fChannelConversion == "fast") fFastSimChannelConvert = true;

        // Correct out the prescaling applied during simulation
        auto const *LarProp = lar::providerFrom<detinfo::LArPropertiesService>();
        fQE = tempfQE / LarProp->ScintPreScale();
        
        if (fQE > 1.0001 ) {
            mf::LogError("DUNE35tonOpDetResponse_service") << "Quantum efficiency set in OpDetResponse_service, " << tempfQE
                                                      << " is too large.  It is larger than the prescaling applied during simulation, "
                                                      << LarProp->ScintPreScale()
                                                      << ".  Final QE must be equalt to or smaller than the QE applied at simulation time.";
            assert(false);
        }

    }


    //--------------------------------------------------------------------
    int  DUNE35tonOpDetResponse::doNOpChannels() const
    {
        art::ServiceHandle<geo::Geometry> geom;
        if (fFastSimChannelConvert || fFullSimChannelConvert)
            return geom->NOpChannels();
        else
            return geom->NOpDets();

    }


    //--------------------------------------------------------------------
    bool DUNE35tonOpDetResponse::doDetected(int OpDet, const sim::OnePhoton& Phot, int &newOpChannel) const
    {
        
        // Find the Optical Detector using the geometry service
        art::ServiceHandle<geo::Geometry> geom;
        const TGeoNode* node = geom->OpDetGeoFromOpDet(OpDet).Node();

        // Identify the photon detector type
        int pdtype;
        std::string detname = node->GetName();
        boost::to_lower(detname);
        if      (detname.find("bar") != std::string::npos )     pdtype = 0;
        else if (detname.find("fiber") != std::string::npos )   pdtype = 1;
        else if (detname.find("plank") != std::string::npos )   pdtype = 2;
        else                                                    pdtype = -1;


        if (fFullSimChannelConvert){
            // Override default number of channels for Fiber and Plank
            float NOpHardwareChannels = geom->NOpHardwareChannels(OpDet);
            if (pdtype == 1) NOpHardwareChannels = 3;
            if (pdtype == 2) NOpHardwareChannels = 2;
            
            int hardwareChannel = (int) ( CLHEP::RandFlat::shoot(1.0) * NOpHardwareChannels );
            newOpChannel = geom->OpChannel(OpDet, hardwareChannel);
        }
        else{
            newOpChannel = OpDet;
        }
        
        // Check QE
        if ( CLHEP::RandFlat::shoot(1.0) > fQE ) return false;

        double wavel = wavelength(Phot.Energy);
        // Check wavelength acceptance
        if (wavel < fWavelengthCutLow) return false;
        if (wavel > fWavelengthCutHigh) return false;

        if (fLightGuideAttenuation) {
            // Get the length of the photon detector
            TGeoBBox *box = (TGeoBBox*)node->GetVolume()->GetShape();
            double opdetHalfLength = 0;
            double sipmDistance = 0;

            if (fLongAxis == 0) {
                opdetHalfLength = box->GetDX();
                sipmDistance = opdetHalfLength - Phot.FinalLocalPosition.x();
            }
            else if (fLongAxis == 1) {
                opdetHalfLength = box->GetDY();
                sipmDistance = opdetHalfLength - Phot.FinalLocalPosition.y();
            }
            else if (fLongAxis == 2) {
                opdetHalfLength = box->GetDZ();
                sipmDistance = opdetHalfLength - Phot.FinalLocalPosition.z();
            }
            else {
                mf::LogError("DUNE35tonOpDetResponse") << "Unknown axis, fLongAxis = " << fLongAxis;
                assert(false);
            }



            if (pdtype == 0) {
                double normalize   =  0.6961; // Normalize mean performance to be the same in all PDs
                double lambdaShort =  5.56; // cm
                double lambdaLong  =  44.13; // cm
                double fracShort   =  0.40 * normalize;
                double fracLong    =  0.60 * normalize;

                double AttenuationProb = fracShort*exp(-sipmDistance/lambdaShort) + fracLong*exp(-sipmDistance/lambdaLong);

                //std::cout <<"GRAPH: " << OpDet << "," << Phot.FinalLocalPosition.y() << "," << sipmDistance << "," << AttenuationProb << std::endl;
                
                //mf::LogVerbatim("DUNE35tonOpDetResponse") << "OpDet: " << OpDet << " is a " << pdtype 
                //                                     << " with length " << opdetHalfLength << " in detector "
                //                                     << box->GetDX() << " x " << box->GetDY()  << " x " << box->GetDZ()
                //                                     << " named " << detname;
                //mf::LogVerbatim("DUNE35tonOpDetResponse") << "   Local Position = (" << Phot.FinalLocalPosition.x() 
                //                                     << ", " << Phot.FinalLocalPosition.y() << ", " << Phot.FinalLocalPosition.z() << ")";
                //mf::LogVerbatim("DUNE35tonOpDetResponse") << "   Distance to SiPM = " << sipmDistance << " along axis " << fLongAxis;
                //mf::LogVerbatim("DUNE35tonOpDetResponse") << "   Attenuation Probability = " << AttenuationProb;
                
                // Throw away some photons based on attenuation
                if ( CLHEP::RandFlat::shoot(1.0) > AttenuationProb ) return false;
            }
            else if (pdtype == 1) {
                double normalize   = 1.0; // Normalize mean performance to be the same in all PDs
                double lambda      = 14.6; // cm

                // Throw away some photons based on attenuation
                double AttenuationProb = normalize*exp(-sipmDistance/lambda);

                //std::cout <<"GRAPH: " << OpDet << "," << Phot.FinalLocalPosition.y() << "," << sipmDistance << "," << AttenuationProb << std::endl;

                //mf::LogVerbatim("DUNE35tonOpDetResponse") << "OpDet: " << OpDet << " is a " << pdtype 
                //                                     << " with length " << opdetHalfLength << " in detector "
                //                                     << box->GetDX() << " x " << box->GetDY()  << " x " << box->GetDZ()
                //                                     << " named " << detname;
                //mf::LogVerbatim("DUNE35tonOpDetResponse") << "   Local Position = (" << Phot.FinalLocalPosition.x() 
                //                                     << ", " << Phot.FinalLocalPosition.y() << ", " << Phot.FinalLocalPosition.z() << ")";
                //mf::LogVerbatim("DUNE35tonOpDetResponse") << "   Distance to SiPM = " << sipmDistance << " along axis " << fLongAxis;
                //mf::LogVerbatim("DUNE35tonOpDetResponse") << "   Attenuation Probability = " << AttenuationProb;
                
                // Throw away some photons based on attenuation
                if ( CLHEP::RandFlat::shoot(1.0) > AttenuationProb ) return false;
            }
            else if (pdtype == 2) {
                double normalize   = 0.4498; // Normalize mean performance to be the same in all PDs
                double lambda      = 48.4; // cm
                double altDistance = 2*opdetHalfLength - sipmDistance;
                double frac        = 0.5 * normalize;

                double AttenuationProb = frac*exp(-sipmDistance/lambda) + frac*exp(-altDistance/lambda);

                //std::cout <<"GRAPH: " << OpDet << "," << Phot.FinalLocalPosition.y() << "," << sipmDistance << "," << AttenuationProb << std::endl;
                
                //mf::LogVerbatim("DUNE35tonOpDetResponse") << "OpDet: " << OpDet << " is a " << pdtype 
                //                                     << " with length " << opdetHalfLength << " in detector "
                //                                     << box->GetDX() << " x " << box->GetDY()  << " x " << box->GetDZ()
                //                                     << " named " << detname;
                //mf::LogVerbatim("DUNE35tonOpDetResponse") << "   Local Position = (" << Phot.FinalLocalPosition.x() 
                //                                     << ", " << Phot.FinalLocalPosition.y() << ", " << Phot.FinalLocalPosition.z() << ")";
                //mf::LogVerbatim("DUNE35tonOpDetResponse") << "   Distance to SiPM = " << sipmDistance << " along axis " << fLongAxis;
                //mf::LogVerbatim("DUNE35tonOpDetResponse") << "   Attenuation Probability = " << AttenuationProb;
                                
                // Throw away some photons based on attenuation
                if ( CLHEP::RandFlat::shoot(1.0) > AttenuationProb ) return false;
            }
            else {
                mf::LogWarning("DUNE35tonOpDetResponse") << "OpDet: " << OpDet << " is an unknown PD type named: " << detname 
                                                    << ". Assuming no attenuation.";
            }

        }

        return true;
    }

    //--------------------------------------------------------------------
    bool DUNE35tonOpDetResponse::doDetectedLite(int OpDet, int &newOpChannel) const
    {
        if (fFastSimChannelConvert){

            // Find the Optical Detector using the geometry service
            art::ServiceHandle<geo::Geometry> geom;
            // Here OpDet must be opdet since we are introducing
            // channel mapping here.
            const TGeoNode* node = geom->OpDetGeoFromOpDet(OpDet).Node();

            
            // Identify the photon detector type
            int pdtype;
            std::string detname = node->GetName();
            boost::to_lower(detname);
            if      (detname.find("bar") != std::string::npos )   pdtype = 0;
            else if (detname.find("fiber") != std::string::npos ) pdtype = 1;
            else if (detname.find("plank") != std::string::npos ) pdtype = 2;
            else                                                  pdtype = -1;

            // Override default number of channels for Fiber and Plank
            float NOpHardwareChannels = geom->NOpHardwareChannels(OpDet);
            if (pdtype == 1) NOpHardwareChannels = 8;
            if (pdtype == 2) NOpHardwareChannels = 2;

            int hardwareChannel = (int) ( CLHEP::RandFlat::shoot(1.0) * NOpHardwareChannels );
            newOpChannel = geom->OpChannel(OpDet, hardwareChannel);
        }
        else{
            newOpChannel = OpDet;
        }
        
        // Check QE
        if ( CLHEP::RandFlat::shoot(1.0) > fQE ) return false;

        
        return true;
    }

} // namespace

DEFINE_ART_SERVICE_INTERFACE_IMPL(opdet::DUNE35tonOpDetResponse, opdet::OpDetResponseInterface)

