// -*- mode: c++; c-basic-offset: 4; -*-
////////////////////////////////////////////////////////////////////////
//
//  \file DUNEOpDetResponse_service.cc
//
////////////////////////////////////////////////////////////////////////


#include "dune/OpticalDetector/DUNEOpDetResponse.h"
#include "TGeoNode.h"
#include "TGeoBBox.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "CLHEP/Random/RandFlat.h"


namespace opdet{


    //--------------------------------------------------------------------
    DUNEOpDetResponse::DUNEOpDetResponse(fhicl::ParameterSet const& pset,
                                         art::ActivityRegistry &/*reg*/)
    {
        this->doReconfigure(pset);
    }

    //--------------------------------------------------------------------
    DUNEOpDetResponse::~DUNEOpDetResponse() throw()
    { }


    //--------------------------------------------------------------------
    void DUNEOpDetResponse::doReconfigure(fhicl::ParameterSet const& pset)
    {
        double tempfQE =         pset.get<double>("QuantumEfficiency");
        fWavelengthCutLow =      pset.get<double>("WavelengthCutLow");
        fWavelengthCutHigh =     pset.get<double>("WavelengthCutHigh");
        fLightGuideAttenuation = pset.get<bool>("LightGuideAttenuation");
        lambdaShort =            pset.get<double>("LambdaShort");
        lambdaLong =             pset.get<double>("LambdaLong");
        fracShort =              pset.get<double>("FracShort");
        fracLong =               pset.get<double>("FracLong");
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
            mf::LogError("DUNEOpDetResponse_service") << "Quantum efficiency set in OpDetResponse_service, " << tempfQE
                                                      << " is too large.  It is larger than the prescaling applied during simulation, "
                                                      << LarProp->ScintPreScale()
                                                      << ".  Final QE must be equalt to or smaller than the QE applied at simulation time.";
            assert(false);
        }

    }


    //--------------------------------------------------------------------
    int  DUNEOpDetResponse::doNOpChannels() const
    {
        art::ServiceHandle<geo::Geometry> geom;
        if (fFastSimChannelConvert || fFullSimChannelConvert)
            return geom->NOpChannels();
        else
            return geom->NOpDets();

    }


    //--------------------------------------------------------------------
    bool DUNEOpDetResponse::doDetected(int OpDet, const sim::OnePhoton& Phot, int &newOpChannel) const
    {

        // Find the Optical Detector using the geometry service
        art::ServiceHandle<geo::Geometry> geom;


        if (fFullSimChannelConvert){
            // Override default number of channels for Fiber and Plank
            float NOpHardwareChannels = geom->NOpHardwareChannels(OpDet);
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
            const TGeoNode* node = geom->OpDetGeoFromOpDet(OpDet).Node();
            TGeoBBox *box = (TGeoBBox*)node->GetVolume()->GetShape();
            double opdetLength = 0;
            double sipmDistance = 0;

            if (fLongAxis == 0) {
                opdetLength = box->GetDX();
                sipmDistance = opdetLength - Phot.FinalLocalPosition.x();
            }
            else if (fLongAxis == 1) {
                opdetLength = box->GetDY();
                sipmDistance = opdetLength - Phot.FinalLocalPosition.y();
            }
            else if (fLongAxis == 2) {
                opdetLength = box->GetDZ();
                sipmDistance = opdetLength - Phot.FinalLocalPosition.z();
            }
            else {
                mf::LogError("DUNEOpDetResponse") << "Unknown axis, fLongAxis = " << fLongAxis;
                assert(false);
            }



            // Throw away some photons based on attenuation
            double AttenuationProb = fracShort*exp(-sipmDistance/lambdaShort) + fracLong*exp(-sipmDistance/lambdaLong);

            //mf::LogVerbatim("DUNEOpDetResponse") << "OpDet: " << OpDet
            //                                     << " has length " << opdetLength << " in detector "
            //                                     << box->GetDX() << " x " << box->GetDY()  << " x " << box->GetDZ();
            //mf::LogVerbatim("DUNEOpDetResponse") << "   Local Position = (" << Phot.FinalLocalPosition.x()
            //                                     << ", " << Phot.FinalLocalPosition.y() << ", " << Phot.FinalLocalPosition.z() << ")";
            //mf::LogVerbatim("DUNEOpDetResponse") << "   Distance to SiPM = " << sipmDistance << " along axis " << fLongAxis;
            //mf::LogVerbatim("DUNEOpDetResponse") << "   Attenuation Probability = " << AttenuationProb;

            if ( CLHEP::RandFlat::shoot(1.0) > AttenuationProb ) return false;


        }

        return true;
    }

    //--------------------------------------------------------------------
    bool DUNEOpDetResponse::doDetectedLite(int OpDet, int &newOpChannel) const
    {
        if (fFastSimChannelConvert){

            // Find the Optical Detector using the geometry service
            art::ServiceHandle<geo::Geometry> geom;
            // Here OpDet must be opdet since we are introducing
            // channel mapping here.
            float NOpHardwareChannels = geom->NOpHardwareChannels(OpDet);
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

    //--------------------------------------------------------------------
    bool DUNEOpDetResponse::doDetectedLiteWithChannel(int OpDet, int &newOpChannel, int &hardwareChannel) const
    {
        if (fFastSimChannelConvert){

            // Find the Optical Detector using the geometry service
            art::ServiceHandle<geo::Geometry> geom;
            // Here OpDet must be opdet since we are introducing
            // channel mapping here.
            float NOpHardwareChannels = geom->NOpHardwareChannels(OpDet);
            hardwareChannel = (int) ( CLHEP::RandFlat::shoot(1.0) * NOpHardwareChannels );
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

DEFINE_ART_SERVICE_INTERFACE_IMPL(opdet::DUNEOpDetResponse, opdet::OpDetResponseInterface)

