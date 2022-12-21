// -*- mode: c++; c-basic-offset: 4; -*-
////////////////////////////////////////////////////////////////////////
//
//  \file ProtoDUNEOpDetResponse_service.cc
//
////////////////////////////////////////////////////////////////////////


#include "duneopdet/OpticalDetector/ProtoDUNEOpDetResponse.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "TGeoNode.h"
#include "TGeoBBox.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "CLHEP/Random/RandFlat.h"
#include <boost/algorithm/string.hpp>

namespace opdet{


    //--------------------------------------------------------------------
    ProtoDUNEOpDetResponse::ProtoDUNEOpDetResponse(fhicl::ParameterSet const& pset,
                                         art::ActivityRegistry &/*reg*/)
    {
        this->doReconfigure(pset);
    }

    //--------------------------------------------------------------------
    ProtoDUNEOpDetResponse::~ProtoDUNEOpDetResponse() throw()
    { }


    //--------------------------------------------------------------------
    void ProtoDUNEOpDetResponse::doReconfigure(fhicl::ParameterSet const& pset)
    {
        double tempfQE =         pset.get<double>("QuantumEfficiency");
        double tempfQEAraBeam =    pset.get<double>("QuantumEfficiencyArapucaBeam");
        double tempfQEAraNonBeam = pset.get<double>("QuantumEfficiencyArapucaNonBeam");
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
        fQEArapucaBeam = tempfQEAraBeam / LarProp->ScintPreScale();
        fQEArapucaNonBeam = tempfQEAraNonBeam / LarProp->ScintPreScale();

        if (fQE > 1.0001 || fQEArapucaBeam > 1.0001 || fQEArapucaNonBeam > 1.0001) {
            mf::LogError("ProtoDUNEOpDetResponse_service") << "Quantum efficiency set in OpDetResponse_service, " << tempfQE << ", " 
                                                        << tempfQEAraBeam << ", or " << tempfQEAraNonBeam
                                                      << " is too large.  It is larger than the prescaling applied during simulation, "
                                                      << LarProp->ScintPreScale()
                                                      << ".  Final QE must be equalt to or smaller than the QE applied at simulation time.";
            std::abort();
        }

    }


    //--------------------------------------------------------------------
    int  ProtoDUNEOpDetResponse::doNOpChannels() const
    {
        art::ServiceHandle<geo::Geometry> geom;
        if (fFastSimChannelConvert || fFullSimChannelConvert)
            return geom->NOpChannels();
        else
            return geom->NOpDets();

    }


    //--------------------------------------------------------------------
    bool ProtoDUNEOpDetResponse::doDetected(int OpDet, const sim::OnePhoton& Phot, int &newOpChannel) const
    {

        // Find the Optical Detector using the geometry service
        art::ServiceHandle<geo::Geometry> geom;

        if (fFullSimChannelConvert){
            // Override default number of channels for Fiber and Plank
            float NOpHardwareChannels = geom->NOpHardwareChannels(OpDet);
            int hardwareChannel = (int) ( CLHEP::RandFlat::shoot(1.0) * NOpHardwareChannels );
            newOpChannel = geom->OpChannel(OpDet, hardwareChannel);
        }else{
            newOpChannel = OpDet;
        }

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
                mf::LogError("ProtoDUNEOpDetResponse") << "Unknown axis, fLongAxis = " << fLongAxis;
                std::abort();
        }

        if(opdetLength > 50){ //i.e. if this is a paddle (this is half length in units of cm)
          // Check QE
          if ( CLHEP::RandFlat::shoot(1.0) > fQE ) return false;

          double wavel = wavelength(Phot.Energy);
          // Check wavelength acceptance
          if (wavel < fWavelengthCutLow) return false;
          if (wavel > fWavelengthCutHigh) return false;

          if (fLightGuideAttenuation) {
            // Throw away some photons based on attenuation
            double AttenuationProb = fracShort*exp(-sipmDistance/lambdaShort) + fracLong*exp(-sipmDistance/lambdaLong);

            if ( CLHEP::RandFlat::shoot(1.0) > AttenuationProb ) return false;
          }
        }else{ //if it is Arapuca
          if(Phot.FinalLocalPosition.x()<0){ //beam side
            // Check QE
            if ( CLHEP::RandFlat::shoot(1.0) > fQEArapucaBeam ) return false;
          }else{ //non-beam side
            // Check QE
            if ( CLHEP::RandFlat::shoot(1.0) > fQEArapucaNonBeam ) return false; 
          }
        }
    return true;
  }

    //--------------------------------------------------------------------
    bool ProtoDUNEOpDetResponse::doDetectedLite(int OpDet, int &newOpChannel) const
    {
        // Find the Optical Detector using the geometry service
        art::ServiceHandle<geo::Geometry> geom;

        if (fFastSimChannelConvert){
            // Here OpDet must be opdet since we are introducing
            // channel mapping here.
            float NOpHardwareChannels = geom->NOpHardwareChannels(OpDet);
            int hardwareChannel = (int) ( CLHEP::RandFlat::shoot(1.0) * NOpHardwareChannels );
            newOpChannel = geom->OpChannel(OpDet, hardwareChannel);
        }else{
            newOpChannel = OpDet;
        }

        double opdetLength = 0;
        const TGeoNode* node = geom->OpDetGeoFromOpDet(OpDet).Node();
        TGeoBBox *box = (TGeoBBox*)node->GetVolume()->GetShape();
        opdetLength = box->GetDX();
        if(opdetLength > 50){ //i.e. if this is a paddle (this is half length in units of cm)
        // Check QE
      	  if ( CLHEP::RandFlat::shoot(1.0) > fQE ) return false;
        }else{ //if it is Arapuca
          auto const xyz = geom->OpDetGeoFromOpChannel(OpDet).GetCenter();
          if(xyz.X()<0){ //beam side
            // Check QE
            if ( CLHEP::RandFlat::shoot(1.0) > fQEArapucaBeam ) return false;
          }else{ //non-beam side
            // Check QE
            if ( CLHEP::RandFlat::shoot(1.0) > fQEArapucaNonBeam ) return false;
          }
        }

        return true;
    }

    //--------------------------------------------------------------------
    bool ProtoDUNEOpDetResponse::doDetectedLiteWithChannel(int OpDet, int &newOpChannel, int &hardwareChannel) const
    {
        // Find the Optical Detector using the geometry service
        art::ServiceHandle<geo::Geometry> geom;

        if (fFastSimChannelConvert){
            // Here OpDet must be opdet since we are introducing
            // channel mapping here.
            float NOpHardwareChannels = geom->NOpHardwareChannels(OpDet);
            hardwareChannel = (int) ( CLHEP::RandFlat::shoot(1.0) * NOpHardwareChannels );
            newOpChannel = geom->OpChannel(OpDet, hardwareChannel);
        }else{
            newOpChannel = OpDet;
        }

        double opdetLength = 0;
        const TGeoNode* node = geom->OpDetGeoFromOpDet(OpDet).Node();
        TGeoBBox *box = (TGeoBBox*)node->GetVolume()->GetShape();
        opdetLength = box->GetDX();
        if(opdetLength > 50){ //i.e. if this is a paddle 
          // Check QE
          if ( CLHEP::RandFlat::shoot(1.0) > fQE ) return false;
        }else{ //if it is Arapuca
          auto const xyz = geom->OpDetGeoFromOpChannel(OpDet).GetCenter();
          if(xyz.X()<0){ //beam side
            // Check QE
            if ( CLHEP::RandFlat::shoot(1.0) > fQEArapucaBeam ) return false;
          }else{ //non-beam side
            // Check QE
            if ( CLHEP::RandFlat::shoot(1.0) > fQEArapucaNonBeam ) return false;
          }
        }

        return true;
    }


} // namespace

DEFINE_ART_SERVICE_INTERFACE_IMPL(opdet::ProtoDUNEOpDetResponse, opdet::OpDetResponseInterface)
