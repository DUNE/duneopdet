//=========================================================
// OpDetDigitizerDUNE_module.cc
// This module produces digitized output 
// (creating OpDetWaveform)
// from photon detectors taking SimPhotonsLite as input.
//
// Gleb Sinev, Duke, 2015
// Based on OpMCDigi_module.cc
//=========================================================

#ifndef OpDetDigitizerDUNE_h
#define OpDetDigitizerDUNE_h 1

// Framework includes

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Utilities/Exception.h"
#include "fhiclcpp/ParameterSet.h"

// ART extensions
#include "artextensions/SeedService/SeedService.hh"

// LArSoft includes

#include "Simulation/sim.h"
#include "Simulation/SimPhotons.h"
#include "Simulation/LArG4Parameters.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/TimeService.h"
#include "OpticalDetector/OpDetResponseInterface.h"
#include "RawData/OpDetWaveform.h"

// CLHEP includes

#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandExponential.h"
#include "CLHEP/Random/RandFlat.h"

// C++ includes

#include <vector>
#include <map>
#include <cmath>

namespace opdet {

  class OpDetDigitizerDUNE : public art::EDProducer{

    public:
      
      OpDetDigitizerDUNE(fhicl::ParameterSet const&);
      // Should the destructor be empty?
//      virtual ~OpDetDigitizerDUNE();
      
      void produce(art::Event&);

    private:

      // The parameters read from the FHiCL file
      std::string fInputModule; // Input tag for OpDet collection
      float fSampleFreq;        // Sampling frequency in MHz
      float fTimeBegin;         // Beginning of sample in us
      float fTimeEnd;           // End of sample in us
      float fVoltageToADC;      // Conversion factor mV to ADC counts
      float fLineNoise;         // Pedestal RMS in ADC counts
      float fDarkNoiseRate;     // In Hz
      float fCrossTalk;         // Probability of SiPM producing 2 PE signal
                                // in response to 1 photon
      short fPedestal;          // In ADC counts

      // Random number engines
      CLHEP::RandGauss       *fRandGauss;
      CLHEP::RandExponential *fRandExponential;
      CLHEP::RandFlat        *fRandFlat;

      // Function that adds n pulses to a waveform
      void AddPulse(int timeBin, int scale, std::vector< float >& waveform);

      // Functional response to one photoelectron (time in ns)
      float Pulse1PE(float time) const;

      // Single photoelectron pulse parameters
      float fPulseLength;   // 1PE pulse length in us
      float fPeakTime;      // Time when the pulse reaches its maximum in us
      float fMaxAmplitude;  // Maximum amplitude of the pulse in mV
      float fFrontTime;     // Constant in the exponential function 
                            // of the leading edge in us
      float fBackTime;      // Constant in the exponential function 
                            // of the tail in us

      std::vector< float > fSinglePEWaveform;
      void CreateSinglePEWaveform();

      // Values in this map determine whether we are going
      // to record the waveform for this optical channel or not
      std::map< int, bool > fRecordMap;

      // Produce waveform on one of the optical detectors
      void CreateOpDetWaveform(sim::SimPhotonsLite const&, 
                               opdet::OpDetResponseInterface const&,
                               std::vector< std::vector< float > >&);

      // Vary the pedestal
      void AddLineNoise(std::vector< std::vector< float > >&);

      void AddDarkNoise(std::vector< std::vector< float > >&, int);

      int CrossTalk() const;
  };

}

#endif

namespace opdet {

  DEFINE_ART_MODULE(OpDetDigitizerDUNE)

}

namespace opdet {
  
  //---------------------------------------------------------------------------
  // Constructor
  OpDetDigitizerDUNE::OpDetDigitizerDUNE(fhicl::ParameterSet const& pset)
  {

    // This module produces (infrastructure piece)
    produces< std::vector< raw::OpDetWaveform > >();

    // Read the fcl-file
    fInputModule   = pset.get< std::string >("InputModule"  );
    fVoltageToADC  = pset.get< float       >("VoltageToADC" );
    fLineNoise     = pset.get< float       >("LineNoise"    );
    fDarkNoiseRate = pset.get< float       >("DarkNoiseRate");
    fCrossTalk     = pset.get< float       >("CrossTalk"    );
    fPedestal      = pset.get< short       >("Pedestal"     );

    // Obtaining parameters from the TimeService
    art::ServiceHandle< util::TimeService > timeService;
    fSampleFreq = timeService->OpticalClock().Frequency();

    // Assume the readout starts at 0
    fTimeBegin  = 0.0;

    // Take the TPC readout window size and covert 
    // to us with the electronics clock frequency
    fTimeEnd    = art::ServiceHandle< util::DetectorProperties >()->
                  ReadOutWindowSize()/timeService->TPCClock().Frequency();

    
    // Initializing random number engines
    unsigned int seed = 
             pset.get< unsigned int >("Seed", sim::GetRandomNumberSeed());
    createEngine(seed);

    art::ServiceHandle< art::RandomNumberGenerator > rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine();
    fRandGauss        = new CLHEP::RandGauss(engine);
    fRandExponential  = new CLHEP::RandExponential(engine);
    fRandFlat         = new CLHEP::RandFlat(engine);

    // Creating a single photoelectron waveform
    // Hardcoded, probably need to read them from the FHiCL file
    fPulseLength  = 4.0;
    fPeakTime     = 0.260;
    fMaxAmplitude = 0.12;
    fFrontTime    = 0.009;
    fBackTime     = 0.476;
    CreateSinglePEWaveform();

  }
/*
  //---------------------------------------------------------------------------
  // Destructor
  OpDetDigitizerDUNE::~OpDetDigitizerDUNE()
  {
  }
 */
  //---------------------------------------------------------------------------
  void OpDetDigitizerDUNE::produce(art::Event& evt)
  {
    
    // A pointer that will store produced OpDetWaveforms
    std::unique_ptr< std::vector< raw::OpDetWaveform > > 
                      pulseVecPtr(new std::vector< raw::OpDetWaveform >);
    
    art::ServiceHandle< sim::LArG4Parameters > lgp;
    bool fUseLitePhotons = lgp->UseLitePhotons();

    // Total number of ticks in our readout
    int nSamples = (fTimeEnd - fTimeBegin)*fSampleFreq;

    // Service for determining opdet responses
    art::ServiceHandle< opdet::OpDetResponseInterface > odResponse;
    // Total number of optical channels
    int nOpChannels = odResponse->NOpChannels();

    // This vector stores waveforms we create for each optical detector
    std::vector< std::vector< float > > 
               opDetWaveforms(nOpChannels, 
                              std::vector< float >(nSamples,float(fPedestal)));

    if (fUseLitePhotons)
    {
      // Get SimPhotonsLite from the event
      art::Handle< std::vector< sim::SimPhotonsLite > > litePhotonHandle;
      evt.getByLabel(fInputModule, litePhotonHandle);

      // For every optical detector:
      for (auto const& opDet : (*litePhotonHandle)) 
        CreateOpDetWaveform(opDet, *odResponse, opDetWaveforms);
    }
    else throw art::Exception(art::errors::UnimplementedFeature)
      << "Sorry, but for now only Lite Photon digitization is implemented!\n";

    // Generate dark noise
    if (fDarkNoiseRate > 0.0) AddDarkNoise(opDetWaveforms, nOpChannels);

    // Vary the pedestal
    if (fLineNoise > 0.0) AddLineNoise(opDetWaveforms);

    for (int channel = 0; channel != nOpChannels; ++channel)
      if (fRecordMap[channel])
      {
        // Produce ADC pulse of integers rather than floats
        raw::OpDetWaveform adcVec(0.0, channel, nSamples);

        // Round each value, then cast to short
        for (float value : opDetWaveforms[channel])
          adcVec.emplace_back(short(std::roundf(value)));

        pulseVecPtr->emplace_back(std::move(adcVec));

      }

    // Push the OpDetWaveforms into the event
    evt.put(std::move(pulseVecPtr));

    // Empty the map
    fRecordMap.clear();

  }

  //---------------------------------------------------------------------------
  void OpDetDigitizerDUNE::AddPulse(int timeBin, 
                                    int scale, 
                                    std::vector< float >& waveform)
  {

    // How many bins will be changed
    size_t pulseLength = fSinglePEWaveform.size();
    if ((timeBin + fSinglePEWaveform.size()) > waveform.size()) 
      pulseLength = (waveform.size() - timeBin);

    // Adding a pulse to the waveform
    for (size_t tick = 0; tick != pulseLength; ++tick)
      waveform[timeBin + tick] += scale*fSinglePEWaveform[tick];

  }

  //---------------------------------------------------------------------------
  float OpDetDigitizerDUNE::Pulse1PE(float time) const
  {

    if (time < fPeakTime) return 
      (fVoltageToADC*fMaxAmplitude*std::exp((time - fPeakTime)/fFrontTime));
    else return 
      (fVoltageToADC*fMaxAmplitude*std::exp(-(time - fPeakTime)/fBackTime));

  }

  //---------------------------------------------------------------------------
  void OpDetDigitizerDUNE::CreateSinglePEWaveform()
  {

    size_t length = size_t(fPulseLength*fSampleFreq + 0.5);
    fSinglePEWaveform.resize(length);
    for (size_t tick = 0; tick != length; ++tick)
      fSinglePEWaveform[tick] = Pulse1PE(float(tick)/fSampleFreq);

  }

  //---------------------------------------------------------------------------
  void OpDetDigitizerDUNE::CreateOpDetWaveform
                             (sim::SimPhotonsLite const& opDet,
                              opdet::OpDetResponseInterface const& odResponse,
                              std::vector< std::vector< float > >& 
                                                                opDetWaveforms)
  {

    int const channel = opDet.OpChannel;
    int readoutCh;
    // For a group of photons arriving at the same time this is a map
    // of < arrival time (in ns), number of photons >
    std::map< int, int > const& photonsMap = opDet.DetectedPhotons;

    // For every pair of (arrival time, number of photons) in the map:
    for (auto const& pulse : photonsMap)
    {
      // Converting ns to us
      float photonTime = float(pulse.first/1000.0);
      for (int i = 0; i < pulse.second; ++i)
      {
        if ((photonTime >= fTimeBegin) && (photonTime < fTimeEnd)) {
          // Sample a random subset according to QE
          if (odResponse.detectedLite(channel, readoutCh))
          {
            // Convert the time of the pulse to ticks
            int timeBin = int((photonTime - fTimeBegin)*fSampleFreq);
            // Add 1 pulse to the waveform
            AddPulse(timeBin, CrossTalk(), opDetWaveforms[readoutCh]);
            // Make sure this waveform is recorded
            fRecordMap[readoutCh] = true;
          }
        }
        else 
          mf::LogInfo("OpDetDigitizerDUNE") 
            << "Throwing away an out-of-time photon at " << photonTime << "\n";
      }
    }

  }


  //---------------------------------------------------------------------------
  void OpDetDigitizerDUNE::AddLineNoise
                           (std::vector< std::vector< float > >& waveforms)
  {

    for(auto& waveform : waveforms)
      for(float& value : waveform) value += 
                                        float(fRandGauss->fire(0, fLineNoise));

  }

  //---------------------------------------------------------------------------
  void OpDetDigitizerDUNE::AddDarkNoise
                           (std::vector< std::vector< float > >& waveforms,
                            int nChannels)
  {

    for (int channel = 0; channel < nChannels; ++channel)
    {
      // Multiply by 10^6 since fDarkNoiseRate is in Hz
      float darkNoiseTime = 
        float(fRandExponential->fire(1.0/fDarkNoiseRate)*1000000.0) 
        + fTimeBegin;
      while (darkNoiseTime < fTimeEnd) 
      {
        int timeBin = int((darkNoiseTime - fTimeBegin)*fSampleFreq);
        AddPulse(timeBin, CrossTalk(), waveforms[channel]);
        darkNoiseTime += 
          float(fRandExponential->fire(1.0/fDarkNoiseRate)*1000000.0);
        // Make sure this waveform is recorded
        fRecordMap[channel] = true;
      }

    }
  }
  
  //---------------------------------------------------------------------------
  int OpDetDigitizerDUNE::CrossTalk() const
  {
    if      (fCrossTalk <= 0.0)                 return 1;
    else if (fRandFlat->fire(1.0) > fCrossTalk) return 1;
    else                                        return 2;
  }

}
