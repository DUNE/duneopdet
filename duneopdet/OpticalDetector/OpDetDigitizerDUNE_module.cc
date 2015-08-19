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
#include "OpticalDetector/AlgoSiPM.h"

// CLHEP includes

#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandExponential.h"
#include "CLHEP/Random/RandFlat.h"

// C++ includes

#include <vector>
#include <map>
#include <cmath>
#include <memory>

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

      // Threshold algorithm
      std::unique_ptr< pmtana::AlgoSiPM > fThreshAlg;

      // Random number engines
      CLHEP::RandGauss       *fRandGauss;
      CLHEP::RandExponential *fRandExponential;
      CLHEP::RandFlat        *fRandFlat;

      // Function that adds n pulses to a waveform
      void AddPulse(unsigned int timeBin, int scale, 
                    std::vector< float >& waveform);

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

      // Produce waveform on one of the optical detectors
      void CreatePDWaveform(sim::SimPhotonsLite const&, 
                            opdet::OpDetResponseInterface const&,
                            geo::Geometry const&,
                            std::vector< std::vector< float > >&);

      // Vary the pedestal
      void AddLineNoise(std::vector< std::vector< float > >&);

      void AddDarkNoise(std::vector< std::vector< float > >&);

      unsigned short CrossTalk() const;

      // Create a vector of shorts from a vector of doubles
      // rounding it properly
      std::vector< short > VectorOfFloatsToVectorOfShorts
                                           (std::vector< float > const&);

      // Make several shorter waveforms out of a long one using a hit finder, 
      // recording also when they start in the long waveform
      std::map< size_t, std::vector< short > > 
                              SplitWaveform(std::vector< short > const&);
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

    fThreshAlg = std::make_unique< pmtana::AlgoSiPM >
                    (pset.get< fhicl::ParameterSet >("algo_threshold"));

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

    if (!fUseLitePhotons)
      throw art::Exception(art::errors::UnimplementedFeature)
        << "Sorry, but for now only Lite Photon digitization is implemented!"
        << '\n';

    // Total number of ticks in our readout
    unsigned int nSamples = (fTimeEnd - fTimeBegin)*fSampleFreq;

    // Geometry service
    art::ServiceHandle< geo::Geometry > geometry;

    // Service for determining optical detector responses
    art::ServiceHandle< opdet::OpDetResponseInterface > odResponse;

    // Get SimPhotonsLite from the event
    art::Handle< std::vector< sim::SimPhotonsLite > > litePhotonHandle;
    evt.getByLabel(fInputModule, litePhotonHandle);

    // For every optical detector:
    for (auto const& litePhotons : (*litePhotonHandle)) 
    {
      // OpChannel in SimPhotonsLite is actually the photon detector number
      unsigned int opDet = litePhotons.OpChannel;

      // Get number of channels in this optical detector
      unsigned int nChannelsPerOpDet = geometry->NOpHardwareChannels(opDet);
      // This vector stores waveforms created for each optical channel
      std::vector< std::vector< float > > pdWaveforms(nChannelsPerOpDet, 
        std::vector< float >(nSamples, static_cast< float >(fPedestal)));

      CreatePDWaveform(litePhotons, *odResponse, *geometry, pdWaveforms);

      // Generate dark noise
      if (fDarkNoiseRate > 0.0) AddDarkNoise(pdWaveforms);

      // Vary the pedestal
      if (fLineNoise > 0.0) AddLineNoise(pdWaveforms);

      // Loop over all the created waveforms, split them into shorter
      // waveforms and use them to initialize OpDetWaveforms
      for (unsigned int hardwareChannel = 0; 
           hardwareChannel < nChannelsPerOpDet; ++hardwareChannel)
      {
        std::vector< short > waveformOfShorts = 
             VectorOfFloatsToVectorOfShorts(pdWaveforms.at(hardwareChannel));

        std::map< size_t, std::vector < short > > mapTickWaveform = 
                                             SplitWaveform(waveformOfShorts);

        unsigned int opChannel = geometry->OpChannel(opDet, hardwareChannel);

        for (auto const& pairTickWaveform : mapTickWaveform)
        {
          double timeStamp = 
            static_cast< double >(pairTickWaveform.first/fSampleFreq) 
                                                               + fTimeBegin;

          raw::OpDetWaveform adcVec(timeStamp, opChannel, 
                                    pairTickWaveform.second.size());

          for (short const& value : pairTickWaveform.second)
            adcVec.emplace_back(value);

          pulseVecPtr->emplace_back(std::move(adcVec));
        }
      }
    }

    // Push the OpDetWaveforms into the event
    evt.put(std::move(pulseVecPtr));

  }

  //---------------------------------------------------------------------------
  void OpDetDigitizerDUNE::AddPulse(unsigned int timeBin, 
                                    int scale, std::vector< float >& waveform)
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
      fSinglePEWaveform[tick] = 
         Pulse1PE(static_cast< float >(tick)/fSampleFreq);

  }

  //---------------------------------------------------------------------------
  void OpDetDigitizerDUNE::CreatePDWaveform
                             (sim::SimPhotonsLite const& litePhotons,
                              opdet::OpDetResponseInterface const& odResponse,
                              geo::Geometry const& geometry,
                              std::vector< std::vector< float > >& 
                                                                pdWaveforms)
  {

    unsigned int const opDet = litePhotons.OpChannel;
    // This is int because otherwise detectedLite doesn't work
    int readoutChannel;
    // For a group of photons arriving at the same time this is a map
    // of < arrival time (in ns), number of photons >
    std::map< int, int > const& photonsMap = litePhotons.DetectedPhotons;

    // For every pair of (arrival time, number of photons) in the map:
    for (auto const& pulse : photonsMap)
    {
      // Converting ns to us
      float photonTime = static_cast< float >(pulse.first/1000.0);
      for (int i = 0; i < pulse.second; ++i)
      {
        if ((photonTime >= fTimeBegin) && (photonTime < fTimeEnd)) 
        {
          // Sample a random subset according to QE
          if (odResponse.detectedLite(opDet, readoutChannel))
          {
            unsigned int hardwareChannel = 
                      geometry.HardwareChannelFromOpChannel(readoutChannel);
            // Convert the time of the pulse to ticks
            unsigned int timeBin = static_cast< unsigned int >
                                     ((photonTime - fTimeBegin)*fSampleFreq);
            // Add 1 pulse to the waveform
            AddPulse(timeBin, CrossTalk(), pdWaveforms.at(hardwareChannel));
          }
        }
        else 
          mf::LogInfo("OpDetDigitizerDUNE") 
            << "Throwing away an out-of-time photon at " << photonTime << '\n';
      }
    }

  }


  //---------------------------------------------------------------------------
  void OpDetDigitizerDUNE::AddLineNoise
                           (std::vector< std::vector< float > >& waveforms)
  {

    for(auto& waveform : waveforms)
      for(float& value : waveform) value += 
                         static_cast< float >(fRandGauss->fire(0, fLineNoise));

  }

  //---------------------------------------------------------------------------
  void OpDetDigitizerDUNE::AddDarkNoise
                           (std::vector< std::vector< float > >& waveforms)
  {

    for (auto& waveform : waveforms)
    {
      // Multiply by 10^6 since fDarkNoiseRate is in Hz
      float darkNoiseTime = static_cast< float >(fRandExponential->
                            fire(1.0/fDarkNoiseRate)*1000000.0) + fTimeBegin;
      while (darkNoiseTime < fTimeEnd) 
      {
        unsigned int timeBin = static_cast< unsigned int >
                               ((darkNoiseTime - fTimeBegin)*fSampleFreq);
        AddPulse(timeBin, CrossTalk(), waveform);
        darkNoiseTime += static_cast< float >
                        (fRandExponential->fire(1.0/fDarkNoiseRate)*1000000.0);
      }

    }
  }
  
  //---------------------------------------------------------------------------
  unsigned short OpDetDigitizerDUNE::CrossTalk() const
  {
    // Sometimes this should produce 3 or more PEs (not implemented)
    if      (fCrossTalk <= 0.0)                 return 1;
    else if (fRandFlat->fire(1.0) > fCrossTalk) return 1;
    else                                        return 2;
  }

  //---------------------------------------------------------------------------
  std::vector< short > OpDetDigitizerDUNE::VectorOfFloatsToVectorOfShorts(
                                    std::vector< float > const& vectorOfFloats)
  {

    std::vector< short > vectorOfShorts;
    vectorOfShorts.reserve(vectorOfFloats.size());

    for (short const& value : vectorOfFloats)
      vectorOfShorts.emplace_back(static_cast< short >(std::roundf(value)));

    return vectorOfShorts;

  }

  //---------------------------------------------------------------------------
  std::map< size_t, std::vector< short > > OpDetDigitizerDUNE::SplitWaveform(
                                          std::vector< short > const& waveform)
  {

    std::map< size_t, std::vector< short > > mapTickWaveform;

    fThreshAlg->RecoPulse(waveform);

    std::vector< pmtana::pulse_param > pulses;
    for (size_t pulseCounter = 0; pulseCounter < fThreshAlg->GetNPulse(); 
                                                          ++pulseCounter)
      pulses.emplace_back(fThreshAlg->GetPulse(pulseCounter));

    // We have to refine this algorithm later
    for (pmtana::pulse_param const& pulse : pulses)
    {
      if (pulse.t_end <= pulse.t_start)
        // Can I call it a logic error?
        throw art::Exception(art::errors::LogicError)
          << "Pulse ends before or at the same time it starts!\n";

      std::vector< short >::const_iterator window_start = 
              waveform.begin() + static_cast< size_t >(pulse.t_start);
      std::vector< short >::const_iterator window_end   = 
              waveform.begin() + static_cast< size_t >(pulse.t_end  );
      mapTickWaveform[static_cast< size_t >(pulse.t_start)] = 
                       std::vector< short >(window_start, window_end);
    }


    return mapTickWaveform;

  }

}
