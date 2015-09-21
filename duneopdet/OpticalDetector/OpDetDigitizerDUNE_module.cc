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
#include "Utilities/LArProperties.h"
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
      std::string fInputModule;   // Input tag for OpDet collection
      float  fSampleFreq;         // Sampling frequency in MHz
      float  fTimeBegin;          // Beginning of waveform in us
      float  fTimeEnd;            // End of waveform in us
      float  fVoltageToADC;       // Conversion factor mV to ADC counts
      float  fLineNoise;          // Pedestal RMS in ADC counts
      float  fDarkNoiseRate;      // In Hz
      float  fCrossTalk;          // Probability of SiPM producing 2 PE signal
                                  // in response to 1 photon
      short  fPedestal;           // In ADC counts
      bool   fDefaultSimWindow;   // Set the start time to -1 drift window and
                                  // the end time to the end time 
                                  // of the TPC readout
      bool   fFullWaveformOutput; // Output full waveforms -- produces large
                                  // output. Mostly for debug purposes
      size_t fReadoutWindow;      // In ticks
      size_t fPreTrigger;         // In ticks

      // Threshold algorithm
      std::unique_ptr< pmtana::AlgoSiPM > fThreshAlg;

      // Random number engines
      std::unique_ptr< CLHEP::RandGauss >       fRandGauss;
      std::unique_ptr< CLHEP::RandExponential > fRandExponential;
      std::unique_ptr< CLHEP::RandFlat >        fRandFlat;

      // Function that adds n pulses to a waveform
      void AddPulse(size_t timeBin, int scale, 
                    std::vector< float >& waveform) const;

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
                            std::vector< std::vector< float > >&) const;

      // Vary the pedestal
      void AddLineNoise(std::vector< std::vector< float > >&) const;

      void AddDarkNoise(std::vector< std::vector< float > >&) const;

      unsigned short CrossTalk() const;

      // Create a vector of shorts from a vector of floats
      // rounding it properly
      std::vector< short > VectorOfFloatsToVectorOfShorts
                                           (std::vector< float > const&) const;

      // Make several shorter waveforms out of a long one using a hit finder, 
      // recording also when they start in the long waveform
      std::map< size_t, std::vector< short > > 
                              SplitWaveform(std::vector< short > const&) const;

      float GetDriftWindow() const;

      // Convert time to ticks or the other way around
      // without any checks
      float  TickToTime(size_t tick) const;
      size_t TimeToTick(float  time) const;

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
    fInputModule        = pset.get< std::string >("InputModule"  );
    fVoltageToADC       = pset.get< float  >("VoltageToADC"      );
    fLineNoise          = pset.get< float  >("LineNoise"         );
    fDarkNoiseRate      = pset.get< float  >("DarkNoiseRate"     );
    fCrossTalk          = pset.get< float  >("CrossTalk"         );
    fPedestal           = pset.get< short  >("Pedestal"          );
    fDefaultSimWindow   = pset.get< bool   >("DefaultSimWindow"  );
    fFullWaveformOutput = pset.get< bool   >("FullWaveformOutput");
    fReadoutWindow      = pset.get< size_t >("ReadoutWindow"     );
    fPreTrigger         = pset.get< size_t >("PreTrigger"        );

    fThreshAlg = std::make_unique< pmtana::AlgoSiPM >
                   (pset.get< fhicl::ParameterSet >("algo_threshold"));

    // Obtaining parameters from the TimeService
    art::ServiceHandle< util::TimeService > timeService;
    fSampleFreq = timeService->OpticalClock().Frequency();

    if (fDefaultSimWindow)
    {
      // Assume the readout starts at -1 drift window
      fTimeBegin = -1*GetDriftWindow();

      // Take the TPC readout window size and convert 
      // to us with the electronics clock frequency
      fTimeEnd   = art::ServiceHandle< util::DetectorProperties >()->
                   ReadOutWindowSize()/timeService->TPCClock().Frequency();
    }
    else
    {
      fTimeBegin = pset.get< float >("TimeBegin");
      fTimeEnd   = pset.get< float >("TimeEnd"  );
    }

    // Initializing random number engines
    unsigned int seed = 
             pset.get< unsigned int >("Seed", sim::GetRandomNumberSeed());
    createEngine(seed);

    art::ServiceHandle< art::RandomNumberGenerator > rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine();
    fRandGauss       = std::make_unique< CLHEP::RandGauss       >(engine);
    fRandExponential = std::make_unique< CLHEP::RandExponential >(engine);
    fRandFlat        = std::make_unique< CLHEP::RandFlat        >(engine);

    // Creating a single photoelectron waveform
    // Hardcoded, probably need to read them from the FHiCL file
    fPulseLength  = 4.0;
    fPeakTime     = 0.260;
    fMaxAmplitude = 0.12;
    fFrontTime    = 0.009;
    fBackTime     = 0.476;
    CreateSinglePEWaveform();

  }
   
  //---------------------------------------------------------------------------
  void OpDetDigitizerDUNE::produce(art::Event& evt)
  {
    
    // A pointer that will store produced OpDetWaveforms
    std::unique_ptr< std::vector< raw::OpDetWaveform > > 
      pulseVecPtr(std::make_unique< std::vector< raw::OpDetWaveform > >());
    
    art::ServiceHandle< sim::LArG4Parameters > lgp;
    bool fUseLitePhotons = lgp->UseLitePhotons();

    if (!fUseLitePhotons)
      throw art::Exception(art::errors::UnimplementedFeature)
        << "Sorry, but for now only Lite Photon digitization is implemented!"
        << '\n';

    // Total number of ticks in the full waveform
    // Including one pretrigger window before the waveform
    // and one readout window - pretrigger window after the waveform
    unsigned int nSamples = (fTimeEnd - fTimeBegin)*fSampleFreq
                          + fReadoutWindow;

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
      if (fLineNoise > 0.0)     AddLineNoise(pdWaveforms);

      // Loop over all the created waveforms, split them into shorter
      // waveforms and use them to initialize OpDetWaveforms
      for (unsigned int hardwareChannel = 0; 
           hardwareChannel < nChannelsPerOpDet; ++hardwareChannel)
      {
        std::vector< short > waveformOfShorts = 
             VectorOfFloatsToVectorOfShorts(pdWaveforms.at(hardwareChannel));

        std::map< size_t, std::vector < short > > mapTickWaveform = 
          (!fFullWaveformOutput) ? 
          SplitWaveform(waveformOfShorts) : 
          std::map< size_t, std::vector< short > >{ std::make_pair(0, 
                                                          waveformOfShorts) };

        unsigned int opChannel = geometry->OpChannel(opDet, hardwareChannel);

        for (auto const& pairTickWaveform : mapTickWaveform)
        {
          double timeStamp = 
            static_cast< double >(TickToTime(pairTickWaveform.first));

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
  void OpDetDigitizerDUNE::AddPulse(size_t timeBin, 
                              int scale, std::vector< float >& waveform) const
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

    size_t length = 
      static_cast< size_t > (std::round(fPulseLength*fSampleFreq));
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
                                                            pdWaveforms) const
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
            size_t timeBin = TimeToTick(photonTime);
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
                        (std::vector< std::vector< float > >& waveforms) const
  {

    for(auto& waveform : waveforms)
      for(float& value : waveform) value += 
                         static_cast< float >(fRandGauss->fire(0, fLineNoise));

  }

  //---------------------------------------------------------------------------
  void OpDetDigitizerDUNE::AddDarkNoise
                        (std::vector< std::vector< float > >& waveforms) const
  {

    for (auto& waveform : waveforms)
    {
      // Multiply by 10^6 since fDarkNoiseRate is in Hz
      float darkNoiseTime = static_cast< float >(fRandExponential->
                            fire(1.0/fDarkNoiseRate)*1000000.0) + fTimeBegin;
      while (darkNoiseTime < fTimeEnd) 
      {
        size_t timeBin = TimeToTick(darkNoiseTime);
        AddPulse(timeBin, CrossTalk(), waveform);
        // Find next time to simulate a single PE pulse
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
  std::vector< short > OpDetDigitizerDUNE::VectorOfFloatsToVectorOfShorts
                            (std::vector< float > const& vectorOfFloats) const
  {

    std::vector< short > vectorOfShorts;
    vectorOfShorts.reserve(vectorOfFloats.size());

    for (short const& value : vectorOfFloats)
      vectorOfShorts.emplace_back(static_cast< short >(std::round(value)));

    return vectorOfShorts;

  }

  //---------------------------------------------------------------------------
  std::map< size_t, std::vector< short > > OpDetDigitizerDUNE::SplitWaveform
                                  (std::vector< short > const& waveform) const
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
      mapTickWaveform.emplace(static_cast< size_t >(pulse.t_start), 
                              std::vector< short >(window_start, window_end));
      // Don't forget to check that the time output by the (new) algortihm
      // is < fTimeEnd and > fTimeBegin!
    }

    return mapTickWaveform;

  }

  //---------------------------------------------------------------------------
  float OpDetDigitizerDUNE::GetDriftWindow() const
  {

    float driftWindow;

    float maxDrift = 0.0;
    for (geo::TPCGeo const& tpc : 
           art::ServiceHandle< geo::Geometry >()->IterateTPCs())
      if (maxDrift < tpc.DriftDistance()) maxDrift = tpc.DriftDistance();

    driftWindow = 
      maxDrift/art::ServiceHandle< util::LArProperties >()->DriftVelocity();

    return driftWindow;

  }

  //---------------------------------------------------------------------------
  float OpDetDigitizerDUNE::TickToTime(size_t tick) const
  {

    if (tick > fPreTrigger)
      return (static_cast< float >(tick - fPreTrigger)/fSampleFreq 
                                                                 + fTimeBegin);
    else
      return (static_cast< float >(fPreTrigger - tick)/fSampleFreq*(-1.0)
                                                                 + fTimeBegin);

  }

  //---------------------------------------------------------------------------
  size_t OpDetDigitizerDUNE::TimeToTick(float time) const
  {

    return static_cast< size_t >(std::round((time - fTimeBegin)*fSampleFreq 
                                                               + fPreTrigger));

  }
  
}
