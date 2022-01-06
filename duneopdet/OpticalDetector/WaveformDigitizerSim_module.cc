//=========================================================
// WaveformDigitizerSim_module.cc
// This module produces digitized waveforms (creating OpDetWaveform)
// from photon detectors taking OpDetDivRec as input.
// The only random numbers thrown are for line noise.
//
// Gleb Sinev, Duke, 2015
// Anne Christensen, CSU, 2019
// Alex Himmel, FNAL, 2021
//=========================================================

#ifndef WaveformDigitizerSim_h
#define WaveformDigitizerSim_h

// Framework includes

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "canvas/Utilities/Exception.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/DelegatedParameter.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Sequence.h"


// ART extensions
#include "nurandom/RandomUtils/NuRandomService.h"

// LArSoft includes

#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/Simulation/OpDetBacktrackerRecord.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larana/OpticalDetector/OpDetResponseInterface.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoSiPM.h"
#include "dune/OpticalDetector/AlgoSSPLeadingEdge.h"
#include "dune/DuneObj/OpDetDivRec.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"

// CLHEP includes

#include "CLHEP/Random/RandGauss.h"

// C++ includes

#include <vector>
#include <map>
#include <cmath>
#include <memory>
#include <limits>
#include <iomanip>

// ROOT includes

#include "TTree.h"


namespace opdet {

  using std::vector;
  using std::pair;
  typedef std::vector< std::pair< size_t, size_t > > Ranges_t;

  class FocusList
  {
  public:
    FocusList(size_t nSamples, size_t padding)
      : fNSamples(nSamples), fPadding(padding) 
    {}

    void AddRange(size_t from, size_t to)
    {
      from -= std::min(from, fPadding);
      to   =  std::min(to+fPadding, fNSamples-1);

      for(size_t i = 0; i < ranges.size(); ++i){
        pair<size_t, size_t>& r = ranges[i];
        // Completely nested, discard
        if(from >= r.first && to <= r.second) return;
        // Extend end
        if(from >= r.first && from <= r.second){
          r.second = to;
          return;
        }
        // Extend front
        if(to >= r.first && to <= r.second){
          r.first = from;
          return;
        }
      }
      // Discontiguous, add
      ranges.emplace_back(from, to);
    }

    void everything()
    {
      ranges.clear();
      ranges.emplace_back(0, fNSamples-1);
    }

    
    Ranges_t ranges;

  protected:
    size_t fNSamples;
    size_t fPadding;
  };

  class WaveformDigitizerSim : public art::EDProducer{

  public:

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      // Inputs
      fhicl::Sequence<art::InputTag> InputTags           { Name("InputTags"), Comment("Input tags for OpDetDivRecs")  };

      // Single photoelectron pulse parameters
      fhicl::Atom<double>            VoltageToADC        { Name("VoltageToADC"),
                                                           Comment("mV per ADC count") };
      fhicl::Atom<double>            PulseLength         { Name("PulseLength"),         
                                                           Comment("1PE pulse length in us") };
      fhicl::Atom<double>            PeakTime            { Name("PeakTime"),
                                                           Comment("Time when the pulse reaches its maximum in us") };
      fhicl::Atom<double>            MaxAmplitude        { Name("MaxAmplitude"), 
                                                           Comment("Maximum amplitude of the pulse in mV") };
      fhicl::Atom<double>            FrontTime           { Name("FrontTime"), 
                                                           Comment("Constant in the exponential function of the leading edge in us") };
      fhicl::Atom<double>            BackTime            { Name("BackTime"), 
                                                           Comment("Constant in the exponential function of the tail in us") };

      // Waveform output Settings
      fhicl::Atom<size_t>            Padding             { Name("Padding"),
                                                           Comment("Minimum ticks around pulses to simulate")};
      fhicl::Atom<size_t>            ReadoutWindow       { Name("ReadoutWindow"),
                                                           Comment("Total size of the triggered readout window in ticks")};
      fhicl::Atom<size_t>            PreTrigger          { Name("PreTrigger"),
                                                           Comment("Length of the ReadoutWindow which comes before the trigger in ticks")};
      fhicl::Atom<double>            Threshold           { Name("Threshold"),
                                                           Comment("Minimum threshold to trigger readout in the CFD (PE)")};
      fhicl::Atom<size_t>            Dwindow             { Name("Dwindow"),
                                                           Comment("Sample time difference in the CFD trigger in ticks")};
      
      // Digitizer properties
      fhicl::Atom<short>             Pedestal            { Name("Pedestal"),
                                                           Comment("Baseline pedestal in ADC counts")};
      fhicl::Atom<double>            LineNoiseRMS        { Name("LineNoiseRMS"),
                                                           Comment("Pedestal RMS in ADC counts")};
      fhicl::OptionalAtom<short>     DynamicBitRange     { Name("DynamicBitRange"),
                                                           Comment("Maximum number of bits in readout if specified.") };

      // Optional debugging settings
      fhicl::OptionalAtom<double>    TimeBegin           { Name("TimeBegin"),
                                                           Comment("Override earliest allowed waveform time, default -1 drift window") };
      fhicl::OptionalAtom<double>    TimeEnd             { Name("TimeEnd"), 
                                                           Comment("Override latest allowed waveform time, default end of TPC readout") };
      fhicl::Atom<bool>              FullWaveformOutput  { Name("FullWaveformOutput"),  
                                                           Comment("Write out the whole waveform, slow with *large* output sizes. Default false."), 
                                                           false };
    };
    using Parameters = art::EDProducer::Table<Config>;

    explicit WaveformDigitizerSim(Parameters const & config);
    void produce(art::Event&);

  private:

    //////////////////////
    // FHICL Parameters //
    //////////////////////

    vector<art::InputTag> fInputTags;
 
    // Single photoelectron pulse parameters
    double  fVoltageToADC;
    double  fPulseLength;
    double  fPeakTime;
    double  fMaxAmplitude;
    double  fFrontTime;
    double  fBackTime;

    // Waveform output settings
    size_t  fPadding;
    size_t  fReadoutWindow;
    size_t  fPreTrigger;
    double  fThresholdPE;
    size_t  fDwindow;

    // Derived from the above

    // Digitizer properties
    short   fPedestal;
    double  fLineNoiseRMS;
    int     fMaxSaturationCutOff; // derived from above

    // Optional debugging settings
    double  fTimeBegin;
    double  fTimeEnd;
    bool    fFullWaveformOutput;


    ////////////////////////////////////
    // Other Members set during setup //
    ////////////////////////////////////
    
    // Random number engines
    CLHEP::HepRandomEngine& fOpDigiEngine;
    CLHEP::RandGauss        fRandGauss;

    double fSampleFreqMHz; 
    size_t fPulseLengthTicks;
    double fThresholdADC;

    // Template for a single PE
    vector< double > fSinglePEWaveform;


    /////////////////////
    // Setup functions //
    /////////////////////

    void   CreateSinglePEWaveform();
    double Pulse1PE(double time_in_us) const;
    void   SetDefaultBeginEndTimes();
    void   CheckFHiCLParameters() const;


    ///////////////////////////////////
    // Support functions for produce //
    ///////////////////////////////////

    // Add photons to an existing waveform

    void AddPEsToWaveform(const sim::OpDetDivRec* dr_p,
                          vector<double> &        pdWaveform,
                          FocusList &             fls) const;

    // Vary the pedestal
    void AddLineNoise(vector< double > & waveform, 
                      const FocusList  & fls);

    // Apply saturation and cast into shorts
    // OpDetWaveform constructor requires this specific type 
    // (probably not for good reasons)
    vector< uint16_t > Digitize(vector<double>::iterator itBegin, 
                                vector<double>::iterator itEnd) const;

    // Constant fraction discriminator single channel trigger
    template<typename T> Ranges_t CFDTrigger(vector<T> const& wf, const FocusList& fls) const;


    // Functions to convert between time and ticks, being explicit about units
    // Times within the event count from fBeginTime, Ticks count from 0
    double  Tick2us(size_t tick) const       { return fTimeBegin + static_cast<double>(tick)/fSampleFreqMHz; };
    double  Tick2ns(size_t tick) const       { return 1000.*Tick2us(tick); };
    size_t  us2Tick(double time_in_us) const { return static_cast<size_t>(std::round( (time_in_us-fTimeBegin)*fSampleFreqMHz) ); }
    size_t  ns2Tick(double time_in_ns) const { return us2Tick(time_in_ns/1000.); };

  };

}

#endif

namespace opdet {

  DEFINE_ART_MODULE(WaveformDigitizerSim)

}

namespace opdet {

  /////////////////////////////////////
  // Constructor and setup functions //
  /////////////////////////////////////

  //---------------------------------------------------------------------------
  WaveformDigitizerSim::WaveformDigitizerSim(Parameters const & config)
    : EDProducer{config}
    , fInputTags{           config().InputTags() }

    , fVoltageToADC{        config().VoltageToADC() }
    , fPulseLength{         config().PulseLength() }
    , fPeakTime{            config().PeakTime() }
    , fMaxAmplitude{        config().MaxAmplitude() }
    , fFrontTime{           config().FrontTime() }
    , fBackTime{            config().BackTime() }

    , fPadding{             config().Padding() }
    , fReadoutWindow{       config().ReadoutWindow() }
    , fPreTrigger{          config().PreTrigger() }
    , fThresholdPE{         config().Threshold() }
    , fDwindow{             config().Dwindow() }
    
    , fPedestal{            config().Pedestal() }
    , fLineNoiseRMS{        config().LineNoiseRMS() }
    , fMaxSaturationCutOff{ std::numeric_limits<int>::max() }

    , fFullWaveformOutput{  config().FullWaveformOutput() }

    , fOpDigiEngine( art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, 
                                                                               "HepJamesRandom",
                                                                               "waveformdigi", 
                                                                               config.get_PSet(), 
                                                                               "SeedWaveformDigi") )
    , fRandGauss(fOpDigiEngine)
  {

    // This module produces (infrastructure piece)
    produces< vector< raw::OpDetWaveform > >();

    for (auto tag: fInputTags) {
      consumes< vector< sim::OpDetDivRec > >(tag);
    }

    // Get the optical clock frequency
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
    fSampleFreqMHz = clockData.OpticalClock().Frequency();
    mf::LogInfo("WaveformDigitizerSim") << "Using a sampling frequency of " << fSampleFreqMHz << " MHz";
    

    // Creating a single photoelectron waveform template
    // based on fhicl configuration
    CreateSinglePEWaveform();


    // Set a dynamic range if given
    short br;
    if ( config().DynamicBitRange(br) ) {
      fMaxSaturationCutOff = pow(2, br) - 1;
      mf::LogInfo("WaveformDigitizerSim") << "Limiting output to " << br << " bits";
    }

    
    // Set earliest and latest allowed times
    std::string tsource;
    if (config().TimeBegin(fTimeBegin) && config().TimeEnd(fTimeEnd))
      tsource = "override";
    else {
      tsource = "default";
      SetDefaultBeginEndTimes();
    }
    mf::LogInfo("WaveformDigitizerSim") << "Using " << tsource << " time limits on PD digitizer: "
                                        << fTimeBegin << " us to " << fTimeEnd << " us";

    // Check for valid configuration
    CheckFHiCLParameters();
  }

  //---------------------------------------------------------------------------
  void WaveformDigitizerSim::CreateSinglePEWaveform()
  {
    double maxADC = 0.;
    fPulseLengthTicks = std::round(fPulseLength*fSampleFreqMHz);
    fSinglePEWaveform.resize(fPulseLengthTicks);
    mf::LogInfo("WaveformDigitizerSim") << "Requested pulse length of, " << fPulseLength << " us "
                                        << "which is " << fPulseLengthTicks << " ticks";

    for (size_t tick = 0; tick < fPulseLengthTicks; ++tick) {
      double val = Pulse1PE(static_cast< double >(tick)/fSampleFreqMHz);
      fSinglePEWaveform[tick] = val;
      if (val > maxADC) maxADC = val;
    }

    // Set the ADC threshold based on the PE threshold
    fThresholdADC = fThresholdPE * maxADC;
    mf::LogInfo("WaveformDigitizerSim") << "Requested PE threshold, " << std::fixed
                                        << std::setprecision(2) << fThresholdPE 
                                        << ", converted to ADC threshold " 
                                        << std::setprecision(0) << fThresholdADC;
  }

  //---------------------------------------------------------------------------
  double WaveformDigitizerSim::Pulse1PE(double time_in_us) const
  {

    if (time_in_us < fPeakTime) 
      return (fVoltageToADC*fMaxAmplitude*std::exp((time_in_us - fPeakTime)/fFrontTime));
    else 
      return (fVoltageToADC*fMaxAmplitude*std::exp(-(time_in_us - fPeakTime)/fBackTime));

  }

  //---------------------------------------------------------------------------
  void WaveformDigitizerSim::SetDefaultBeginEndTimes()
  {
    // A little wasteful to get clockData again, but only happens once per job
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
    auto const detProp   = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob(clockData);

    // Take the TPC readout window size and convert to us with the electronics
    // clock frequency, adding 5 us of padding
    fTimeEnd   = detProp.ReadOutWindowSize() / clockData.TPCClock().Frequency() + 5;

    // Assume the readout is symmetrical around 0
    fTimeBegin = -1.*fTimeEnd;

  }


  //---------------------------------------------------------------------------
  void WaveformDigitizerSim::CheckFHiCLParameters() const
  {
    // Check that we have input tags
    if (fInputTags.size() == 0) 
      throw art::Exception(art::errors::Configuration)
        << "No input tags were given to WaveformDigitizerSim.\n";

    // Sanity check the line noise
    if (fLineNoiseRMS < 0.0)
      throw art::Exception(art::errors::Configuration)
        << "fLineNoiseRMS: " << fLineNoiseRMS << '\n'
        << "Line noise RMS should be non-negative!\n";

    // Sanity check beginning and end times
    if (fTimeBegin >= fTimeEnd) {
      throw art::Exception(art::errors::Configuration)
        << "TimeBegin: " << fTimeBegin << " and " << "TimeEnd: "   << fTimeEnd   << '\n'
        << "TimeBegin should be less than TimeEnd!\n";
    }
  }




  /////////////////////////////////////////////////
  // produce data products and support functions //
  /////////////////////////////////////////////////


  //---------------------------------------------------------------------------
  void WaveformDigitizerSim::produce(art::Event& event)
  {
    // A pointer that will store produced OpDetWaveforms
    auto wave_forms_p = std::make_unique< vector< raw::OpDetWaveform > >();

    // Total number of ticks in the whole event
    unsigned int nSamples = (fTimeEnd - fTimeBegin)*fSampleFreqMHz;

    // First, pull all of the DivRec handles from the event, and collect up
    // all DivRecs on the same channel into a single vector
    std::map< int, vector<const sim::OpDetDivRec *> > DivRecsByChannel;
    for (auto tag: fInputTags) {
      auto dr_handle = event.getHandle< vector< sim::OpDetDivRec > >(tag);
      if (!dr_handle) {
        mf::LogWarning("WaveformDigitizerSim") << "Could not load OpDetDivRecs " << tag << ". Skipping.";
        continue;
      }
      for (auto const& dr : *dr_handle) {
        DivRecsByChannel[dr.OpDetNum()].push_back( &dr );
      }
    }

    // Now, loop through channels, treating all photons on a cha
    for (auto const& [opDet, vDivRecs]: DivRecsByChannel)
    {
      // Create the empty waveform vector and focus list
      vector< double > pdWaveform(nSamples, fPedestal);
      FocusList fls(nSamples, fPadding);

      // Add a PE template to the waveform for each true photon
      for (auto dr_p: vDivRecs) AddPEsToWaveform(dr_p, pdWaveform, fls);

      // So that line noise is added to all ticks in full output mode
      if (fFullWaveformOutput)  fls.everything(); 

      // Add line noise
      AddLineNoise(pdWaveform, fls);

      if (fFullWaveformOutput) {
        wave_forms_p->emplace_back(Tick2us(0), opDet, Digitize(pdWaveform.begin(), pdWaveform.end()));
      }
      else {
        // Checking for tiggers on floats, rather than shorts.
        // This is an approximation, but it saves making an extra copy
        // of the waveform and makes the code easier to follow.
        for ( auto t: CFDTrigger(pdWaveform, fls) ) {

          // Digitize and store
          auto shortWF = Digitize(pdWaveform.begin()+t.first, pdWaveform.begin()+t.second+1);
          wave_forms_p->emplace_back(Tick2us(t.first), opDet,  shortWF);
        }
      }
    }

    // Push the OpDetWaveforms into the event
    event.put(std::move(wave_forms_p));

  }

  //---------------------------------------------------------------------------
  void WaveformDigitizerSim::AddPEsToWaveform(sim::OpDetDivRec const* dr_p,
                                              vector<double>&         pdWaveform,
                                              FocusList&              fls) const
  {
    // Vector of DivRec time bins (struct OpDet_Time_Chans)
    for (auto odtc: dr_p->GetTimeChans()) {

      // Extract time for this odtc within the event
      double photonTime_ns = odtc.time;
      size_t timeBin       = ns2Tick(photonTime_ns);

      // Check if the photon is inside the digitization range. If not, skip it.
      if ( timeBin < 0 || timeBin >= pdWaveform.size() ) {
        mf::LogWarning("WaveformDigitizerSim") << "Skipping a photon at " << photonTime_ns/1000. << " us, outside digitization window of " << fTimeBegin << " to " << fTimeEnd;
        continue;
      }

      // Loop through records at this time and count photons
      int nPE = 0;
      for (auto const& sdp : odtc.phots)
        nPE += sdp.phot;

      // Add ticks until the end of the single PE waveform or end of the whole pdWaveform
      size_t stop = std::min(fPulseLengthTicks, pdWaveform.size()-timeBin);

      // Add this range to the focus list
      fls.AddRange(timeBin, timeBin+stop-1);

      // Add the nPE pulse to the waveform
      for (size_t tick = 0; tick < stop; ++tick)
        pdWaveform[timeBin+tick] += fSinglePEWaveform[tick]*nPE;
    }
  }

  //---------------------------------------------------------------------------
  void WaveformDigitizerSim::AddLineNoise(vector< double >& waveform, const FocusList& fls) 
  {
    if (fLineNoiseRMS <= 0.0) return;

    for (auto p: fls.ranges) {
      for(auto k = p.first; k <= p.second; ++k){
        waveform[k] += fRandGauss.fire(0, fLineNoiseRMS);
      }
    }
  }

  //---------------------------------------------------------------------------
  vector< uint16_t > WaveformDigitizerSim::Digitize(vector<double>::iterator itBegin, vector<double>::iterator itEnd) const
  {
    for(auto it = itBegin; it != itEnd; ++it) {
      if(*it > fMaxSaturationCutOff) 
        *it = fMaxSaturationCutOff; 
    }

    // Don't bother to round properly, it's faster this way
    return vector< uint16_t >(itBegin, itEnd);
  }

  //---------------------------------------------------------------------------
  template<typename T> 
  Ranges_t WaveformDigitizerSim::CFDTrigger(vector<T> const& wf, const FocusList& fls) const
  {

    Ranges_t readouts;

    for (auto range: fls.ranges) {
      size_t  wstart = -1;
      size_t  wend   = -1;
      bool fire   = false;

      for (size_t tick = range.first; tick < range.second - fDwindow; ++tick) {

        // Fire CFD
        if (wf[tick+fDwindow] - wf[tick] > fThresholdADC) {

          if (!fire) {
            // Start a new readout window 
            fire   = true;
            wstart = tick-fPreTrigger;
            wend   = tick-fPreTrigger+fReadoutWindow;
          }
          else {
            // Extend the current readout window
            // Simplest to implement here, update to the
            // actual DAPHNE algorithm once known
            wend = tick-fPreTrigger+fReadoutWindow;
          }

        }
        else if (fire && tick >= wend) {
          // We've reached the end of a window, save it
          fire = false;
          readouts.emplace_back(wstart, wend);
        }

      }

      // Check for lingering window, add a final window
      // up to the end of the waveform if so.
      if (fire == true) {
        readouts.emplace_back(wstart, range.second-1);
      }
    }

    return readouts;
  }



} // end namespace

