//=========================================================
// OpDetDigitizerDUNE_module.cc
// This module produces digitized output
// (creating OpDetWaveform)
// from photon detectors taking SimPhotonsLite as input.
//
// Gleb Sinev, Duke, 2015
// Maritza Delgado, 2022 SPE testbench option
// Based on OpMCDigi_module.cc
//=========================================================

#ifndef OpDetDigitizerDUNE_h
#define OpDetDigitizerDUNE_h

// Framework includes

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/Exception.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"



// ART extensions
#include "nurandom/RandomUtils/NuRandomService.h"

// LArSoft includes

#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/Simulation/OpDetBacktrackerRecord.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larana/OpticalDetector/OpDetResponseInterface.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoSiPM.h"
#include "duneopdet/OpticalDetector/AlgoSSPLeadingEdge.h"
#include "dunecore/DuneObj/OpDetDivRec.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"

// CLHEP includes

#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandExponential.h"
#include "CLHEP/Random/RandFlat.h"

// C++ includes

#include <vector>
#include <map>
#include <cmath>
#include <memory>
#include <string>
#include <sstream>
#include <fstream>

// ROOT includes

#include "TTree.h"
#include "TFile.h"


namespace opdet {

  class FocusList
  {
  public:
      FocusList(int nSamples, int padding)
        : fNSamples(nSamples), fPadding(padding) {}

      void AddRange(int from, int to)
      {
        from -= fPadding;
        to += fPadding;

        if(from < 0) from = 0;
        if(to >= fNSamples) to = fNSamples-1;

        for(unsigned int i = 0; i < ranges.size(); ++i){
          std::pair<int, int>& r = ranges[i];
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
    
      std::vector<std::pair<int, int>> ranges;

    protected:
      int fNSamples;
      int fPadding;
  };

  class OpDetDigitizerDUNE : public art::EDProducer{

  public:

      OpDetDigitizerDUNE(fhicl::ParameterSet const&);
      // Should the destructor be empty?
      //      virtual ~OpDetDigitizerDUNE();
    
      void produce(art::Event&);

    private:

      // The parameters read from the FHiCL file
      std::vector<std::string> vInputModules;
      std::set<std::string> fInputModules;   // Input tag for OpDet collection
      double  fSampleFreq;                   // Sampling frequency in MHz
      double  fTimeBegin;                    // Beginning of waveform in us
      double  fTimeEnd;                      // End of waveform in us
      double  fVoltageToADC;                 // Conversion factor mV to ADC counts
      double  fLineNoiseRMS;                 // Pedestal RMS in ADC counts
      double  fDarkNoiseRate;                // In Hz
      double  fCrossTalk;                    // Probability of SiPM producing 2 PE signal
                                             // in response to 1 photon
      short  fPedestal;                      // In ADC counts
      bool   fDefaultSimWindow;              // Set the start time to -1 drift window and
                                             // the end time to the end time
                                             // of the TPC readout
      bool   fFullWaveformOutput;            // Output full waveforms -- produces large
                                             // output. Mostly for debug purposes
      size_t fReadoutWindow;                 // In ticks
      size_t fPreTrigger;                    // In ticks

      int    fPadding;                       // In ticks
    
      bool   fDigiTree_SSP_LED;              // To create a analysis Tree for SSP LED
      bool   fUseSDPs;                       // = pset.get< bool   >("UseSDPs", true);
      bool   TestbenchSinglePE;                 // File testbench

      double  fQEOverride;
      double  fRefQEOverride;
      std::string fSPEDataFile;

      //-----------------------------------------------------
      // Trigger analysis variables
      std::vector<double> t_photon; // vitor
      std::vector<int>    op_photon;

      TTree *arvore2;
      //-----------------------------------------------------

      // Threshold algorithm
      std::unique_ptr< pmtana::AlgoSSPLeadingEdge > fThreshAlg;

      // Random number engines
      std::unique_ptr< CLHEP::RandGauss       > fRandGauss;
      std::unique_ptr< CLHEP::RandExponential > fRandExponential;
      std::unique_ptr< CLHEP::RandFlat        > fRandFlat;

      // Function that adds n pulses to a waveform
      void AddPulse(size_t timeBin, int scale,
                    std::vector< double >& waveform,
                    FocusList& fl) const;

      // Functional response to one photoelectron (time in ns)
      double Pulse1PE(double time) const;
      // Single photoelectron pulse parameters
      double fPulseLength;   // 1PE pulse length in us
      double fPeakTime;      // Time when the pulse reaches its maximum in us
      double fMaxAmplitude;  // Maximum amplitude of the pulse in mV
      double fFrontTime;     // Constant in the exponential function
                             // of the leading edge in us
      double fBackTime;      // Constant in the exponential function
                             // of the tail in us

      // Make sure the FHiCL parameters make sense
      void CheckFHiCLParameters() const;

      std::vector< double > fSinglePEWaveform;
      void CreateSinglePEWaveform();
    
      // Produce waveform on one of the optical detectors
      void CreatePDWaveform(art::Ptr<sim::OpDetBacktrackerRecord> const& btr_p,
                            geo::Geometry const& geometry,
                            std::vector< std::vector< double > >& pdWaveforms,
                            std::vector<FocusList>& fls,
                            sim::OpDetDivRec& DivRec,
                            bool Reflected);

      // Vary the pedestal
      void AddLineNoise(std::vector< std::vector< double > >&,
                        const std::vector<FocusList>& fls) const;

      void AddDarkNoise(std::vector< std::vector< double > >&,
                        std::vector<FocusList>& fls) const;
    
      unsigned short CrossTalk() const;

      // Create a vector of shorts from a vector of doubles
      // rounding it properly
      std::vector< short > VectorOfDoublesToVectorOfShorts(std::vector< double > const&) const;

      // Make several shorter waveforms out of a long one using a hit finder,
      // recording also when they start in the long waveform
      std::map< size_t, std::vector< short > > SplitWaveform(std::vector< short > const&,
                                                             const FocusList&);

      double GetDriftWindow(detinfo::DetectorPropertiesData const& detProp) const;
    
      // Convert time to ticks or the other way around
      // without any checks
      double  TickToTime(size_t tick) const;
      size_t TimeToTick(double  time) const;
    

      bool Detected(int opDet, int &readoutChannel, bool Reflected);
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
    : EDProducer{pset}
    , vInputModules(pset.get< std::vector<std::string> >("InputModules"))
    , fInputModules(vInputModules.begin(), vInputModules.end())
  {


    // Read the fcl-file
//    auto tempvec        = 
//           = 
    fVoltageToADC       = pset.get< double  >("VoltageToADC"      );
    fLineNoiseRMS       = pset.get< double  >("LineNoiseRMS"      );
    fDarkNoiseRate      = pset.get< double  >("DarkNoiseRate"     );
    fCrossTalk          = pset.get< double  >("CrossTalk"         );
    fPedestal           = pset.get< short  >("Pedestal"          );
    fDefaultSimWindow   = pset.get< bool   >("DefaultSimWindow"  );
    fFullWaveformOutput = pset.get< bool   >("FullWaveformOutput");
    fReadoutWindow      = pset.get< size_t >("ReadoutWindow"     );
    fPreTrigger         = pset.get< size_t >("PreTrigger"        );
    
    fPadding            = pset.get< int    >("Padding"           );
    
    fPulseLength        = pset.get< double >("PulseLength"       );
    fPeakTime           = pset.get< double >("PeakTime"          );
    fMaxAmplitude       = pset.get< double >("MaxAmplitude"      );
    fFrontTime          = pset.get< double >("FrontTime"         );
    fBackTime           = pset.get< double >("BackTime"          );
    fDigiTree_SSP_LED   = pset.get< bool   >("SSP_LED_DigiTree"  );
    fUseSDPs            = pset.get< bool   >("UseSDPs", true     );
    fSPEDataFile        = pset.get< std::string >("SPEDataFile");
    TestbenchSinglePE       = pset.get< bool   >("TestbenchSinglePE"     );

    if (!fUseSDPs) {
      throw art::Exception(art::errors::UnimplementedFeature) << "SimPhotonsLite is now deprecated in favor SDPs. If you do not have SDPs because your input file is old, use an older version of dunetpc to run this digitizer";
    }

    // Replace QE with effective area. Get area from OpDet geometry
    double tempQE       = pset.get< double >("QEOverride",    0     );
    double tempRefQE    = pset.get< double >("QERefOverride", 0     );

    fQEOverride         = 0;
    fRefQEOverride      = 0;
    
    if (tempQE > 0) {
      mf::LogInfo("OpDetDigitizerDUNE") << "Overriding QE. All the functions of OpDetResponseService being ignored.";
        
      // Correct out the prescaling applied during simulation
      auto const *LarProp = lar::providerFrom<detinfo::LArPropertiesService>();
      fQEOverride = tempQE / LarProp->ScintPreScale();
      
      if (fQEOverride > 1.0001 ) {
        mf::LogError("OpDetDigitizerDUNE") << "Quantum efficiency set as OpDetDigitizerDUNE.QEOverride, " << tempQE
                                           << " is too large.  It is larger than the prescaling applied during simulation, "
                                           << LarProp->ScintPreScale()
                                           << " (fQEOverride = "
                                           << fQEOverride
                                           << ").  Final QE must be equal to or smaller than the QE applied at simulation time.";
        std::abort();
      }
    }
    if (tempRefQE > 0) {
      mf::LogInfo("OpDetDigitizerDUNE") << "Overriding QE. All the functions of OpDetResponseService being ignored.";
        
      // Correct out the prescaling applied during simulation
      auto const *LarProp = lar::providerFrom<detinfo::LArPropertiesService>();
      fRefQEOverride = tempRefQE / LarProp->ScintPreScale();
      
      if (fRefQEOverride > 1.0001 ) {
        mf::LogError("OpDetDigitizerDUNE") << "Quantum efficiency set as OpDetDigitizerDUNE.QERefOverride, " << tempRefQE
                                           << " is too large.  It is larger than the prescaling applied during simulation, "
                                           << LarProp->ScintPreScale()
                                           << " (fRefQEOverride = "
                                           << fRefQEOverride
                                           << ").  Final QE must be equal to or smaller than the QE applied at simulation time.";
        std::abort();
      }
    }

    fThreshAlg = std::make_unique< pmtana::AlgoSSPLeadingEdge >
      (pset.get< fhicl::ParameterSet >("algo_threshold"));

    if(fDigiTree_SSP_LED){
      art::ServiceHandle< art::TFileService > tfs;
      arvore2 = tfs->make<TTree>("PhotonData", "Photon_analysis");
      arvore2->Branch("photon_opCh",&op_photon);
      arvore2->Branch("photon_pulse",&t_photon);
    }

    // This module produces (infrastructure piece)
    produces< std::vector< raw::OpDetWaveform > >();
    produces<std::vector<sim::OpDetDivRec> > (); 


    // Obtaining parameters from the DetectorClocksService
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
    fSampleFreq = clockData.OpticalClock().Frequency();

    if (fDefaultSimWindow)
    {
      auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob(clockData);

      // Assume the readout starts at -1 drift window
      fTimeBegin = -1*GetDriftWindow(detProp);

      // Take the TPC readout window size and convert
      // to us with the electronics clock frequency
      fTimeEnd   = detProp.ReadOutWindowSize() / clockData.TPCClock().Frequency();
    }
    else
    {
      fTimeBegin = pset.get< double >("TimeBegin");
      fTimeEnd   = pset.get< double >("TimeEnd"  );
    }

    CheckFHiCLParameters();
    
    // Initializing random number engines
    unsigned int seed = pset.get< unsigned int >("Seed", sim::GetRandomNumberSeed());
    auto& engine = createEngine(seed);
    fRandGauss       = std::make_unique< CLHEP::RandGauss       >(engine);
    fRandExponential = std::make_unique< CLHEP::RandExponential >(engine);
    fRandFlat        = std::make_unique< CLHEP::RandFlat        >(engine);
    
    // Creating a single photoelectron waveform
    
      // Hardcoded, probably need to read them from the FHiCL file
      //fPulseLength  = 4.0;
      //fPeakTime     = 0.260;
      //fMaxAmplitude = 0.12;
      //fFrontTime    = 0.009;
      //fBackTime     = 0.476;
      CreateSinglePEWaveform();
    
   
  }

  

  //---------------------------------------------------------------------------
  void OpDetDigitizerDUNE::produce(art::Event& evt)
  {
    
    if(fDigiTree_SSP_LED) art::ServiceHandle< art::TFileService > tfs;

    // A pointer that will store produced OpDetWaveforms
    //    std::unique_ptr< std::vector< raw::OpDetWaveform > >
    //      std::make_unique< std::vector< raw::OpDetWaveform > >());
    auto wave_forms_p = std::make_unique< std::vector< raw::OpDetWaveform > >();
    auto bt_DivRec_p  = std::make_unique< std::vector< sim::OpDetDivRec > >();


    // Total number of ticks in the full waveform
    // Including one pretrigger window before the waveform
    // and one readout window - pretrigger window after the waveform
    unsigned int nSamples = (fTimeEnd - fTimeBegin)*fSampleFreq
      + fReadoutWindow;

    // Geometry service
    art::ServiceHandle< geo::Geometry > geometry;
    
    auto const btr_handles = evt.getMany<std::vector<sim::OpDetBacktrackerRecord>>();

    if (btr_handles.size() == 0)
      throw art::Exception(art::errors::ProductNotFound)<<"No OpDetBacktrackerRecords retrieved.";

    for (auto btr_handle: btr_handles) {
      // Do some checking before we proceed
      if (!btr_handle.isValid()) continue;
      if (!fInputModules.count(btr_handle.provenance()->moduleLabel())) continue;

      bool Reflected = (btr_handle.provenance()->productInstanceName() == "Reflected");

      
      std::vector<art::Ptr<sim::OpDetBacktrackerRecord>> btr_vec;
      art::fill_ptr_vector(btr_vec, btr_handle);
    

      // For every optical detector:
      for (auto const& btr : btr_vec)
      {
        int opDet = btr->OpDetNum();
        //unsigned int opDet = btr->OpDetNum();
      
        // Get number of channels in this optical detector
        unsigned int nChannelsPerOpDet = geometry->NOpHardwareChannels(opDet);
      
        std::vector<FocusList> fls(nChannelsPerOpDet, FocusList(nSamples, fPadding));
        //std::vector<sim::OpDetDivRec> DivRec;
        sim::OpDetDivRec DivRec(opDet);
        //DivRec.chans.resize(nChannelsPerOpDet);
      
        // This vector stores waveforms created for each optical channel
        std::vector< std::vector< double > > pdWaveforms(nChannelsPerOpDet,
                                                         std::vector< double >(nSamples, static_cast< double >(fPedestal)));
      
        CreatePDWaveform(btr, *geometry, pdWaveforms, fls, DivRec, Reflected);
        //DivRec comes out with all of the ticks filled correctly, with each channel filled in it's map.
        //Break here to investigate div recs as they are made and compare them to btrs
      
      
        // Generate dark noise //I will not at this time include dark noise in my split backtracking records.
        if (fDarkNoiseRate > 0.0) AddDarkNoise(pdWaveforms, fls);
      
        // Uncomment to undo the effect of FocusLists. Replaces the accumulated
        // lists with ones asserting we need to look at the whole trace.
        // for(FocusList& fl: fls){
        //        fl.ranges.clear();
        //        fl.ranges.emplace_back(0, nSamples-1);
        // }
      
        // Vary the pedestal
        if (fLineNoiseRMS > 0.0)  AddLineNoise(pdWaveforms, fls);

        // Loop over all the created waveforms, split them into shorter
        // waveforms and use them to initialize OpDetWaveforms
        for (unsigned int hardwareChannel = 0;
             hardwareChannel < nChannelsPerOpDet; ++hardwareChannel)
        {
          for(const std::pair<int, int>& p: fls[hardwareChannel].ranges){
            // It's a shame we copy here. We could actually avoid by making the
            // functions below take a begin()/end() pair.
            const std::vector<double> sub(pdWaveforms[hardwareChannel].begin()+p.first,
                                          pdWaveforms[hardwareChannel].begin()+p.second+1);
          
            std::vector< short > waveformOfShorts = VectorOfDoublesToVectorOfShorts(sub);
          
            std::map< size_t, std::vector < short > > mapTickWaveform =
              (!fFullWaveformOutput) ?
              SplitWaveform(waveformOfShorts, fls[hardwareChannel]) :
              std::map< size_t, std::vector< short > >{ std::make_pair(0,
                                                                       waveformOfShorts) };

            unsigned int opChannel = geometry->OpChannel(opDet, hardwareChannel);

            for (auto const& pairTickWaveform : mapTickWaveform)
            {
              double timeStamp =
                static_cast< double >(TickToTime(pairTickWaveform.first+p.first));

              raw::OpDetWaveform adcVec(timeStamp, opChannel,
                                        pairTickWaveform.second.size());
              
              for (short const& value : pairTickWaveform.second){
                adcVec.emplace_back(value);
              }

              //              wave_forms_p->emplace_back(std::move(adcVec));
              wave_forms_p->push_back(std::move(adcVec));
              //art::Ptr< raw::OpDetWaveform > wave_form_p = make_wave_form_ptr(wave_forms_p->size()-1);
              //bt_wave_assns->addSingle(btr, wave_form_p, DivRec/*DIVREC*/);
            }
          }
        }
        bt_DivRec_p->push_back(std::move(DivRec));
      }
    }
    
    if(fDigiTree_SSP_LED){
      arvore2->Fill();
      t_photon.clear();
      op_photon.clear();
    }


    // Push the OpDetWaveforms into the event
    evt.put(std::move(wave_forms_p));
    if(fUseSDPs){ evt.put(std::move(bt_DivRec_p));}


  }

  //---------------------------------------------------------------------------
  void OpDetDigitizerDUNE::AddPulse(size_t timeBin,
      int scale, std::vector< double >& waveform,
      FocusList& fl) const
  {

    // How many bins will be changed
    size_t pulseLength = fSinglePEWaveform.size();
    if ((timeBin + fSinglePEWaveform.size()) > waveform.size())
      pulseLength = (waveform.size() - timeBin);

    fl.AddRange(timeBin, timeBin+pulseLength-1);

    // Adding a pulse to the waveform
    for (size_t tick = 0; tick != pulseLength; ++tick)
      waveform[timeBin + tick] += scale*fSinglePEWaveform[tick];

  }

  //---------------------------------------------------------------------------
 double OpDetDigitizerDUNE::Pulse1PE(double time) const   
  {
  
    if (time < fPeakTime) return
      (fVoltageToADC*fMaxAmplitude*std::exp((time - fPeakTime)/fFrontTime));
    else return
      (fVoltageToADC*fMaxAmplitude*std::exp(-(time - fPeakTime)/fBackTime));

  }

  //---------------------------------------------------------------------------
  void OpDetDigitizerDUNE::CreateSinglePEWaveform()
  {

    if (TestbenchSinglePE) {
      std::ifstream SPEData;
      SPEData.open(fSPEDataFile);
      if (SPEData.is_open()) {
       mf::LogDebug("OpDetDigitizerDUNE") << " using testbench pe response";
       std::vector< double > SinglePEVec_x;   //1 column
       Double_t  x; 
       while (SPEData >> x ) { SinglePEVec_x.push_back(x); } 
       fSinglePEWaveform = SinglePEVec_x;
             
       std::cout << " out "<<" using TESTbench spe "<< std ::endl;
        
       fPulseLength = fSinglePEWaveform.size();
       SPEData.close(); 
       return;    
      }
       else {
      throw cet::exception("OpDetDigitizerDUNE") << "No Waveform File: Cannot open SPE template file.\n"; 
      }       
   }
    else {
      //shape of single pulse 
       mf::LogDebug("OpDetDigitizerDUNE") << " ideal pe response";
       size_t length = static_cast< size_t > (std::round(fPulseLength*fSampleFreq));
       fSinglePEWaveform.resize(length);
       for (size_t tick = 0; tick != length; ++tick){
        fSinglePEWaveform[tick] =
        Pulse1PE(static_cast< double >(tick)/fSampleFreq);
       }
      std::cout << " out "<<" using ideal spe "<< std ::endl;
   } 
 }

  //---------------------------------------------------------------------------
  void OpDetDigitizerDUNE::CreatePDWaveform
    (art::Ptr<sim::OpDetBacktrackerRecord> const& btr_p,
     geo::Geometry const& geometry,
     std::vector< std::vector< double > >& pdWaveforms,
     std::vector<FocusList>& fls,
     sim::OpDetDivRec& DivRec,
     bool Reflected)
  {

      int const opDet = btr_p->OpDetNum();
      //unsigned int const opDet = btr_p->OpDetNum();
      // This is int because otherwise detectedLite doesn't work
      int readoutChannel;
      // For a group of photons arriving at the same time this is a map
      // of < arrival time (in ns), number of photons >
      //      std::map< int, int > const& photonsMap = litePhotons.DetectedPhotons;
      //sim::timePDclockSDPs_t time_sdps_vector = btr_p->timePDclockSDPsMap();
      //int time_sdps_vector = btr_p->timePDclockSDPsMap();
      auto time_sdps_vector = btr_p->timePDclockSDPsMap();
      /*
         for(auto& divchan : DivRec.chans){
         if(divchan.tick_photons_frac.size()<time_sdps_vector.size())
         divchan.tick_photons_frac.resize(time_sdps_vector.size(), 0.0);
         }*/

      // For every pair of (arrival time, number of photons) in the map:
      //for (auto const& pulse : photonsMap)
      //for (auto const& time_sdps : time_sdps_vector)
      for (size_t i = 0; i<time_sdps_vector.size(); ++i) //I is now the time bin.
      {
        auto time_sdps = time_sdps_vector[i];
        // Converting ns to us
        //        double photonTime = static_cast< double >(time_sdps.first)/1000.0;
        //Really need to do something better here. This conversion is fragile, and makes the populating of the DivRecs confusing.
        double photonTime = time_sdps.first/1000.0;//This should be done with the timing service
        //for (int i = 0; i < pulse.second; ++i)
        //for (size_t i = 0; i < time_sdps.second.size(); ++i)
        for (auto const& sdp : time_sdps.second)
        {
          int tid = sdp.trackID;
          for(int j=0; j<sdp.numPhotons;++j)
          {
            if ((photonTime >= fTimeBegin) && (photonTime < fTimeEnd))
            {
              // Sample a random subset according to QE
              if ( Detected(opDet, readoutChannel, Reflected) )
              {
                unsigned int hardwareChannel =
                  geometry.HardwareChannelFromOpChannel(readoutChannel);
                // Convert the time of the pulse to ticks
                size_t timeBin = TimeToTick(photonTime);
                // Add 1 pulse to the waveform
                AddPulse(timeBin, CrossTalk(), pdWaveforms.at(hardwareChannel), fls[hardwareChannel]);

                unsigned int opChannel = geometry.OpChannel(opDet, hardwareChannel);
                //Set/find tick. Set/find Channel
                sim::OpDet_Time_Chans::stored_time_t tmp_time=time_sdps.first;
                DivRec.AddPhoton(opChannel, tid, tmp_time);
                if(fDigiTree_SSP_LED){
                  op_photon.emplace_back(opChannel);
                  t_photon.emplace_back(photonTime); //vitor: devo usar o time ou o tick?
                }
              }//else{
              /*  unsigned int hardwareChannel =
                  geometry.HardwareChannelFromOpChannel(readoutChannel);
                  unsigned int opChannel = geometry.OpChannel(opDet, hardwareChannel);
                  DivRec.tick_chans[i].DivChans.IfNotInit(opChannel);
                  }*/
            }
            // else
            // remove this as it fills up logfiles for cosmic-ray runs
            //mf::LogInfo("OpDetDigitizerDUNE")
            //  << "Throwing away an out-of-time photon at " << photonTime << '\n';
          }
        }//end of this sdp.
        /*double total=0.0; //This section is no longer needed. I store the un-normalized values, and normalize them at retrieval.
          for(auto it=DivRec.tick_chans[i].DivChans.bit(); it!= DivRec.tick_chans[i].DivChans.lit(); ++it)
          { //With a bit of work, we could avoid this loop by remembering this number from above.
          total+=it->second.tick_photons_frac;
          }
          for(auto chan : DivRec.tick_chans[i].DivChans.chans())
          {
          DivRec.tick_chans[i].DivChans.NormPhotonFrac(chan, total);
          }
          total=0.0; // This isn't actually needed.*/
      }//End this time tick.

  }

  //---------------------------------------------------------------------------
  void OpDetDigitizerDUNE::
    AddLineNoise(std::vector< std::vector< double > >& waveforms,
        const std::vector<FocusList>& fls) const
    {
      int i = 0;
      for(auto& waveform : waveforms){
        for(unsigned int j = 0; j < fls[i].ranges.size(); ++j){
          const std::pair<int, int>& p = fls[i].ranges[j];
          for(int k = p.first; k <= p.second; ++k){
            waveform[k] += fRandGauss->fire(0, fLineNoiseRMS);
          }
        }

        ++i;
      }
    }

  //---------------------------------------------------------------------------
  void OpDetDigitizerDUNE::
    AddDarkNoise(std::vector< std::vector< double > >& waveforms,
        std::vector<FocusList>& fls) const
    {
      int i = 0;
      for (auto& waveform : waveforms)
      {
        // Multiply by 10^6 since fDarkNoiseRate is in Hz
        double darkNoiseTime = static_cast< double >(fRandExponential->
            fire(1.0/fDarkNoiseRate)*1000000.0) + fTimeBegin;
        while (darkNoiseTime < fTimeEnd)
        {
          size_t timeBin = TimeToTick(darkNoiseTime);
          AddPulse(timeBin, CrossTalk(), waveform, fls[i]);
          // Find next time to simulate a single PE pulse
          darkNoiseTime += static_cast< double >
            (fRandExponential->fire(1.0/fDarkNoiseRate)*1000000.0);
        }

        ++i;
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
  std::vector< short > OpDetDigitizerDUNE::VectorOfDoublesToVectorOfShorts
    (std::vector< double > const& vectorOfDoubles) const
    {
      // Don't bother to round properly, it's faster this way
      return std::vector<short>(vectorOfDoubles.begin(), vectorOfDoubles.end());

      /*
         std::vector< short > vectorOfShorts;
         vectorOfShorts.reserve(vectorOfDoubles.size());

         for (short const& value : vectorOfDoubles)
         vectorOfShorts.emplace_back(static_cast< short >(std::round(value)));

         return vectorOfShorts;
         */
    }

  //---------------------------------------------------------------------------
  std::map< size_t, std::vector< short > > OpDetDigitizerDUNE::
    SplitWaveform(std::vector< short > const& waveform,
        const FocusList& fls)
    {

      std::map< size_t, std::vector< short > > mapTickWaveform;

      ::pmtana::PedestalMean_t  ped_mean (waveform.size(),0);
      ::pmtana::PedestalSigma_t ped_sigma(waveform.size(),0);

      fThreshAlg->Reconstruct(waveform,ped_mean,ped_sigma);

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
  double OpDetDigitizerDUNE::GetDriftWindow(detinfo::DetectorPropertiesData const& detProp) const
  {

    double driftWindow;

    double maxDrift = 0.0;
    for (geo::TPCGeo const& tpc :
           art::ServiceHandle< geo::Geometry >()->Iterate<geo::TPCGeo>())
      if (maxDrift < tpc.DriftDistance()) maxDrift = tpc.DriftDistance();

    driftWindow = maxDrift/detProp.DriftVelocity();

    return driftWindow;

  }

  //---------------------------------------------------------------------------
  double OpDetDigitizerDUNE::TickToTime(size_t tick) const
  {

    if (tick > fPreTrigger)
      return (static_cast< double >(tick - fPreTrigger)/fSampleFreq
          + fTimeBegin);
    else
      return (static_cast< double >(fPreTrigger - tick)/fSampleFreq*(-1.0)
          + fTimeBegin);

  }

  //---------------------------------------------------------------------------
  size_t OpDetDigitizerDUNE::TimeToTick(double time) const
  {

    return static_cast< size_t >(std::round((time - fTimeBegin)*fSampleFreq
          + fPreTrigger));

  }

  //---------------------------------------------------------------------------
  void OpDetDigitizerDUNE::CheckFHiCLParameters() const
  {

    // Are all these logic errors?

    if (fLineNoiseRMS < 0.0)
      throw art::Exception(art::errors::LogicError)
        << "fLineNoiseRMS: " << fLineNoiseRMS << '\n'
        << "Line noise RMS should be non-negative!\n";

    if (fDarkNoiseRate < 0.0)
      throw art::Exception(art::errors::LogicError)
        << "fDarkNoiseRate: " << fDarkNoiseRate << '\n'
        << "Dark noise rate should be non-negative!\n";

    if (fPreTrigger >= fReadoutWindow)
      throw art::Exception(art::errors::LogicError)
        << "PreTrigger: "    << fPreTrigger    << " and "
        << "ReadoutWindow: " << fReadoutWindow << '\n'
        << "Pretrigger window has to be shorter than readout window!\n";

    if (fTimeBegin >= fTimeEnd)
      throw art::Exception(art::errors::LogicError)
        << "TimeBegin: " << fTimeBegin << " and "
        << "TimeEnd: "   << fTimeEnd   << '\n'
        << "TimeBegin should be less than TimeEnd!\n";

  }

  //---------------------------------------------------------------------------
  bool OpDetDigitizerDUNE::Detected(int OpDet, int &readoutChannel, bool Reflected) {

    if (Reflected && !fRefQEOverride) {
      cet::exception("OpDetDigitizerDUNE") << "Cannot handle reflected light if QERefOverride isn't set";
    }
    if ( (!Reflected && fQEOverride > 0) || Reflected ) {
      // Find the Optical Detector using the geometry service
      art::ServiceHandle<geo::Geometry> geom;
      
      // Here OpDet must be opdet since we are introducing
      // channel mapping here.
      float NOpHardwareChannels = geom->NOpHardwareChannels(OpDet);
      int hardwareChannel = (int) ( CLHEP::RandFlat::shoot(1.0) * NOpHardwareChannels );
      readoutChannel = geom->OpChannel(OpDet, hardwareChannel);

      // Check QE
      if (!Reflected) {
        if ( CLHEP::RandFlat::shoot(1.0) > fQEOverride ) return false;
      }
      else {
        if ( CLHEP::RandFlat::shoot(1.0) > fRefQEOverride ) return false;        
      }
      return true;
    }
    else {
      art::ServiceHandle< opdet::OpDetResponseInterface > odResponse;
      // Service for determining optical detector responses
      return odResponse->detectedLite(OpDet, readoutChannel);
    }
  }
  
  
}
