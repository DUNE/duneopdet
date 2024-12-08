//=========================================================
// OpDetDigitizerProtoDUNEHD_module.cc
// This module produces takes many OpDetBacktracker records
// and combined them to create OpDetWaveforms.
// It allows to select selftrigger or fullstreaming channels.
//
// J.Soto
// Based on OpMCDigi_module.cc
//=========================================================

#ifndef OpDetDigitizerProtoDUNEHD_h
#define OpDetDigitizerProtoDUNEHD_h

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
#include "larcore/Geometry/WireReadout.h"
#include "larcore/Geometry/Geometry.h"
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
#include <iostream>
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
    void Reset()
    {
      ranges=std::vector<std::pair<int,int>>({{0,fNSamples-1}});
    }
    
    std::vector<std::pair<int, int>> ranges;
    void Print()
    {
      
      std::cout << "ranges.size():" << ranges.size() << std::endl;
      for(size_t i=0;i<ranges.size();i++) std::cout << "\t" << i << "\t" << ranges[i].first <<"\t" << ranges[i].second << std::endl;
      }

    protected:
      int fNSamples;
      int fPadding;
  };

  class OpDetDigitizerProtoDUNEHD : public art::EDProducer{

  public:

      OpDetDigitizerProtoDUNEHD(fhicl::ParameterSet const&);
      // Should the destructor be empty?
      //      virtual ~OpDetDigitizerProtoDUNEHD();
    
      void produce(art::Event&);

    private:

      // The parameters read from the FHiCL file
      std::vector<std::string> vInputModules;
      std::set<std::string> fInputModules;   // Input tag for OpDet collection
      double  fSampleFreq;                   // Sampling frequency in MHz
      double  fTimeBegin;                    // Beginning of waveform in us, not used if fDefaultSimWindow=true
      double  fTimeEnd;                      // End of waveform in us, not used if fDefaultSimWindow=true
      double  fLineNoiseRMS;                 // Pedestal RMS in ADC counts
      double  fDarkNoiseRate;                // In Hz
      double  fCrossTalk;                    // Probability of SiPM producing 2 PE signal
                                             // in response to 1 photon
      short  fPedestal;                      // In ADC counts
      bool   fDefaultSimWindow;              // Set the start time to -1 drift window and
                                             // the end time to the end time
                                             // of the TPC readout
      size_t fReadoutWindow;                 // In ticks
      size_t fPreTrigger;                    // In ticks

      int    fPadding;                       // In ticks, used by the focus list
    
      bool   fDigiTree_SSP_LED;              // To create a analysis Tree for SSP LED
      bool   fUseSDPs;                       // = pset.get< bool   >("UseSDPs", true);

      bool   fCustomPDE;
      bool   fCustomSPEAmplitude;
      double fPDE;
      std::vector<double> fPDEVector;
      double fSPEAmplitude;
      std::vector<double> fSPEAmplitudeVector;

      double  getPDE(int OpDet) const;
      double  getSPEAmplitude(int OpDet) const;
      double fSPEAmplitudesError;
      double fSPEWidth;
      short fDynamicRangeSaturation;
      bool fNegativeSignal;
      
      std::vector<int> fFullStreamingChannels;
      std::vector<std::pair<std::string,std::vector<int>>>fSPETemplateMap;

      std::vector<std::vector<double>> fSinglePEWaveforms;
      std::map<int,int> ChannelToSPE;
      
      //-----------------------------------------------------
      // Trigger analysis variables
      std::vector<double> t_photon; // vitor
      std::vector<int>    op_photon;

      TTree *arvore2;
      //-----------------------------------------------------


      // Random number engines
      std::unique_ptr< CLHEP::RandGauss       > fRandGauss;
      std::unique_ptr< CLHEP::RandExponential > fRandExponential;
      std::unique_ptr< CLHEP::RandFlat        > fRandFlat;

      // Function that adds n pulses to a waveform
      void AddPulse(size_t timeBin, int scale,
                    std::vector< double >& waveform,
                    FocusList& fl,int opDet) const;

      // Make sure the FHiCL parameters make sense
      void CheckFHiCLParameters() const;

//      std::vector< double > fSinglePEWaveform;
      void CreateSinglePEWaveforms();
    
      // Produce waveform on one of the optical detectors
      void CreatePDWaveform(art::Ptr<sim::OpDetBacktrackerRecord> const& btr_p,
                            geo::WireReadoutGeom const& wireReadout,
                            std::vector< std::vector< double > >& pdWaveforms,
                            std::vector<FocusList>& fls,
                            sim::OpDetDivRec& DivRec);

      // Vary the pedestal
      void AddLineNoise(std::vector< std::vector< double > >&,
                        const std::vector<FocusList>& fls) const;

      void AddDarkNoise(std::vector< std::vector< double > >&,
                        std::vector<FocusList>& fls, int opDet) const;
    
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

      void DynamicRangeSaturation(std::vector< short > &wvf);
    

      
      // Waveform tree variables
      bool fExportWaveformTree;
      TTree *fWaveformTree;
      TTree *fRunInfo;

      int Run;
      int SubRun;
      int Event;
      double TimeStamp;
      int OpChannelTree;
      unsigned int nSamplesTree;
      double SampleSize;
      std::vector<short> adc_value;
      unsigned int nOpDet;
      int TriggerType;
      
      int fSelfTrigger_Pretrigger;
      int fSelfTrigger_ReadoutWindow;
      int fSelfTrigger_DaphneThreshold;
      std::vector<std::pair<double,std::vector<int>>>fSelfTrigger_DaphneThresholdMapV;
      std::map<int,double>fSelfTrigger_DaphneThresholdMap;

      std::map< size_t, std::vector< short > > DAPHNESelfTriggerSplit(std::vector< short > const& waveform,
        const FocusList& fls, int opdet);

      double getDaphneThreshold(int OpDet) const ;
  };
}

#endif

namespace opdet {
  
  DEFINE_ART_MODULE(OpDetDigitizerProtoDUNEHD)

}

namespace opdet {

  //---------------------------------------------------------------------------	
  // Constructor
  OpDetDigitizerProtoDUNEHD::OpDetDigitizerProtoDUNEHD(fhicl::ParameterSet const& pset)
    : EDProducer{pset}
    , vInputModules(pset.get< std::vector<std::string> >("InputModules"))
    , fInputModules(vInputModules.begin(), vInputModules.end())
  {

    // Read the fcl-file
//    auto tempvec        = 
//           = 
    fLineNoiseRMS       = pset.get< double  >("LineNoiseRMS"      );
    fDarkNoiseRate      = pset.get< double  >("DarkCountRate"     );
    fCrossTalk          = pset.get< double  >("CrossTalk"         );
    fPedestal           = pset.get< short  >("Pedestal"          );
    fDefaultSimWindow   = pset.get< bool   >("DefaultSimWindow"  ); //long readout window (2 drift windows).
    fPreTrigger         = pset.get< size_t >("PreTrigger"        ,0); //set to zero if TimeBegin<
    fPadding            = pset.get< int    >("Padding"           );
    fDynamicRangeSaturation = pset.get< short    >("DynamicRangeSaturation");
    fNegativeSignal = pset.get< bool    >("NegativeSignal");
    
    fDigiTree_SSP_LED   = pset.get< bool   >("SSP_LED_DigiTree"  );
    fUseSDPs            = pset.get< bool   >("UseSDPs", true     );
    fExportWaveformTree = pset.get<bool>("ExportWaveformTree",false);

    fCustomPDE           = pset.get< bool   >("CustomPDEperOpDet",  false    ); //apply a custom PDE per channel
    fCustomSPEAmplitude         = pset.get< bool   >("CustomSPEAmplitudeperOpDet",false    ); //apply a custom SPEAmplitude per channel


    if(fCustomPDE) fPDEVector = pset.get<std::vector<double>>("PDEVector"); //only used if CustomPDEperOpDet is true
    else fPDE = pset.get< double >("PDE",3);
    if(fCustomSPEAmplitude) fSPEAmplitudeVector = pset.get<std::vector<double>>("SPEAmplitudeVector"); //only used if CustomSPEAmplitudeperOpDet is true
    else fSPEAmplitude          = pset.get< double >("SPEAmplitude",15);///SPE amplitude in ADC per PE Not used at the moment
    
    fFullStreamingChannels =pset.get< std::vector<int>>("FullStreamingChannels"); //channels that are not here, are self-trigger.
    fSelfTrigger_Pretrigger = pset.get<int>("SelfTrigger_Pretrigger",200); //Pretrigger samples in Selftrigger channels.
    fSelfTrigger_ReadoutWindow = pset.get<int>("SelfTrigger_ReadoutWindow",1000); //channels that are not here, are self-trigger.
    fSelfTrigger_DaphneThreshold = pset.get<int>("SelfTrigger_DaphneThreshold",65); //Daphne selftrigger threshold default value.
    fSelfTrigger_DaphneThresholdMapV = pset.get<std::vector<std::pair<double,std::vector<int>>>>("SelfTrigger_DaphneThresholdMap"); //Daphne self-trigger threshold custom value.
    for (auto p : fSelfTrigger_DaphneThresholdMapV)
    {
      for (auto ch : p.second){fSelfTrigger_DaphneThresholdMap[ch]=p.first;}
    }

    fSPETemplateMap     = pset.get< std::vector<std::pair<std::string,std::vector<int>>>>("SPETemplateMap");
    //std::map SPE file -> Vector of channels using this map

    

    if (!fUseSDPs) {
      throw art::Exception(art::errors::UnimplementedFeature) << "SimPhotonsLite is now deprecated in favor SDPs. If you do not have SDPs because your input file is old, use an older version of dunetpc to run this digitizer";
    }


    art::ServiceHandle< art::TFileService > tfs;
    if(fDigiTree_SSP_LED){
      arvore2 = tfs->make<TTree>("PhotonData", "Photon_analysis");
      arvore2->Branch("photon_opCh",&op_photon);
      arvore2->Branch("photon_pulse",&t_photon);
    }
    

    // This module produces (infrastructure piece)
    produces< std::vector< raw::OpDetWaveform > >();
    produces<std::vector<sim::OpDetDivRec> > (); 

    CreateSinglePEWaveforms();
    // Obtaining parameters from the DetectorClocksService
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();

    fSampleFreq = clockData.OpticalClock().Frequency(); //in MHz
    SampleSize = 1000.0/fSampleFreq; //in ns

    if (fDefaultSimWindow)
    {
      auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob(clockData);

      // Assume the readout starts at -1 drift window
      fTimeBegin = -1*GetDriftWindow(detProp);

      // Take the TPC readout window size and convert
      // to us with the electronics clock frequency
      fTimeEnd   = detProp.ReadOutWindowSize() / clockData.TPCClock().Frequency();
      fReadoutWindow = (fTimeEnd- fTimeBegin)*fSampleFreq;
      fPreTrigger =0.0;
      
    }
    else
    {
      fTimeBegin = pset.get< double >("TimeBegin"); //in us
      fTimeEnd   = pset.get< double >("TimeEnd"  ); //in us
      fReadoutWindow = (fTimeEnd- fTimeBegin)*clockData.OpticalClock().Frequency(); //ticks
      fPreTrigger =0.0;

    }
    CheckFHiCLParameters();
    
    // Initializing random number engines
    unsigned int seed = pset.get< unsigned int >("Seed", sim::GetRandomNumberSeed());
    auto& engine = createEngine(seed);
    fRandGauss       = std::make_unique< CLHEP::RandGauss       >(engine);
    fRandExponential = std::make_unique< CLHEP::RandExponential >(engine);
    fRandFlat        = std::make_unique< CLHEP::RandFlat        >(engine);
    
//    art::ServiceHandle< geo::Geometry > geometry;
//    nOpDet=geometry->NOpDets();
    nOpDet = art::ServiceHandle<geo::WireReadout>()->Get().NOpChannels();
    if(fExportWaveformTree)
    {
      double *vaux; vaux = (double*)malloc(sizeof(double)*nOpDet);
      double *vaux2; vaux2 = (double*)malloc(sizeof(double)*nOpDet);
      fRunInfo = tfs->make<TTree>("RunInfo","MonteCarlo Run Info");
      fRunInfo->Branch("SampleSize"   , &SampleSize   , "SampleSize/D"     );
      fRunInfo->Branch("nOpDet"        , &nOpDet   , Form("nOpDet/I")    );
      fRunInfo->Branch("PDE"        , vaux   , Form("PDE[nOpDet]/D")    );
      fRunInfo->Branch("SPEAmplitude"        , vaux2   , Form("SPEAmplitude[nOpDet]/D")    );
      for (unsigned int i=0; i<nOpDet; i++) vaux[i]=getPDE(i);
      for (unsigned int i=0; i<nOpDet; i++) vaux2[i]=getSPEAmplitude(i);
      for (unsigned int i=0; i<nOpDet; i++) std::cout << " PDE [ " << i << " ] = " << vaux[i] << " percent." << std::endl;
      for (unsigned int i=0; i<nOpDet; i++) std::cout << " SPEAmplitude [ " << i << " ] = " << vaux2[i] << " ADC" << std::endl;
      for (unsigned int i=0; i<nOpDet; i++) std::cout << " DaphneThreshold [ " << i << " ] = " << getDaphneThreshold(i) << " ADC" << std::endl;

      fRunInfo->Fill();

      fWaveformTree = tfs->make<TTree>("WaveformTree","Waveforms Tree");
      fWaveformTree->Branch("Run"       , &Run       , "Run/I"       );
      fWaveformTree->Branch("SubRun"    , &SubRun    , "SubRun/I"    );
      fWaveformTree->Branch("Event"     , &Event     , "Event/I"     );
      fWaveformTree->Branch("Trigger"     , &TriggerType     , "Trigger/I"     );
      fWaveformTree->Branch("TimeStamp" , &TimeStamp     , "TimeStamp/D"     );
      fWaveformTree->Branch("NSamples"     , &nSamplesTree     , "NSamples/I"     );
      fWaveformTree->Branch("OpChannel"     , &OpChannelTree     , "OpChannel/I"     );
      fWaveformTree->Branch("adc", &adc_value);
    }
    std::cout << "Generating waveforms of " << fTimeEnd-fTimeBegin << "us = "<< fReadoutWindow <<" Samples"<< std::endl;
    std::cout << "\tTimeBegin: " << fTimeBegin <<" "<< std::endl;
    std::cout << "\tfTimeEnd: " << fTimeEnd <<" "<< std::endl;
    std::cout << "\tSampleFreq: " << fSampleFreq <<" MHz"<< std::endl;
    std::cout << "\tReadoutWindow: " << fReadoutWindow <<" ticks"<< std::endl;
    std::cout << "\tPreTrigger: " << fPreTrigger <<" ticks"<< std::endl;
    std::cout << "\tSampleSize: " << SampleSize <<" ns"<< std::endl;

//    std::cout << "\fSinglePEWaveform: " << fSinglePEWaveform.size() <<std::endl;

    std::cout << "Full streaming channels: "; for (auto i : fFullStreamingChannels) std::cout << i <<" "; std::cout << std::endl; 

    std::cout << "SPE Templates: " << std::endl;
    
    for (auto p : fSPETemplateMap)
    {
      std::cout << p.first <<" ";
      for (auto ch : p.second){std::cout << ch << " ";}
      std::cout << std::endl; 
    }
  }
  

  //---------------------------------------------------------------------------
  void OpDetDigitizerProtoDUNEHD::produce(art::Event& evt)
  {
    // Geometry service
    auto const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
    
    if(fExportWaveformTree)
    {
      Run    = evt.run();
      SubRun = evt.subRun();
      Event  = evt.event();
    }

    auto wave_forms_p = std::make_unique< std::vector< raw::OpDetWaveform > >();
    auto bt_DivRec_p  = std::make_unique< std::vector< sim::OpDetDivRec > >();

    //We can have several handles of OpDetBacktrackerRecords, thus first we do a map ODBMap to locate them
    //on a per channel logic to facilitate the waveform creation.

    std::vector<std::vector<const art::Ptr<sim::OpDetBacktrackerRecord>*>> ODBMap(nOpDet,std::vector<const art::Ptr<sim::OpDetBacktrackerRecord>*>());
    auto const btr_handles = evt.getMany<std::vector<sim::OpDetBacktrackerRecord>>();
    if (btr_handles.size() == 0)
      throw art::Exception(art::errors::ProductNotFound)<<"No OpDetBacktrackerRecords retrieved.";
    std::vector<std::vector<art::Ptr<sim::OpDetBacktrackerRecord>>> mybtr_vec;

    for (size_t hh=0; hh<btr_handles.size(); hh++)
    {
      if (!btr_handles[hh].isValid()) continue;
      if (!fInputModules.count(btr_handles[hh].provenance()->moduleLabel())) continue; //avoid modules not added to the list
      mybtr_vec.push_back(std::vector<art::Ptr<sim::OpDetBacktrackerRecord>>());
      art::fill_ptr_vector(mybtr_vec[mybtr_vec.size()-1], btr_handles[hh]);
      for (auto const& btr : mybtr_vec[mybtr_vec.size()-1])
      {
        ODBMap[btr.get()->OpDetNum()].push_back(&btr);
      }
    }
    
    //now we start creating waveforms    
    for (size_t opDet=0; opDet<ODBMap.size();opDet++)
    {
      unsigned int nChannelsPerOpDet = wireReadout.NOpHardwareChannels(opDet);      
      std::vector<FocusList> fls(nChannelsPerOpDet, FocusList(fReadoutWindow, fPadding));
      sim::OpDetDivRec DivRec(opDet);

      // This vector stores waveforms created for all hardware channels in OpDet.
      std::vector< std::vector< double > > pdWaveforms(nChannelsPerOpDet,
                                                 std::vector< double >(fReadoutWindow, static_cast< double >(fPedestal)));

      //We add all photons from ODBs to the waveform.
      for (size_t jj=0; jj<ODBMap[opDet].size();jj++)//(auto btr : btr_vec)
      {        
        auto btr = ODBMap[opDet][jj];
        if(btr->get()->OpDetNum()!=(int)opDet)
          throw art::Exception(art::errors::LogicError)
            << "Memory issues in the ODBMap!\n";
        CreatePDWaveform(*btr, wireReadout, pdWaveforms, fls, DivRec);
      }

      bool SelfTrigger=false;
      if(std::find(fFullStreamingChannels.begin(), fFullStreamingChannels.end(), opDet)
          ==fFullStreamingChannels.end()) SelfTrigger=true;

      if(!SelfTrigger) for (unsigned int hardwareChannel = 0;
         hardwareChannel < nChannelsPerOpDet; ++hardwareChannel) fls[hardwareChannel].Reset();
      
      // Generate dark noise
      if (fDarkNoiseRate > 0.0) AddDarkNoise(pdWaveforms, fls,opDet);
            
      // Vary the pedestal
      if (fLineNoiseRMS > 0.0)  AddLineNoise(pdWaveforms, fls);

      // Loop over all the created waveforms, split them into shorter
      // waveforms and use them to initialize OpDetWaveforms
      for (unsigned int hardwareChannel = 0;
           hardwareChannel < nChannelsPerOpDet; ++hardwareChannel)
      {
        for(const std::pair<int, int>& p: fls[hardwareChannel].ranges){
          const std::vector<double> sub(pdWaveforms[hardwareChannel].begin()+p.first,
                                        pdWaveforms[hardwareChannel].begin()+p.second+1);
          
          std::vector< short > waveformOfShorts = VectorOfDoublesToVectorOfShorts(sub);
          DynamicRangeSaturation(waveformOfShorts);
          
          std::map< size_t, std::vector < short > > mapTickWaveform =
            (SelfTrigger) ?
            DAPHNESelfTriggerSplit(waveformOfShorts, fls[hardwareChannel],opDet) :
            std::map< size_t, std::vector< short > >{ std::make_pair(0,
                                                                     waveformOfShorts) };
                                                                
          unsigned int opChannel = wireReadout.OpChannel(opDet, hardwareChannel);
          for (auto const& pairTickWaveform : mapTickWaveform)
          {
            double timeStamp =
              static_cast< double >(TickToTime(pairTickWaveform.first+p.first));
            raw::OpDetWaveform adcVec(timeStamp, opChannel,
                                      pairTickWaveform.second.size());
              
            for (short const& value : pairTickWaveform.second){
                adcVec.emplace_back(value);
            }
            wave_forms_p->push_back(std::move(adcVec));
          }
        }
        if(fUseSDPs) bt_DivRec_p->push_back(std::move(DivRec));
      }
    }
    
    if(fDigiTree_SSP_LED){
      arvore2->Fill();
      t_photon.clear();
      op_photon.clear();
    }
    if(fExportWaveformTree)
    {	
      for(unsigned int jj=0; jj < wave_forms_p->size(); jj++)
      {
        TriggerType=0;
        OpChannelTree = wave_forms_p->at(jj).ChannelNumber();
        nSamplesTree=wave_forms_p->at(jj).Waveform().size();
        TimeStamp=wave_forms_p->at(jj).TimeStamp();
        adc_value.resize(nSamplesTree);
        for (unsigned int ii = 0; ii< wave_forms_p->at(jj).Waveform().size() ; ii++)
        {                    
          adc_value[ii] = wave_forms_p->at(jj).Waveform()[ii];
        }
        fWaveformTree->Fill();
      }
    }


    // Push the OpDetWaveforms into the event
    evt.put(std::move(wave_forms_p));
    if(fUseSDPs){ evt.put(std::move(bt_DivRec_p));}
  }

  //---------------------------------------------------------------------------
  void OpDetDigitizerProtoDUNEHD::AddPulse(size_t timeBin,
      int scale, std::vector< double >& waveform,
      FocusList& fl, int opCh) const
  {

    // How many bins will be changed
    auto & thisSPE=fSinglePEWaveforms[ChannelToSPE.at(opCh)];
    size_t pulseLength = thisSPE.size();
    if ((timeBin + thisSPE.size()) > waveform.size())
      pulseLength = (waveform.size() - timeBin);

    fl.AddRange(timeBin, timeBin+pulseLength-1);

    // Adding a pulse to the waveform
    for (size_t tick = 0; tick != pulseLength; ++tick)
      waveform[timeBin + tick] += scale*thisSPE[tick]*getSPEAmplitude(opCh);

  }

  //---------------------------------------------------------------------------
  void OpDetDigitizerProtoDUNEHD::CreateSinglePEWaveforms()
  {
      double Sign=1.0; //positive by default
      if(fNegativeSignal) Sign=-1.0;
      for(auto v : fSPETemplateMap)
      {
        std::cout << "Using custom SPE response taken from " << v.first << std::endl;
        std::string datafile;
        cet::search_path sp("FW_SEARCH_PATH");
        // taking the file name as the first argument,
        // the second argument is the local variable where to store the full path - both are std::string objects
        sp.find_file(v.first, datafile);
        std::ifstream SPEData;
        SPEData.open(datafile);
        if (SPEData.is_open()) {
          mf::LogDebug("OpDetDigitizerProtoDUNEHD") << " using testbench pe response";
          std::vector< double > SinglePEVec_x;   //1 column
          Double_t  x; Double_t xmax=0;
          while (SPEData >> x ) {
            if(x>xmax)xmax=x;
            SinglePEVec_x.push_back(Sign*x);
          }
          for(size_t i=0; i<SinglePEVec_x.size(); i++) SinglePEVec_x[i]/=xmax; 
          fSinglePEWaveforms.push_back(SinglePEVec_x);
          SPEData.close();
          for( auto ch : v.second) ChannelToSPE[ch]=fSinglePEWaveforms.size()-1;
        }
        else {
          throw cet::exception("OpDetDigitizerProtoDUNEHD") << "No Waveform File: Cannot open SPE template file.\n"; 
        }
    }

 }

   void OpDetDigitizerProtoDUNEHD::DynamicRangeSaturation(std::vector< short > &wvf)
  {
     for(size_t i=0; i<wvf.size();i++)
     {
       if(wvf[i]>fDynamicRangeSaturation) wvf[i]=fDynamicRangeSaturation;
       if(wvf[i]<0) wvf[i]=0;
     }
  }
  //---------------------------------------------------------------------------
  void OpDetDigitizerProtoDUNEHD::CreatePDWaveform
    (art::Ptr<sim::OpDetBacktrackerRecord> const& btr_p,
     geo::WireReadoutGeom const& wireReadout,
     std::vector< std::vector< double > >& pdWaveforms,
     std::vector<FocusList>& fls,
     sim::OpDetDivRec& DivRec)
  {

      int const opDet = btr_p->OpDetNum();
      //unsigned int const opDet = btr_p->OpDetNum();
      // This is int because otherwise detectedLite doesn't work
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
            if ((photonTime+1.0/fSampleFreq >= fTimeBegin) && (photonTime < fTimeEnd-1.0/fSampleFreq))            
            {
              // Sample a random subset according to PDE
              if (CLHEP::RandFlat::shoot(1.0) <getPDE(opDet))
              {
                 float NOpHardwareChannels = wireReadout.NOpHardwareChannels(opDet);
                 int hardwareChannel = (int) ( CLHEP::RandFlat::shoot(1.0) * NOpHardwareChannels );
//                 int readoutChannel = geometry.OpChannel(opDet, hardwareChannel);

                // Convert the time of the pulse to ticks
                size_t timeBin = TimeToTick(photonTime);
                // Add 1 pulse to the waveform
                if(timeBin>=fReadoutWindow) continue;
                AddPulse(timeBin, CrossTalk(), pdWaveforms.at(hardwareChannel), fls[hardwareChannel],hardwareChannel);

                unsigned int opChannel = wireReadout.OpChannel(opDet, hardwareChannel);
                //Set/find tick. Set/find Channel
                sim::OpDet_Time_Chans::stored_time_t tmp_time=time_sdps.first;
                DivRec.AddPhoton(opChannel, tid, tmp_time);
                if(fDigiTree_SSP_LED){
                  op_photon.emplace_back(opChannel);
                  t_photon.emplace_back(photonTime); //vitor: devo usar o time ou o tick?
                }
              }
            }
          }
        }//end of this sdp.
      }//End this time tick.
  }

  //---------------------------------------------------------------------------
  void OpDetDigitizerProtoDUNEHD::
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
  void OpDetDigitizerProtoDUNEHD::
    AddDarkNoise(std::vector< std::vector< double > >& waveforms,
        std::vector<FocusList>& fls, int opDet) const
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
          if(timeBin>=fReadoutWindow) continue;
          AddPulse(timeBin, CrossTalk(), waveform, fls[i],opDet);
          // Find next time to simulate a single PE pulse
          darkNoiseTime += static_cast< double >
            (fRandExponential->fire(1.0/fDarkNoiseRate)*1000000.0);
        }

        ++i;
      }
    }

  //---------------------------------------------------------------------------
  unsigned short OpDetDigitizerProtoDUNEHD::CrossTalk() const
  {
    // Sometimes this should produce 3 or more PEs (not implemented)
    if      (fCrossTalk <= 0.0)                 return 1;
    else if (fRandFlat->fire(1.0) > fCrossTalk) return 1;
    else                                        return 2;
  }

  //---------------------------------------------------------------------------
  std::vector< short > OpDetDigitizerProtoDUNEHD::VectorOfDoublesToVectorOfShorts
    (std::vector< double > const& vectorOfDoubles) const
    {
      return std::vector<short>(vectorOfDoubles.begin(), vectorOfDoubles.end());
    }


  //---------------------------------------------------------------------------
  std::map< size_t, std::vector< short > > OpDetDigitizerProtoDUNEHD::
    DAPHNESelfTriggerSplit(std::vector< short > const& waveform,
        const FocusList& fls, int opDet)
    {
      std::map< size_t, std::vector< short > > mapTickWaveform;
      for(size_t i=0;i< waveform.size();i++)
      {
        if(waveform[i]>getDaphneThreshold(opDet)+fPedestal)
        {

          size_t t_start;
          if(i>(size_t)fSelfTrigger_Pretrigger) t_start=i-fSelfTrigger_Pretrigger;
          else t_start = 0;
          size_t t_end = t_start+fSelfTrigger_ReadoutWindow;
          if(t_end>=waveform.size()) t_end=waveform.size()-1;
          std::vector< short >::const_iterator window_start =
            waveform.begin() + (t_start);
          std::vector< short >::const_iterator window_end   =
            waveform.begin() + (t_end  );
            mapTickWaveform.emplace((t_start),
          std::vector< short >(window_start, window_end));
          mapTickWaveform.emplace((t_start),
             std::vector< short >(window_start, window_end));
          i+=fSelfTrigger_ReadoutWindow-1;
        }
      }
      return mapTickWaveform;
    }


  double OpDetDigitizerProtoDUNEHD::GetDriftWindow(detinfo::DetectorPropertiesData const& detProp) const
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
  double OpDetDigitizerProtoDUNEHD::TickToTime(size_t tick) const
  {

    if (tick > fPreTrigger)
      return (static_cast< double >(tick - fPreTrigger)/fSampleFreq
          + fTimeBegin);
    else
      return (static_cast< double >(fPreTrigger - tick)/fSampleFreq*(-1.0)
          + fTimeBegin);

  }

  //---------------------------------------------------------------------------
  size_t OpDetDigitizerProtoDUNEHD::TimeToTick(double time) const
  {

    return static_cast< size_t >(std::round((time - fTimeBegin)*fSampleFreq
          + fPreTrigger));

  }

  //---------------------------------------------------------------------------
  void OpDetDigitizerProtoDUNEHD::CheckFHiCLParameters() const
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
  double  OpDetDigitizerProtoDUNEHD::getPDE(int OpDet) const
  {
    if(!fCustomPDE) return fPDE;
    else return fPDEVector[OpDet];
  }

  double  OpDetDigitizerProtoDUNEHD::getSPEAmplitude(int OpDet) const
  {
    if(!fCustomSPEAmplitude) return fSPEAmplitude;
    else return fSPEAmplitudeVector[OpDet];
  }

  double OpDetDigitizerProtoDUNEHD::getDaphneThreshold(int OpDet) const
  {
    if(fSelfTrigger_DaphneThresholdMap.find(OpDet)==fSelfTrigger_DaphneThresholdMap.end()) return fSelfTrigger_DaphneThreshold;
    else return fSelfTrigger_DaphneThresholdMap.at(OpDet);
  }
}
