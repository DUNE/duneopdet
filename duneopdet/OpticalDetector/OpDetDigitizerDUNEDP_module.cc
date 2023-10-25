//=========================================================
// OpDetDigitizerDUNEDP_module.cc
// This module produces digitized output
// (creating OpDetWaveform)
// from photon detectors taking SimPhotonsLite as input.
//
// J. Soto, CIEMAT, 2018
// Based on OpMCDigi_module.cc
//=========================================================

#ifndef OpDetDigitizerDUNEDP_h
#define OpDetDigitizerDUNEDP_h 1

// Framework includes

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/Exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "CLHEP/Random/RandFlat.h"


// ART extensions
#include "nurandom/RandomUtils/NuRandomService.h"

// LArSoft includes

#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larana/OpticalDetector/OpDetResponseInterface.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoSiPM.h"
#include "duneopdet/OpticalDetector/AlgoSSPLeadingEdge.h"
#include "larana/OpticalDetector/OpDigiProperties.h"

// CLHEP includes

#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandExponential.h"
#include "CLHEP/Random/RandFlat.h"

// C++ includes

#include <vector>
#include <map>
#include <cmath>
#include <memory>

// ROOT includes

#include "TTree.h"


namespace opdet {

  class FocusList
  {
  public:
    FocusList(int nSamples, int padding)
      : fNSamples(nSamples), fPadding(padding) {}

    //Warning, this FocusList feature is actually not used to split the waveform,
    //the use a standard hit finding algorithm instead.
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

  class OpDetDigitizerDUNEDP : public art::EDProducer{

    public:

      OpDetDigitizerDUNEDP(fhicl::ParameterSet const&);
      // Should the destructor be empty?
//      virtual ~OpDetDigitizerDUNEDP();

      void produce(art::Event&);

    private:

      // The parameters read from the FHiCL file
      std::vector < std::string > fInputModule;   // Input tag for OpDet collection
      double  fSampleFreq;         // Sampling frequency in MHz
      double  fTimeBegin;          // Beginning of waveform in us
      double  fTimeEnd;            // End of waveform in us
      double  fVoltageToADC;       // Conversion factor mV to ADC counts
      double  fLineNoiseRMS;       // Pedestal RMS in ADC counts
      double  fDarkNoiseRate;      // In Hz
      double  fCrossTalk;          // Probability of SiPM producing 2 PE signal
                                  // in response to 1 photon
      double  fPedestal;           // In ADC counts
      double  fGain;               // Gain of the PMT
      bool   fDefaultSimWindow;   // Set the start time to -1 drift window and
                                  // the end time to the end time
                                  // of the TPC readout
      bool   fFullWaveformOutput; // Output full waveforms -- produces large
                                  // output. Mostly for debug purposes
      size_t fReadoutWindow;      // In ticks
      size_t fPreTrigger;         // In ticks

      int    fPadding;            // In ticks

      bool   fDigiTree_SSP_LED;   // To create a analysis Tree for SSP LED

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
                    FocusList& fl, double Gain) const;
      double fQE;
      bool fCustomQE;
      std::vector<double> fQEVector;
      bool fCustomGain;
      std::vector<double> fGainVector;
      double  getQE(int OpDet);
      double  getGain(int OpDet);
      double fGainsError;
      double fSPEWidth;
      // Make sure the FHiCL parameters make sense
      void CheckFHiCLParameters() const;

      std::vector< double > fSinglePEWaveform;
      void CreateSinglePEWaveform();
      bool fNegativeSignal; // negative signal if true (as real pmts)

      // Produce waveform on one of the optical detectors
      void CreatePDWaveform(sim::SimPhotonsLite const&,
                            opdet::OpDetResponseInterface const&,
                            geo::Geometry const&,
                            std::vector< std::vector< double > >&,
                            std::vector<FocusList>&);

      // Produce waveform on one of the optical detectors
      void CreatePDWaveform(sim::SimPhotons const&,
                            opdet::OpDetResponseInterface const&,
                            geo::Geometry const&,
                            std::vector< std::vector< double > >&,
                            std::vector<FocusList>&);

      // Vary the pedestal
      void AddLineNoise(std::vector< std::vector< double > >&,
                        const std::vector<FocusList>& fls) const;

      void AddDarkNoise(std::vector< std::vector< double > >&,
                        std::vector<FocusList>& fls,  int opDet);

      unsigned short CrossTalk() const;

      // Create a vector of shorts from a vector of doubles
      // rounding it properly
      std::vector< short > VectorOfDoublesToVectorOfShorts
                                           (std::vector< double > const&) const;

      // Make several shorter waveforms out of a long one using a hit finder,
      // recording also when they start in the long waveform
      std::map< size_t, std::vector< short > >
      SplitWaveform(std::vector< short > const&,
                    const FocusList&);

      double GetDriftWindow(detinfo::DetectorPropertiesData const& detProp) const;

      // Convert time to ticks or the other way around
      // without any checks
      double  TickToTime(size_t tick) const;
      size_t TimeToTick(double  time) const;

      int PMTSaturationFunction(int);
      double SumOfElements(std::vector<double>);
      bool fExportWaveformTree;
      TTree *fWaveformTree;
      TTree *fRunInfo;

      int Run;
      int SubRun;
      int Event;
      unsigned int nSamples;
      double SampleSize;
      std::vector<short*> adc_value;
      unsigned int nOpDet;
  };

}

#endif

namespace opdet {

  DEFINE_ART_MODULE(OpDetDigitizerDUNEDP)

}

namespace opdet {

  //---------------------------------------------------------------------------
  // Constructor
  OpDetDigitizerDUNEDP::OpDetDigitizerDUNEDP(fhicl::ParameterSet const& pset)
  : EDProducer{pset}
  {

    // This module produces (infrastructure piece)
    produces< std::vector< raw::OpDetWaveform > >();

    // Read the fcl-file
    fInputModule = pset.get<std::vector<std::string>>("InputModule",{"largeant"});
    fVoltageToADC       = pset.get< double >("VoltageToADC"      );
    fLineNoiseRMS       = pset.get< double >("LineNoiseRMS"      );
    fCrossTalk          = pset.get< double >("CrossTalk"         );
    fPedestal           = pset.get< double  >("Pedestal"          );
    fDefaultSimWindow   = pset.get< bool   >("DefaultSimWindow"  );
    fFullWaveformOutput = pset.get< bool   >("FullWaveformOutput");
    fReadoutWindow      = pset.get< size_t >("ReadoutWindow"     );
    fPreTrigger         = pset.get< size_t >("PreTrigger"        );
    fPadding            = pset.get< int    >("Padding"           );
    fGain               = pset.get< double >("Gain"              );
    fDigiTree_SSP_LED   = pset.get< bool   >("SSP_LED_DigiTree"  );
    fNegativeSignal     = pset.get< bool   >("NegativeSignal"    );
    fCustomQE           = pset.get< bool   >("CustomQEperOpDet",false    );
    fCustomGain         = pset.get< bool   >("CustomGainperOpDet",false    );
    fGainsError         = pset.get< double >("GainsError",0.0    ); //Factor of gain smearing in units of the gain value
    fSPEWidth           = pset.get< double >("SPEWidth",0.0    ); //Smear the SPE amplitude by the SPE width, in units of the gain value.
    fExportWaveformTree = pset.get< bool   >("ExportWaveformTree", false);

    if(fCustomQE) fQEVector = pset.get<std::vector<double>>("QEVector");
    if(fCustomGain) fGainVector = pset.get<std::vector<double>>("GainVector");

    fThreshAlg = std::make_unique< pmtana::AlgoSSPLeadingEdge >
                   (pset.get< fhicl::ParameterSet >("algo_threshold"));

    if(fDigiTree_SSP_LED){
    	art::ServiceHandle< art::TFileService > tfs;
    	arvore2 = tfs->make<TTree>("PhotonData", "Photon_analysis");
    	arvore2->Branch("photon_opCh",&op_photon);
    	arvore2->Branch("photon_pulse",&t_photon);    
    }
//
    art::ServiceHandle<OpDigiProperties> odp;

    // Obtaining parameters from the DetectorClocksService
    auto const clockData = art::ServiceHandle< detinfo::DetectorClocksService const>()->DataForJob();
    fSampleFreq    = clockData.OpticalClock().Frequency();
    fDarkNoiseRate =  odp->DarkRate();
    fQE = odp->QE();

//    fSampleFreq = odp->SampleFreq(); //Sample Freq in MHz

    if (fDefaultSimWindow)
    {
      auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob(clockData);
      // Assume the readout starts at -1 drift window
      fTimeBegin = -1*GetDriftWindow(detProp);

      // Take the TPC readout window size and convert
      // to us with the electronics clock frequency
      fTimeEnd   = detProp.ReadOutWindowSize() / clockData.TPCClock().Frequency();

      fPreTrigger = 0;//Since we are using negative times as the pretrigger, PreTrigger is zero.
      //fPreTrigger=GetDriftWindow()*clockData.OpticalClock().Frequency();
      fReadoutWindow = (fTimeEnd- fTimeBegin)*clockData.OpticalClock().Frequency();
    }
    else
    {
      fTimeBegin = odp->TimeBegin();
      fTimeEnd   = odp->TimeEnd();

    }

    CheckFHiCLParameters();

    // Initializing random number engines
    unsigned int seed = pset.get< unsigned int >("Seed", sim::GetRandomNumberSeed());
    auto& engine = createEngine(seed);
    fRandGauss       = std::make_unique< CLHEP::RandGauss       >(engine);
    fRandExponential = std::make_unique< CLHEP::RandExponential >(engine);
    fRandFlat        = std::make_unique< CLHEP::RandFlat        >(engine);


    art::ServiceHandle< geo::Geometry > geometry;
    nOpDet=geometry->NOpDets();

    if(fGainsError!=0)
    {
      std::cout << "Smearing gains vector by " << fGainsError << " percent. " << std::endl;
      std::vector<double> newGains;
      for (unsigned int pm=0; pm<nOpDet; pm++) newGains.push_back(fRandGauss->fire(getGain(pm), fGainsError*getGain(pm)));
      std::cout << "Old gains vector: "; for (unsigned int pm=0; pm<nOpDet; pm++) std::cout << getGain(pm) << " "; std::cout << std::endl;
      fCustomGain=true; fGainVector=newGains;
      std::cout << "New gains vector: "; for (unsigned int pm=0; pm<nOpDet; pm++) std::cout << getGain(pm) << " "; std::cout << std::endl;
    }

    fSinglePEWaveform = odp->SinglePEWaveform();
    if(fSampleFreq!=250)
    {
       int rebin=(int)(1.0*250/fSampleFreq);
       std::vector <double> aux;

       std::cout << "\tbefore: SinglePEWaveform: " << fSinglePEWaveform.size() <<"bins, rebin: "<< rebin<<std::endl;
        if (rebin<(int)fSinglePEWaveform.size())
        {
           for(unsigned int i=0; i<fSinglePEWaveform.size();i++)
           {
              if(i%rebin==0) aux.push_back(fSinglePEWaveform[i]);
              else aux[aux.size()-1]+=fSinglePEWaveform[i];
           }
        }
        else
        {
           aux.resize(1);aux[0]=1.0;
        }
        fSinglePEWaveform=aux;
    }
    double normalization =0;
    for(unsigned int i=0; i<fSinglePEWaveform.size();i++) normalization+=fSinglePEWaveform[i];
    std::cout << "Generating waveforms of " << fTimeEnd - fTimeBegin << "us = "<< (fTimeEnd - fTimeBegin)*fSampleFreq <<" Samples"<< std::endl;
    std::cout << "\tTimeBegin" << fTimeBegin <<" "<< std::endl;
    std::cout << "\tfTimeEnd" << fTimeEnd <<" "<< std::endl;
    std::cout << "\tSampleFreq" << fSampleFreq <<" "<< std::endl;
    std::cout << "\tReadoutWindow" << fReadoutWindow <<" "<< std::endl;
    std::cout << "\tfPreTrigger" << fPreTrigger <<" "<< std::endl;
    std::cout << "\fSinglePEWaveform: " << fSinglePEWaveform.size() <<"bins, norm: "<< normalization<<std::endl;

    nSamples = (fTimeEnd - fTimeBegin)*fSampleFreq;
    SampleSize = 1000.0/fSampleFreq;

    art::ServiceHandle<art::TFileService> tfs;

    if(fExportWaveformTree)
    {
      adc_value.resize(nOpDet);
      for (unsigned int i=0; i<nOpDet; i++)adc_value[i] = (short*)malloc(sizeof(short)*nSamples);

      double *vaux; vaux = (double*)malloc(sizeof(double)*nOpDet);
      fRunInfo = tfs->make<TTree>("RunInfo","MonteCarlo Run Info");
      fRunInfo->Branch("SampleSize"   , &SampleSize   , "SampleSize/D"     );
      fRunInfo->Branch("nOpDet"        , &nOpDet   , Form("nOpDet/I")    );
      fRunInfo->Branch("Gains"        , vaux   , Form("Gains[nOpDet]/D")    );
      for (unsigned int i=0; i<nOpDet; i++) vaux[i]=getGain(i);
      for (unsigned int i=0; i<nOpDet; i++) std::cout << " Gain [ " << i << " ] = " << vaux[i] << " ADCxticks/PE" << std::endl;
      fRunInfo->Fill();


      fWaveformTree = tfs->make<TTree>("WaveformTree","Waveforms Tree");
      fWaveformTree->Branch("Run"       , &Run       , "Run/I"       );
      fWaveformTree->Branch("SubRun"    , &SubRun    , "SubRun/I"    );
      fWaveformTree->Branch("Event"     , &Event     , "Event/I"     );
      fWaveformTree->Branch("NSamples"     , &nSamples     , "NSamples/I"     );
      for (unsigned int i=0; i<nOpDet; i++) fWaveformTree->Branch(Form("adc_channel_%i",i), adc_value[i] , Form("adc_value_%i[NSamples]/S",i)     );

    }

  }

  //---------------------------------------------------------------------------
  void OpDetDigitizerDUNEDP::produce(art::Event& evt)
  {
    if(fExportWaveformTree)
    {
      Run    = evt.run();
      SubRun = evt.subRun();
      Event  = evt.event();
    }
    if(fDigiTree_SSP_LED) art::ServiceHandle< art::TFileService > tfs;

    // A pointer that will store produced OpDetWaveforms
    std::unique_ptr< std::vector< raw::OpDetWaveform > >
    pulseVecPtr(std::make_unique< std::vector< raw::OpDetWaveform > >());

    art::ServiceHandle< sim::LArG4Parameters > lgp;
    bool fUseLitePhotons = lgp->UseLitePhotons();

    if (!fUseLitePhotons)
      throw art::Exception(art::errors::UnimplementedFeature)
        << "Sorry, but for now only Lite Photon digitization is implemented!"
        << '\n';

    // Geometry service
    art::ServiceHandle< geo::Geometry > geometry;

    // Service for determining optical detector responses
    art::ServiceHandle< opdet::OpDetResponseInterface > odResponse;

//    std::vector<art::Handle< std::vector< sim::SimPhotonsLite > > > litePhotonHandle;

    int modulecounter=0;


//     unsigned int nChannelsPerOpDet=1;

     unsigned int nChannelsPerOpDet = geometry->NOpHardwareChannels(nOpDet);


   // This vector stores waveforms created for each optical channel
    
   for(unsigned int opDet=0; opDet<nOpDet; opDet++)
   {

    std::vector< std::vector< double > > pdWaveforms(nChannelsPerOpDet, std::vector< double >(nSamples, static_cast< double >(fPedestal)));
    std::vector<std::vector<FocusList>> fls(nOpDet,std::vector<FocusList>(nChannelsPerOpDet, FocusList(nSamples, fPadding)));

    for(auto mod : fInputModule){
      // For every module collection
      //std::cout << "ANALYZING " << mod<<std::endl;                                       
      modulecounter++;

      art::Handle< std::vector< sim::SimPhotonsLite > > litePhotonHandle;
      if(fUseLitePhotons) litePhotonHandle = evt.getHandle< std::vector< sim::SimPhotonsLite > >(mod);

      art::Handle< std::vector< sim::SimPhotons > > PhotonHandle;
      if(!fUseLitePhotons) PhotonHandle = evt.getHandle< std::vector< sim::SimPhotons > >(mod);

      if(fUseLitePhotons)
      {

       // For every optical detector:
       for (auto const& litePhotons : (*litePhotonHandle))
       {
          // OpChannel in SimPhotonsLite is actually the photon detector number
          //unsigned int opDet = litePhotons.OpChannel;
         if(opDet == (unsigned)litePhotons.OpChannel)
         {
          // Get number of channels in this optical detector
          //unsigned int nChannelsPerOpDet = geometry->NOpHardwareChannels(opDet);
  	  // it is one by default in the dual phase geometry, but let's keep it to have compatible functions with single phase.

          CreatePDWaveform(litePhotons, *odResponse, *geometry, pdWaveforms, fls[opDet]);
          if((unsigned)modulecounter<fInputModule.size()) continue;//==fInputModule.size()

          if(fExportWaveformTree) for (unsigned int hardwareChannel = 0;
             hardwareChannel < nChannelsPerOpDet; ++hardwareChannel) fls[opDet][hardwareChannel].AddRange(0,SampleSize);


          // Generate dark noise
          if (fDarkNoiseRate > 0.0) AddDarkNoise(pdWaveforms, fls[opDet], opDet);

          // Uncomment to undo the effect of FocusLists. Replaces the accumulated
          // lists with ones asserting we need to look at the whole trace.
          // for(FocusList& fl: fls){
          //        fl.ranges.clear();
          //        fl.ranges.emplace_back(0, nSamples-1);
          // }

          // Vary the pedestal

          if (fLineNoiseRMS > 0.0)  AddLineNoise(pdWaveforms, fls[opDet]);

          // Loop over all the created waveforms, split them into shorter
          // waveforms and use them to initialize OpDetWaveforms

          for (unsigned int hardwareChannel = 0;
             hardwareChannel < nChannelsPerOpDet; ++hardwareChannel)
          {
            for(const std::pair<int, int>& p: fls[opDet][hardwareChannel].ranges){
            // It's a shame we copy here. We could actually avoid by making the
            // functions below take a begin()/end() pair.
              const std::vector<double> sub(pdWaveforms[hardwareChannel].begin()+p.first,
                                        pdWaveforms[hardwareChannel].begin()+p.second+1);

              std::vector< short > waveformOfShorts =
                 VectorOfDoublesToVectorOfShorts(sub);
              //std::cout << "waveformOfShorts " << waveformOfShorts.size()<< std::endl;

              std::map< size_t, std::vector < short > > mapTickWaveform =
                (!fFullWaveformOutput) ?
                SplitWaveform(waveformOfShorts, fls[opDet][hardwareChannel]) :
                std::map< size_t, std::vector< short > >{ std::make_pair(0,
                                                                       waveformOfShorts) };

              //std::cout << "mapTickWaveform " << mapTickWaveform.size()<< std::endl;
              unsigned int opChannel = geometry->OpChannel(opDet, hardwareChannel);
              for (auto const& pairTickWaveform : mapTickWaveform)
              {
                double timeStamp =
                static_cast< double >(TickToTime(pairTickWaveform.first+p.first));
                //std::cout << "\tp " << p<< std::endl;
                //std::cout << "\ttimeStamp " << timeStamp<< ", pairTickWaveform.first " << pairTickWaveform.first <<", p.first " << p.first<< std::endl;
                // std::cout << "\tpairTickWaveform.second.size() " << pairTickWaveform.second.size()<< std::endl;

                raw::OpDetWaveform adcVec(timeStamp, opChannel,
                                          pairTickWaveform.second.size());
                int counter=0;
                for (short const& value : pairTickWaveform.second)
                {
                  //std::cout <<"\t\tvalue " << value<< std::endl;
                  adcVec.emplace_back(value); counter++;
	        } 
                //std::cout <<"\t\tcounter " << counter<< std::endl;
       	   
                pulseVecPtr->emplace_back(std::move(adcVec));

	      }//endloop per 
	    }//endloop per Focus List <fls>
          }//endloop per Hardware channel
         }//endif pmt 
        }//endloop per OpDet
      }//endif litephotons
      else //if SimPhotons
      {  
       // For every optical detector:
       for (auto const& Photons : (*PhotonHandle))
       {
         if(opDet == (unsigned)Photons.OpChannel())
         {
          CreatePDWaveform(Photons, *odResponse, *geometry, pdWaveforms, fls[opDet]);
          if((unsigned)modulecounter<fInputModule.size()) continue;//==fInputModule.size()

          if (fDarkNoiseRate > 0.0) AddDarkNoise(pdWaveforms, fls[opDet], opDet);

          if (fLineNoiseRMS > 0.0)  AddLineNoise(pdWaveforms, fls[opDet]);

          // Loop over all the created waveforms, split them into shorter
          // waveforms and use them to initialize OpDetWaveforms

          for (unsigned int hardwareChannel = 0;
             hardwareChannel < nChannelsPerOpDet; ++hardwareChannel)
          {
            for(const std::pair<int, int>& p: fls[opDet][hardwareChannel].ranges){
            // It's a shame we copy here. We could actually avoid by making the
            // functions below take a begin()/end() pair.
              const std::vector<double> sub(pdWaveforms[hardwareChannel].begin()+p.first,
                                        pdWaveforms[hardwareChannel].begin()+p.second+1);

              std::vector< short > waveformOfShorts =
                 VectorOfDoublesToVectorOfShorts(sub);

              std::map< size_t, std::vector < short > > mapTickWaveform =
                (!fFullWaveformOutput) ?
                SplitWaveform(waveformOfShorts, fls[opDet][hardwareChannel]) :
                std::map< size_t, std::vector< short > >{ std::make_pair(0,
                                                                       waveformOfShorts) };

              //std::cout << "mapTickWaveform " << mapTickWaveform.size()<< std::endl;
              unsigned int opChannel = geometry->OpChannel(opDet, hardwareChannel);
              for (auto const& pairTickWaveform : mapTickWaveform)
              {
                double timeStamp =
                static_cast< double >(TickToTime(pairTickWaveform.first+p.first));
                //std::cout << "\tp " << p<< std::endl;
                //std::cout << "\ttimeStamp " << timeStamp<< ", pairTickWaveform.first " << pairTickWaveform.first <<", p.first " << p.first<< std::endl;
                // std::cout << "\tpairTickWaveform.second.size() " << pairTickWaveform.second.size()<< std::endl;

                raw::OpDetWaveform adcVec(timeStamp, opChannel,
                                          pairTickWaveform.second.size());
                int counter=0;
                for (short const& value : pairTickWaveform.second)
                {
                  //std::cout <<"\t\tvalue " << value<< std::endl;
                  adcVec.emplace_back(value); counter++;
	        } 
                //std::cout <<"\t\tcounter " << counter<< std::endl;
       	   
                pulseVecPtr->emplace_back(std::move(adcVec));
	      }//endloop per 
	    }//endloop per Focus List <fls>
          }//endloop per Hardware channel
         }//endif pmt 
        }//endloop per OpDet
      } //endif litePhotons
      if(fUseLitePhotons) litePhotonHandle.clear();
      if(!fUseLitePhotons) PhotonHandle.clear();
    }//endloop per Input Module (S1 and S2 light)

   }//endloop per pmt
    if(fDigiTree_SSP_LED)
    {
      arvore2->Fill();
      t_photon.clear();
      op_photon.clear();
    }

    if(fExportWaveformTree)
    {
       for(unsigned int j=0; j < pulseVecPtr->size(); j++)
       {
         unsigned int opDet= pulseVecPtr->at(j).ChannelNumber();
//         std::cout <<" one - pmt - " <<  j << "\t" << opDet  << "\t"<< pulseVecPtr->at(j).Waveform().size() << " " << pulseVecPtr->at(j).Waveform()[0] <<  std::endl;
         for (unsigned int i = 0; i< pulseVecPtr->at(j).Waveform().size() ; i++) adc_value[opDet][i] = pulseVecPtr->at(j).Waveform()[i];
       }
       fWaveformTree->Fill();

    }

    // Push the OpDetWaveforms into the event
    evt.put(std::move(pulseVecPtr));




  }

  //---------------------------------------------------------------------------
  void OpDetDigitizerDUNEDP::AddPulse(size_t timeBin,
                                    int scale, std::vector< double >& waveform,
                                    FocusList& fl, double Gain) const
  {
	if(timeBin>waveform.size()) return; //photon out of time

    size_t pulseLength = fSinglePEWaveform.size()+1;
    if ((timeBin + fSinglePEWaveform.size()) > waveform.size())
      pulseLength = (waveform.size() - timeBin);

    fl.AddRange(timeBin, timeBin+pulseLength-1);


    double interbinallignment = fRandFlat->fire(1.0);//shifting the deposited charge among consecutive bins.
    double customGain = fRandGauss->fire(Gain, fSPEWidth*Gain); //smear the SPE amplitude by the SPEwidth
    // Adding a pulse to the waveform
    for (size_t tick = 0; tick != pulseLength; ++tick)
    {
      double bincharge;
      if (tick==0) bincharge=interbinallignment*scale*customGain*fSinglePEWaveform[tick];
      else bincharge = interbinallignment*scale*customGain*fSinglePEWaveform[tick] + (1-interbinallignment)*scale*customGain*fSinglePEWaveform[tick-1];

      if(!fNegativeSignal) waveform[timeBin + tick] += bincharge;
      else waveform[timeBin + tick] -= (double)bincharge;

    }
  }

  //---------------------------------------------------------------------------
  void OpDetDigitizerDUNEDP::CreatePDWaveform
                             (sim::SimPhotonsLite const& litePhotons,
                              opdet::OpDetResponseInterface const& odResponse,
                              geo::Geometry const& geometry,
                              std::vector< std::vector< double > >& pdWaveforms,
                              std::vector<FocusList>& fls)
  {

    unsigned int const opDet = litePhotons.OpChannel;
    // This is int because otherwise detectedLite doesn't work
    int readoutChannel;
    // For a group of photons arriving at the same time this is a map
    // of < arrival time (in ns), number of photons >
    std::map< int, int > const& photonsMap = litePhotons.DetectedPhotons;

    int counter=0;
    // For every pair of (arrival time, number of photons) in the map:

    for (auto const& pulse : photonsMap)
    {
      // Converting ns to us
      double photonTime = static_cast< double >(pulse.first)/1000.0;

      int NumberOfPEs = PMTSaturationFunction(pulse.second);
//      std::cout << "Adding " <<NumberOfPEs << "PEs @ " << photonTime <<  " in ch" << opDet << std::endl;
      for (int i = 0; i < NumberOfPEs; ++i)
      {
        if ((photonTime >= fTimeBegin) && (photonTime < fTimeEnd))
        {
          // Sample a random subset according to QE
          if (CLHEP::RandFlat::shoot(1.0) <getQE(opDet))
          {
	    odResponse.detectedLite(opDet, readoutChannel);
            unsigned int hardwareChannel = geometry.HardwareChannelFromOpChannel(readoutChannel);
            // Convert the time of the pulse to ticks
            size_t timeBin = TimeToTick(photonTime);
            // Add 1 pulse to the waveform
            AddPulse(timeBin, CrossTalk(), pdWaveforms.at(hardwareChannel), fls[hardwareChannel], getGain(opDet));
            counter++;
	    unsigned int opChannel = geometry.OpChannel(opDet, hardwareChannel);
	    if(fDigiTree_SSP_LED){
	    	op_photon.emplace_back(opChannel);
	    	t_photon.emplace_back(photonTime);
	    }
	  }
	  //else std::cout << "photon not detected "<< std::endl;
        }
      }
    }
//   std::cout << "Created waveform for channel " << opDet << " with " << counter << " photons." << std::endl;
  
  }

  //---------------------------------------------------------------------------
  void OpDetDigitizerDUNEDP::CreatePDWaveform
                             (sim::SimPhotons const& Photons,
                              opdet::OpDetResponseInterface const& odResponse,
                              geo::Geometry const& geometry,
                              std::vector< std::vector< double > >& pdWaveforms,
                              std::vector<FocusList>& fls)
  {

    unsigned int const opDet = Photons.OpChannel();
    int readoutChannel;
    // For a group of photons arriving at the same time this is a map
    // of < arrival time (in ns), number of photons >
    std::map< int, int > photonsMap;
    for(unsigned int i=0; i<Photons.size();i++)photonsMap[Photons[i].Time]=0;
    for(unsigned int i=0; i<Photons.size();i++)photonsMap[Photons[i].Time]++;
    int counter=0;
    // For every pair of (arrival time, number of photons) in the map:

    for (auto const& pulse : photonsMap)
    {
      // Converting ns to us
      double photonTime = static_cast< double >(pulse.first)/1000.0;

      int NumberOfPEs = PMTSaturationFunction(pulse.second);
//      std::cout << "Adding " <<NumberOfPEs << "PEs @ " << photonTime <<  " in ch" << opDet << std::endl;
      for (int i = 0; i < NumberOfPEs; ++i)
      {
        if ((photonTime >= fTimeBegin) && (photonTime < fTimeEnd))
        {
          // Sample a random subset according to QE
          if (CLHEP::RandFlat::shoot(1.0) <getQE(opDet))
          {
	    odResponse.detectedLite(opDet, readoutChannel);
            unsigned int hardwareChannel = geometry.HardwareChannelFromOpChannel(readoutChannel);
            // Convert the time of the pulse to ticks
            size_t timeBin = TimeToTick(photonTime);
            // Add 1 pulse to the waveform
            AddPulse(timeBin, CrossTalk(), pdWaveforms.at(hardwareChannel), fls[hardwareChannel], getGain(opDet));
            counter++;
	    unsigned int opChannel = geometry.OpChannel(opDet, hardwareChannel);
	    if(fDigiTree_SSP_LED){
	    	op_photon.emplace_back(opChannel);
	    	t_photon.emplace_back(photonTime);
	    }
	  }
	  //else std::cout << "photon not detected "<< std::endl;
        }
      }
    }
//   std::cout << "Created waveform for channel " << opDet << " with " << counter << " photons." << std::endl;
  
  }


  //---------------------------------------------------------------------------
  void OpDetDigitizerDUNEDP::
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
  void OpDetDigitizerDUNEDP::
  AddDarkNoise(std::vector< std::vector< double > >& waveforms,
               std::vector<FocusList>& fls, int optDet)
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
        AddPulse(timeBin, CrossTalk(), waveform, fls[i], getGain(optDet) );
        // Find next time to simulate a single PE pulse
        darkNoiseTime += static_cast< double >
                        (fRandExponential->fire(1.0/fDarkNoiseRate)*1000000.0);
      }

      ++i;
    }
  }

  //---------------------------------------------------------------------------
  unsigned short OpDetDigitizerDUNEDP::CrossTalk() const
  {
    // Sometimes this should produce 3 or more PEs (not implemented)
    if      (fCrossTalk <= 0.0)                 return 1;
    else if (fRandFlat->fire(1.0) > fCrossTalk) return 1;
    else                                        return 2;
  }

  //---------------------------------------------------------------------------
  std::vector< short > OpDetDigitizerDUNEDP::VectorOfDoublesToVectorOfShorts
                            (std::vector< double > const& vectorOfDoubles) const
  {
    // Don't bother to round properly, it's faster this way
    return std::vector<short>(vectorOfDoubles.begin(), vectorOfDoubles.end());

    /*
    std::vector< short > vectorOfShorts;
    vectorOfShorts.reserve(vectorOfDoubles.size());
    int i=0;
    for (short const& value : vectorOfDoubles)
    {
      vectorOfShorts.emplace_back(static_cast< short >(std::round(value)));
    }
    return vectorOfShorts;
 //   */
  }

  //---------------------------------------------------------------------------
  std::map< size_t, std::vector< short > > OpDetDigitizerDUNEDP::
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
  double OpDetDigitizerDUNEDP::GetDriftWindow(detinfo::DetectorPropertiesData const& detProp) const
  {

    double driftWindow;

    double maxDrift = 0.0;
    for (geo::TPCGeo const& tpc :
           art::ServiceHandle< geo::Geometry >()->Iterate<geo::TPCGeo>())
      if (maxDrift < tpc.DriftDistance()) maxDrift = tpc.DriftDistance();

    driftWindow = maxDrift / detProp.DriftVelocity();

    return driftWindow;

  }

  //---------------------------------------------------------------------------
  double OpDetDigitizerDUNEDP::TickToTime(size_t tick) const
  {

//std::cout << "tick "<< tick<< " " << (tick > fPreTrigger) << " " << (static_cast< double >(tick - fPreTrigger)/fSampleFreq
  //                                                               + fTimeBegin) << " " << (static_cast< double >(fPreTrigger - tick)/fSampleFreq*(-1.0)
//                                                                 + fTimeBegin)<< std::endl;
    if (tick > fPreTrigger)
      return (static_cast< double >(tick - fPreTrigger)/fSampleFreq
                                                                 + fTimeBegin);
    else
      return (static_cast< double >(fPreTrigger - tick)/fSampleFreq*(-1.0)
                                                                 + fTimeBegin);

  }

  //---------------------------------------------------------------------------
  size_t OpDetDigitizerDUNEDP::TimeToTick(double time) const
  {

    return static_cast< size_t >(std::round((time - fTimeBegin)*fSampleFreq
                                                               + fPreTrigger));

  }

  //---------------------------------------------------------------------------
  void OpDetDigitizerDUNEDP::CheckFHiCLParameters() const
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

  int  OpDetDigitizerDUNEDP::PMTSaturationFunction(int photons){ return photons;}
  double OpDetDigitizerDUNEDP::SumOfElements( std::vector<double> aa)
  {
    double sum=0;
    for (auto& n : aa)    sum+=n;
    return sum;
  }
  double  OpDetDigitizerDUNEDP::getQE(int OpDet)
  {
    if(!fCustomQE) return fQE;
    else return fQEVector[OpDet];
  }

  double  OpDetDigitizerDUNEDP::getGain(int OpDet)
  {
    if(!fCustomGain) return fGain;
    else return fGainVector[OpDet];
  }

}
