//==========================================================
// OpDetDigiAnaDUNE_module.cc
// This analysis module creates histograms
// with information from OpDetWaveforms
//
// Gleb Sinev, Duke, 2015
// Based on OpDigiAna_module.cc
//==========================================================

#ifndef OpDetDigiAnaDUNE_H
#define OpDetDigiAnaDUNE_H 1

// Framework includes

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "fhiclcpp/ParameterSet.h"

// LArSoft includes

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RawData/OpDetWaveform.h"

// ROOT includes

#include "TH1.h"
#include "TTree.h" //Vitor Luzio

// C++ includes

#include <vector>
#include <map>
#include <string>
#include <sstream>

namespace opdet {

  class OpDetDigiAnaDUNE : public art::EDAnalyzer {

    public:

      // Standard constructor and destructor for an ART module
      OpDetDigiAnaDUNE(fhicl::ParameterSet const&);
      virtual ~OpDetDigiAnaDUNE();

      // The analyzer routine, called once per event
      void analyze(art::Event const&);
   
    private:

      // Parameters we'll read from the fcl-file
      std::string fInputModule;  // Module used to create OpDetWaveforms
      std::string fInstanceName; // Input tag for OpDetWaveforms collection
      double fSampleFreq;        // Sampling frequency in MHz 
      bool fAnaTree_SSPLED;	 // Parameter to write a SSP LED Tree


      std::vector<double> wf_start;
      std::vector<int> canal_op;
      std::vector<double> wf_inicio;
      std::vector<double> wf_end;
      std::vector<int> wf_size;
      TTree *arvore;

  };

}

#endif 

namespace opdet {

  DEFINE_ART_MODULE(OpDetDigiAnaDUNE)

}

namespace opdet {

  //---------------------------------------------------------------------------
  // Constructor
  OpDetDigiAnaDUNE::OpDetDigiAnaDUNE(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
  {

    // Read the fcl-file
    fInputModule  =   pset.get< std::string >("InputModule");
    fInstanceName =   pset.get< std::string >("InstanceName");
    fAnaTree_SSPLED = pset.get< bool        >("SSP_LED_AnaTree"); 
    // Obtain parameters from DetectorClocksService
    auto const *timeService = lar::providerFrom< detinfo::DetectorClocksService >();
    fSampleFreq = timeService->OpticalClock().Frequency();
  
    art::ServiceHandle< art::TFileService > tfs;

//	std::cout<<"Valor bool = "<<fAnaTree_SSPLED<<std::endl;

    if(fAnaTree_SSPLED){
	arvore = tfs->make<TTree>("TriggerData","TriggerInfo");
    	arvore->Branch("OpticalCh",&canal_op);
    	arvore->Branch("WfFirstTime",&wf_inicio);
    	arvore->Branch("WfSize",&wf_size);
    	arvore->Branch("Trg_PreTrg",&wf_start);
    	arvore->Branch("WfEnd",&wf_end);
    }

  }

  //---------------------------------------------------------------------------
  // Destructor
  OpDetDigiAnaDUNE::~OpDetDigiAnaDUNE()
  {
  }

  //---------------------------------------------------------------------------
  void OpDetDigiAnaDUNE::analyze(art::Event const& evt)
  {

    // Map to store how many waveforms are on one optical channel
    std::map< int, int > mapChannelWF;

    // Get OpDetWaveforms from the event
    art::Handle< std::vector< raw::OpDetWaveform > > waveformHandle;
    evt.getByLabel(fInputModule, fInstanceName, waveformHandle);

    // Access ART's TFileService, which will handle creating and writing
    // histograms for us
    art::ServiceHandle< art::TFileService > tfs;


    double firstWaveformTime = -9999;

    for (auto const& waveform : *waveformHandle)
    {
      if (firstWaveformTime < 0) firstWaveformTime = waveform.TimeStamp();
    }


    std::cout << "Event #" << evt.id().event()  << ", firstTime: " << firstWaveformTime << std::endl;
    
    for (auto const& waveform : *waveformHandle)
    {

      int channel = waveform.ChannelNumber();
      // Make a name for the histogram
      std::stringstream histName;
      histName << "event_"      << evt.id().event() 
               << "_opchannel_" << channel
               << "_waveform_"  << mapChannelWF[channel];
      // Increase counter for number of waveforms on this optical channel
      mapChannelWF[channel]++;
 
      // Implement different end time for waveforms of variable length
      double startTime = waveform.TimeStamp() - firstWaveformTime;
      double endTime = double(waveform.size())/fSampleFreq + startTime;

      if(fAnaTree_SSPLED){
      	canal_op.emplace_back(channel);
      	wf_inicio.emplace_back(firstWaveformTime);
      	wf_start.emplace_back(startTime);
      	wf_end.emplace_back(endTime);
      	wf_size.emplace_back(waveform.size());
      }
      // Create a new histogram
      TH1D *waveformHist = tfs->make< TH1D >(histName.str().c_str(),                        TString::Format(";t - %f (#mus);",firstWaveformTime),                           waveform.size(), startTime, endTime);

     // Copy values from the waveform into the histogram
      for (size_t tick = 0; tick < waveform.size(); tick++){
        waveformHist->SetBinContent(tick + 1, waveform[tick]);}
 

    }
	if(fAnaTree_SSPLED){
		arvore->Fill();
        	canal_op.clear();
        	wf_inicio.clear();
        	wf_start.clear();
        	wf_end.clear();
        	wf_size.clear();
	}

  }

}
