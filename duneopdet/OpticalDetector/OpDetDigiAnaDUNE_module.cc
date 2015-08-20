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
#include "art/Persistency/Common/Ptr.h"
#include "fhiclcpp/ParameterSet.h"

// LArSoft includes

#include "Utilities/DetectorProperties.h"
#include "Utilities/TimeService.h"
#include "RawData/OpDetWaveform.h"

// ROOT includes

#include "TH1.h"

// C++ includes

#include <vector>
#include <map>
#include <cstring>

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
      std::string fInputModule; // Module used to create OpDetWaveforms
      std::string fInstanceName;// Input tag for OpDetWaveforms collection
      float fSampleFreq;        // Sampling frequency in MHz 
      float fTimeBegin;         // Beginning of sample in us

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
    fInputModule  = pset.get< std::string >("InputModule");
    fInstanceName = pset.get< std::string >("InstanceName");

    // Obtain parameters from TimeService
    art::ServiceHandle< util::TimeService > timeService;
    fSampleFreq = timeService->OpticalClock().Frequency();

    // Assume starting at 0
    fTimeBegin  = 0;

  }

  //---------------------------------------------------------------------------
  // Destructor
  OpDetDigiAnaDUNE::~OpDetDigiAnaDUNE()
  {
  }

  //---------------------------------------------------------------------------
  void OpDetDigiAnaDUNE::analyze(art::Event const& evt)
  {

    // Create a string for histogram names
    char histName[50];

    // Map to store how many waveforms are on one optical channel
    std::map< int, int > mapChannelWF;

    // Get OpDetWaveforms from the event
    art::Handle< std::vector< raw::OpDetWaveform > > waveformHandle;
    evt.getByLabel(fInputModule, fInstanceName, waveformHandle);

    // Access ART's TFileService, which will handle creating and writing
    // histograms for us
    art::ServiceHandle< art::TFileService > tfs;

    for (size_t i = 0; i < waveformHandle->size(); i++)
    {

      // This was probably required to overcome the "const" problem 
      // with OpDetPulse::Waveform()
      art::Ptr< raw::OpDetWaveform > waveformPtr(waveformHandle, i);
      raw::OpDetWaveform pulse = *waveformPtr;
      int channel = pulse.ChannelNumber();
      // Make a name for the histogram
      sprintf(histName, "event_%i_opchannel_%i_waveform_%i", 
                        evt.id().event(), channel, mapChannelWF[channel]);
      // Increase counter for number of waveforms on this optical channel
      mapChannelWF[channel]++;

      TH1F *waveformHist = nullptr;

      // Implement different end time for waveforms of variable length
      float endTime = float((pulse.size())/fSampleFreq) + pulse.TimeStamp();

      waveformHist = tfs->make< TH1F >(histName, ";t (us);",
                                 pulse.size(), pulse.TimeStamp(), endTime);

      for (size_t tick = 0; tick < pulse.size(); tick++)
        waveformHist->SetBinContent(tick + 1, (float) pulse[tick]);

    }

  }

}
