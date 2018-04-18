//==========================================================
// AverageWaveform_module.cc
// This module finds the average waveform on each channel.
//
// Alex Himmel, ahimmel@fnal.gov
// Based on OpDigiAna_module.cc
//==========================================================

#ifndef AVERAGEWAVEFORM_H
#define AVERAGEWAVEFORM_H 1

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

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RawData/OpDetWaveform.h"

// ROOT includes

#include "TH1.h"

// C++ includes

#include <vector>
#include <map>
#include <cstring>

namespace opdet {

    class AverageWaveform : public art::EDAnalyzer {

    public:

        // Standard constructor and destructor for an ART module
        AverageWaveform(fhicl::ParameterSet const&);
        virtual ~AverageWaveform();

        // The analyzer routine, called once per event
        void analyze(art::Event const&) override;
   
    private:
        void beginJob() override;
        void endJob  () override;

        // Parameters we'll read from the fcl-file
        std::string fInputModule; // Module used to create OpDetWaveforms
        std::string fInstanceName;// Input tag for OpDetWaveforms collection
        double fSampleFreq;        // Sampling frequency in MHz 
        double fTimeBegin;         // Beginning of sample in us

        // Map to store how many waveforms are on one optical channel
        std::map< int, TH1D* > averageWaveforms;
        std::map< int, int   > waveformCount;


    };

}

#endif 

namespace opdet {

    DEFINE_ART_MODULE(AverageWaveform)

}

namespace opdet {

    //---------------------------------------------------------------------------
    // Constructor
    AverageWaveform::AverageWaveform(fhicl::ParameterSet const& pset)
        : EDAnalyzer(pset)
    {

        // Read the fcl-file
        fInputModule  = pset.get< std::string >("InputModule");
        fInstanceName = pset.get< std::string >("InstanceName");

        // Obtain parameters from TimeService
        auto const* timeService = lar::providerFrom<detinfo::DetectorClocksService>();
        fSampleFreq = timeService->OpticalClock().Frequency();

        // Assume starting at 0
        fTimeBegin  = 0;

    }

    //---------------------------------------------------------------------------
    // Destructor
    AverageWaveform::~AverageWaveform()
    {
    }

    
    //---------------------------------------------------------------------------
    void AverageWaveform::beginJob()
    {
    }

    //---------------------------------------------------------------------------
    void AverageWaveform::endJob()
    {
        for (auto iter = averageWaveforms.begin(); iter != averageWaveforms.end(); iter++)
        {
            mf::LogInfo("Scaling down channel ") << iter->first << " by 1./" << waveformCount[iter->first] << std::endl;
            iter->second->Scale(1./waveformCount[iter->first]);
        }
    }


    //---------------------------------------------------------------------------
    void AverageWaveform::analyze(art::Event const& evt)
    {


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

            // Create the TH1 if it doesn't exist
            auto waveform = averageWaveforms.find( channel );
            if ( waveform == averageWaveforms.end() ) {
                TString histName = TString::Format("avgwaveform_channel_%03i", channel);
                averageWaveforms[channel] =  tfs->make< TH1D >(histName, ";t (us);", pulse.size(), 0, double(pulse.size()) / fSampleFreq);
            }

            // Add this waveform to this histogram
            for (size_t tick = 0; tick < pulse.size(); tick++) {
                averageWaveforms[channel]->Fill(double(tick)/fSampleFreq, pulse[tick]);
            }

            // Count number of waveforms on each channel
            waveformCount[channel]++;
        }

    }

}
