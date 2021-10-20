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
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes

#include "larcore/CoreUtils/ServiceUtil.h"
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
        int    fBaselineSubtract;  // Baseline to subtract from each waveform
        int    fNticks;

        // Map to store how many waveforms are on one optical channel
        std::map< int, TH1D* > averageWaveforms;
        std::map< int, int   > waveformCount;
        TH1D* averageWaveformAll;
        int eventCount;


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
        , fInputModule      { pset.get< std::string >("InputModule") }
        , fInstanceName     { pset.get< std::string >("InstanceName") }
        , fTimeBegin        { 0 }
        , fBaselineSubtract { pset.get< int >("BaselineSubtract", 0) }
        , fNticks           { pset.get< int >("Nticks", 0) }
        , averageWaveformAll{ NULL }
        , eventCount        { 0 }
    {

        // Obtain parameters from TimeService
        auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
        fSampleFreq = clockData.OpticalClock().Frequency();

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
        double waveformCountAll = 0;
        for (auto iter = averageWaveforms.begin(); iter != averageWaveforms.end(); iter++)
        {
            waveformCountAll += waveformCount[iter->first];
            mf::LogInfo("AverageWaveform") << "Scaling down channel " << iter->first << " by 1./" << waveformCount[iter->first] << std::endl;
            iter->second->Scale(1./waveformCount[iter->first]);
        }
        mf::LogInfo("AverageWaveform") << "Scaling down all channels by 1./" << eventCount << std::endl;
        averageWaveformAll->Scale(1./eventCount);
    }


    //---------------------------------------------------------------------------
    void AverageWaveform::analyze(art::Event const& evt)
    {

        // Get OpDetWaveforms from the event
        art::InputTag itag1(fInputModule, fInstanceName);
        auto waveformHandle = evt.getHandle< std::vector< raw::OpDetWaveform > >(itag1);

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

            if (fNticks == 0) fNticks = pulse.size();


            // Create the TH1 if it doesn't exist
            auto waveform = averageWaveforms.find( channel );
            if ( waveform == averageWaveforms.end() ) {
                TString histName = TString::Format("avgwaveform_channel_%03i", channel);
                averageWaveforms[channel] =  tfs->make< TH1D >(histName, ";t (us);", fNticks, 0, float(fNticks+1) / fSampleFreq);
            }
            if (!averageWaveformAll) {
                averageWaveformAll =  tfs->make< TH1D >("avgwaveform_channel_all", ";t (us);", fNticks, 0, float(fNticks+1) / fSampleFreq);
            }

            // Add this waveform to this histogram
            for (size_t tick = 0; tick < pulse.size(); tick++) {
                averageWaveforms[channel]->Fill(double(tick)/fSampleFreq, pulse[tick] - fBaselineSubtract);
                averageWaveformAll->Fill(double(tick)/fSampleFreq, pulse[tick] - fBaselineSubtract);
            }

            // Count number of waveforms on each channel
            waveformCount[channel]++;
        }
        eventCount++;

    }

}
