//===========================================================================
// RawWaveformAna_module.cc
// This analyzer module for raw data creates simple ROOT TTree
// with information from raw::OpDetWaveform
//
// @author     : Viktor Pec
// @created    : Oct, 2025
//===========================================================================

#ifndef RawWaveformAna_h
#define RawWaveformAna_h

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpWaveform.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "detdataformats/trigger/TriggerCandidateData2.hpp"
#include "detdataformats/trigger/TriggerCandidateData.hpp"


// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT includes
#include "TH1.h"
#include "THStack.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TTree.h"

// C++ Includes
#include <map>
#include <vector>
#include <iostream>
#include <cstring>
#include <sstream>
#include "math.h"
#include <climits>

namespace opdet {
    class RawWaveformAna : public art::EDAnalyzer{
    public:
        // Standard constructor and destructor for an ART module.
        RawWaveformAna(const fhicl::ParameterSet&);
        virtual ~RawWaveformAna();
        // The analyzer routine, called once per event.
        void analyze (const art::Event&);

    private:
        unsigned int GetTrigType(const art::Event& evt);
        void GetTrigTime(const art::Event& evt, uint64_t &ts_uint64, double &ts_double);


        // The parameters we'll read from the .fcl file.
        std::string fInputModuleLabel;          // Input tag for OpDetWaveform collection
        std::string fTriggerModuleLabel;          // Input tag for raw::Trigger
        std::string fTriggerDataFormat;          // Which trigger data format to use. PDHD uses namespace dunedaq::trgdataformats, PDVD uses dunedaq::trgdataformats2

        TTree* fWaveformTree;

        int Run;
        int SubRun;
        int Event;
        uint64_t TimeStamp_uint64;
        double TimeStamp_double;
        int OpChannel;
        unsigned int nSamples;
        //double SampleSize;
        std::vector<short> adc_value;
        //unsigned int nOpDet;
        unsigned int TriggerType;
        uint64_t TriggerTime_uint64;
        double TriggerTime_double;
    };
}

#endif
namespace opdet {
    DEFINE_ART_MODULE(RawWaveformAna)
}

namespace opdet {
    //-----------------------------------------------------------------------
    // Constructor
    RawWaveformAna::RawWaveformAna(fhicl::ParameterSet const& pset)
        : EDAnalyzer(pset){
        // Read the fcl-file
        fInputModuleLabel  =   pset.get< std::string >("InputModule");
        fTriggerModuleLabel  =   pset.get< std::string >("TriggerModule");
        fTriggerDataFormat  =   pset.get< std::string >("TriggerDataFormat");

        art::ServiceHandle< art::TFileService > tfs;


        fWaveformTree = tfs->make<TTree>("WaveformTree","Waveforms Tree");
        fWaveformTree->Branch("Run"       , &Run       , "Run/I"       );
        fWaveformTree->Branch("SubRun"    , &SubRun    , "SubRun/I"    );
        fWaveformTree->Branch("Event"     , &Event     , "Event/I"     );
        fWaveformTree->Branch("Trigger"     , &TriggerType     , "Trigger/I"     );
        fWaveformTree->Branch("TriggerTime_uint64" , &TriggerTime_uint64, "TriggerTime_uint64/l" );
        fWaveformTree->Branch("TriggerTime_double" , &TriggerTime_double, "TriggerTime_double/D" );
        fWaveformTree->Branch("TimeStamp_uint64" , &TimeStamp_uint64     , "TimeStamp_uint64/l"     );
        fWaveformTree->Branch("TimeStamp_double" , &TimeStamp_double     , "TimeStamp_double/D"     );
        fWaveformTree->Branch("NSamples"     , &nSamples     , "NSamples/I"     );
        fWaveformTree->Branch("OpChannel"     , &OpChannel     , "OpChannel/I"     );
        fWaveformTree->Branch("adc", &adc_value);


    }

    //-----------------------------------------------------------------------
    // Destructor
    RawWaveformAna::~RawWaveformAna() {
    }

    //-----------------------------------------------------------------------
    void RawWaveformAna::analyze(const art::Event& evt){
        auto mfi = mf::LogInfo("RawWaveformAna::analyze()");


        // Map to store how many waveforms are on one optical channel
        std::map< int, int > mapChannelWF;

        // Get deconvolved "OpWaveforms" from the event
        art::Handle< std::vector< raw::OpDetWaveform >> wfmHndl;
        evt.getByLabel(fInputModuleLabel, wfmHndl);

        Run = evt.id().run();
        SubRun = evt.id().subRun();
        Event = evt.id().event();

        TriggerType = GetTrigType(evt);

        GetTrigTime(evt, TriggerTime_uint64, TriggerTime_double);

        for (auto &wfm: *wfmHndl) {
            OpChannel = wfm.ChannelNumber();
            nSamples  = wfm.Waveform().size();

            // Waveform timestamp convetion is dependending on the decoder.
            // For PD HD and VD, the decoder stored this in ticks (16 ns, PDS clock) since the epoch. However, the type double cannot hold its precision.
            // The latest fix to the PD VD digitizer trims the original raw timestamp (uint64) to 40 least significant bits and converts it to a double in microseconds.
            // The uint64 form here stays for backward compatibility. Use the double version if the time stamp is in microseconds.
            TimeStamp_uint64 = (uint64_t)wfm.TimeStamp();
            TimeStamp_double = wfm.TimeStamp(); 
            
            adc_value.resize(nSamples);

            for (unsigned int ii = 0; ii< nSamples ; ii++) {
                adc_value[ii] = wfm.Waveform()[ii];
            }


            fWaveformTree->Fill();
        }


    }


    unsigned int RawWaveformAna::GetTrigType(const art::Event& evt) {
        if (fTriggerDataFormat == "PDVD") {
            auto trigHndl = evt.getHandle< std::vector<dunedaq::trgdataformats2::TriggerCandidateData>>(fTriggerModuleLabel);
            auto &trig = trigHndl->at(0);

            return (unsigned int)trig.type;
        } else if (fTriggerDataFormat == "PDHD") {
            auto trigHndl = evt.getHandle< std::vector<dunedaq::trgdataformats::TriggerCandidateData>>(fTriggerModuleLabel);
            auto &trig = trigHndl->at(0);

            return (unsigned int)trig.type;
        }

        return 0;
    }

    void RawWaveformAna::GetTrigTime(const art::Event& evt, uint64_t &ts_uint64, double &ts_double) {
        if (fTriggerDataFormat == "PDVD") {
            auto trigHndl = evt.getHandle<
                std::vector<dunedaq::trgdataformats2::TriggerCandidateData>
                >(fTriggerModuleLabel);
            auto &trig = trigHndl->at(0);

            // TriggerCandidateData::time_candidate ... time of the trigger signal?
            // TriggerCandidateData::time_start ... time of the DAQ window opened?
            ts_uint64 = trig.time_candidate; // in ticks (16 ns, time system clock) since epoch
            ts_double = (double)trig.time_candidate;

        } else if (fTriggerDataFormat == "PDHD") {
            auto trigHndl = evt.getHandle<
                std::vector<dunedaq::trgdataformats::TriggerCandidateData>
                >(fTriggerModuleLabel);
            auto &trig = trigHndl->at(0);

            ts_uint64 = trig.time_candidate; // in ticks (16 ns, time system clock) since epoch
            ts_double = (double)trig.time_candidate;
        }

        return;
    }
} // namespace opdet
