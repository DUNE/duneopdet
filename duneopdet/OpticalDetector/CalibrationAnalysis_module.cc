//==========================================================
// CalibrationAnalysis_module.cc
// This module performs analysis for calibration runs.
//
// Jonathan Insler, jti3@fnal.gov
// Based on AverageWaveform_module.cc
//==========================================================

#ifndef CALIBRATIONANALYSIS_H
#define CALIBRATIONANALYSIS_H 1

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

    class CalibrationAnalysis : public art::EDAnalyzer {

    public:

        // Standard constructor and destructor for an ART module
        CalibrationAnalysis(fhicl::ParameterSet const&);
        virtual ~CalibrationAnalysis();

        // The analyzer routine, called once per event
        void analyze(art::Event const&);

   
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

      TH1D*       fPedestalMeanPerChannel;
      TH1D*       fPedestalSigmaPerChannel;
      TH1D*       fIntegratedSignalMeanPerChannel;
      TH1D*       fFractionSamplesNearMaximumPerChannel;
      TH1D*       fNumberOfWaveformsProcessedPerChannel;
      

    };

}

#endif 

namespace opdet {

    DEFINE_ART_MODULE(CalibrationAnalysis)

}

namespace opdet {

    //---------------------------------------------------------------------------
    // Constructor
    CalibrationAnalysis::CalibrationAnalysis(fhicl::ParameterSet const& pset)
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
    CalibrationAnalysis::~CalibrationAnalysis()
    {
    }

    
    //---------------------------------------------------------------------------
    void CalibrationAnalysis::beginJob()
    {
     
      // Access ART's TFileService, which will handle creating and writing
      // histograms for us
      art::ServiceHandle< art::TFileService > tfs;


      fPedestalMeanPerChannel =  tfs->make< TH1D >("PedestalMeanPerChannel", "Pedestal Means;Channel Number;Pedestal Mean (ADC)", 711, 0.5, 711.5);
      fPedestalSigmaPerChannel =  tfs->make< TH1D >("PedestalSigmaPerChannel", "Pedestal Sigma;Channel Number;Pedestal Sigma (ADC)", 711, 0.5, 711.5);
      fIntegratedSignalMeanPerChannel =  tfs->make< TH1D >("IntegratedSignalMeanPerChannel", "Integrated Signal Means;Channel Number; Integrated SignalMean (ADC)", 711, 0.5, 711.5);
      fIntegratedSignalMeanPerChannel->Sumw2();
      fFractionSamplesNearMaximumPerChannel =  tfs->make< TH1D >("FractionSamplesNearMaximumPerChannel", "Fraction of Samples Near Maximum;Channel Number;Fraction Near Maximum", 711, 0.5, 711.5);

      fNumberOfWaveformsProcessedPerChannel =  tfs->make< TH1D >("NumberOfWaveformsProcessedPerChannel", "Number of Waveforms Processed per Channel;Channel Number;Fraction Near Maximum", 711, 0.5, 711.5);

    }

    //---------------------------------------------------------------------------
    void CalibrationAnalysis::endJob()
    {

        for (auto iter = averageWaveforms.begin(); iter != averageWaveforms.end(); iter++)
        {
            mf::LogInfo("Scaling down channel ") << iter->first << " by 1./" << waveformCount[iter->first] << std::endl;
            iter->second->Scale(1./waveformCount[iter->first]);

	    //Also scale down pedestal mean, sigma, and integrated signal mean and sigma histograms

	    double adjustedPedestalMean = fPedestalMeanPerChannel->GetBinContent(iter->first)/waveformCount[iter->first];
	    fPedestalMeanPerChannel->SetBinContent(iter->first,adjustedPedestalMean);
	    double adjustedPedestalSigma = fPedestalSigmaPerChannel->GetBinContent(iter->first)/waveformCount[iter->first];
	    fPedestalSigmaPerChannel->SetBinContent(iter->first,adjustedPedestalSigma);

	    double adjustedFractionSamplesNearMaximumPerChannel = fFractionSamplesNearMaximumPerChannel->GetBinContent(iter->first)/waveformCount[iter->first];
	    fFractionSamplesNearMaximumPerChannel->SetBinContent(iter->first,adjustedFractionSamplesNearMaximumPerChannel);
	   
	    fNumberOfWaveformsProcessedPerChannel->SetBinContent(iter->first,waveformCount[iter->first]);

	    double adjustedIntegratedSignalMeanPerChannel = fIntegratedSignalMeanPerChannel->GetBinContent(iter->first)/waveformCount[iter->first];
	    //double adjustedIntegratedSignalErrorPerChannel = fIntegratedSignalMeanPerChannel->GetBinError(iter->first)/waveformCount[iter->first];

	    fIntegratedSignalMeanPerChannel->SetBinContent(iter->first,adjustedIntegratedSignalMeanPerChannel);
	    //fIntegratedSignalMeanPerChannel->SetBinError(iter->first,adjustedIntegratedSignalErrorPerChannel);

        }

	

    }


    //---------------------------------------------------------------------------
    void CalibrationAnalysis::analyze(art::Event const& evt)
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

	    double PedestalMean = 0;
	    double PedestalVariance = 0;
	    double IntegratedSignalMean = 0;
	    // double IntegratedSignalSigma = 0;
	    

            // Add this waveform to this histogram
            for (size_t tick = 0; tick < pulse.size(); tick++) {
                averageWaveforms[channel]->Fill(double(tick)/fSampleFreq, pulse[tick]);

		if(tick < 30)
		  PedestalMean += pulse[tick];

		if(pulse[tick] > 4000)
		  fFractionSamplesNearMaximumPerChannel->Fill(channel,1.0/pulse.size());
            }

	    PedestalMean /= 30;


	    
	    for (size_t tick = 0; tick < pulse.size(); tick++) {
	      
	      PedestalVariance += ( pulse[tick] - PedestalMean)*( pulse[tick] - PedestalMean);
	      
	      IntegratedSignalMean += pulse[tick] - PedestalMean;
            }
	    

	    PedestalVariance /= (pulse.size() - 1);
	    double PedestalSigma = sqrt(PedestalVariance);
	    
	    fPedestalMeanPerChannel->Fill(channel,PedestalMean);
	    fPedestalSigmaPerChannel->Fill(channel,PedestalSigma);

	    fIntegratedSignalMeanPerChannel->Fill(channel,IntegratedSignalMean);

            // Count number of waveforms on each channel
            waveformCount[channel]++;

	    
        }

    }

}
