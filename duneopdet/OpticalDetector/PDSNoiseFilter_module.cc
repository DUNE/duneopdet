//=========================================================
// OpDetDigitizerDUNE_module.cc
// This module produces digitized output
// (creating OpDetWaveform)
// from photon detectors taking SimPhotonsLite as input.
//
// Gleb Sinev, Duke, 2015
// Based on OpMCDigi_module.cc


// This module was modified to filter the raw waveforms
// from protoDUNE SP using two combined techniques:
// 1) mobile average
// 2) denoising code
//
// Marina Reggiani-Guzzo & Dante Totani, October 2018
//=========================================================


// Framework includes

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/Exception.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


// ART extensions
#include "nurandom/RandomUtils/NuRandomService.h"

// LArSoft includes
#include "lardataobj/RawData/OpDetWaveform.h"

// C++ includes
#include <vector>
#include <string>

// ROOT includes
#include <Rtypes.h> 

using std::string;
using std::vector;



// Define the configuration

namespace opdet {

  class PDSNoiseFilter : public art::EDProducer{

    public:
      // Define the configuration
      struct Config {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;

        fhicl::Atom<string>     InputModule { Name("InputModule"), Comment("Module which produced raw waveforms") };
        fhicl::Sequence<string> InputLabels { Name("InputLabels"), Comment("Labels of the raw waveforms. Will be reused for output.") };
      // Additional fhicl parameters here
      };
      using Parameters = art::EDProducer::Table<Config>;


      PDSNoiseFilter(Parameters const & config);

      void produce(art::Event&);

    private:

      // The parameters read from the FHiCL file
      string fInputModule;   // Input tag for OpDet collection
      vector<string> fInputLabels;   // Input tag for OpDet collection
      // Additional fhicl parameters here

    std::vector<float> filter( std::vector<short unsigned int> vec, Int_t nBins);
    void TV1D_denoise(vector<float> input, vector<float>& secondFilter, const int width, const float lambda);
  };

}


namespace opdet {

  //---------------------------------------------------------------------------
  // Constructor
  PDSNoiseFilter::PDSNoiseFilter(Parameters const & config)
    : EDProducer{config},
      fInputModule(config().InputModule()),
      fInputLabels(config().InputLabels())
      // Additional fhicl parameters here
  {
    // Tell ART what we intend to produce
    for (auto label : fInputLabels) {
      produces< vector< raw::OpDetWaveform > >(label);
    }
  }

  //---------------------------------------------------------------------------
  void PDSNoiseFilter::produce(art::Event& evt)
  {

    // Loop over the input labels
    for (auto label : fInputLabels) {

      // Read the waveforms in from the event record
      art::InputTag itag1(fInputModule, label);
      auto wfHandle = evt.getHandle< vector< raw::OpDetWaveform > >(itag1);

      // Check that they are valid
      if (!wfHandle.isValid()) {
        mf::LogWarning("PDSNoiseFilter") << "No waveforms for label" << label;
        continue;
      }

      // Handle -> data
      vector< raw::OpDetWaveform > in_waveforms = *wfHandle;

      // Create a unique_ptr to store the output
      auto out_waveforms = std::make_unique< vector< raw::OpDetWaveform > >();

      // Loop through the waveforms applying filtering.
      for (auto const& in_wave : in_waveforms) {

        //raw::OpDetWaveform out_wave(in_wave.TimeStamp(), in_wave.ChannelNumber(), in_wave.size());
	
	Int_t nBins = in_wave.size();                     // size of the input vector
	std::vector<short unsigned int > out_wave(nBins); // vector in which the filtered waveform will be saved
        std::vector<short unsigned int > out_wave_double(nBins);
	  for(Int_t i=0; i<nBins; i++) out_wave_double[i]=2*in_wave[i];
	    // careful, it doesn't work for vector<short>, but it does for vector<short unsigned int>

        //####################################################################
        //                    Beginning of the Main Code
        //####################################################################

        //std:vector<short> out_wave(nBins); // don't worry about the shifted baseline due to the integer of the output
        float convert = 0;  // variable used to convert from float to short
        short convert2 = 0; // sum 0.5 to the final result to try to recover the lost information from the float-short conversion 
        std::vector<float> filtered = filter(out_wave_double,nBins); // filtered waveform
        
        // loop to convert the filtered output from float to short
        // and save it into the out_wave vector
        for(Int_t i=0; i<nBins; i++){
            convert = filtered[i];
            convert2 = (short) convert + 0.5;
            out_wave[i] = convert2;
        }

        //####################################################################
        //                      End of the Main Code
        //####################################################################

	raw::OpDetWaveform out_waveFinal(in_wave.TimeStamp(), in_wave.ChannelNumber(), out_wave);

        out_waveforms->emplace_back(std::move(out_waveFinal));

      }

      evt.put(std::move(out_waveforms),label);
    }

  }
  
  //####################################################################
  //                       Mobile Average Function
  //####################################################################

  std::vector<float> PDSNoiseFilter::filter( std::vector<short unsigned int> vec, Int_t nBins){
      
      Int_t mobileAVG = 5;                    // stablish how big the mobile average will be
      Double_t avg = 0;                       // saves the average value for each bin
      vector<float> firstFilter(nBins);       // values after mobile average
      vector<float> secondFilter(nBins);      // values after mobile average
      Int_t nComp = vec.size();               // size of the vector
      Int_t  c1=0;                            // counter for n>mobileAVG and n<nComp-mobileAVG
      Int_t  c2=0;                            // counter for n<=mobileAVG
      Int_t  c0=0;                            // counter for n>=nComp-mobileAVG
      
      for(Int_t n=0; n<nComp ; n++){
          if(n>mobileAVG && n<nComp-mobileAVG){
              for(Int_t i=n-mobileAVG; i<=n+mobileAVG; i++){avg = avg + vec[i]; c0=c0+1;}
              avg = avg/c0;
              firstFilter[n] = avg;
              avg=0;
              c0=0;
          }
          else{
              if(n<=mobileAVG){
                  for(Int_t i=0; i<=n+mobileAVG; i++){ avg = avg + vec[i]; c1=c1+1;}
                  avg = avg/c1;
                  firstFilter[n] = avg;
                  avg = 0;
                  c1=0;
              }
              else if(n>=nComp-mobileAVG){
                  for(Int_t i=n-mobileAVG; i<nComp; i++) {avg = avg + vec[i]; c2=c2+1;}
                  avg = avg/c2;
                  firstFilter[n] = avg;
                  avg = 0;
                  c2=0;
              }
          }
      }
      
      TV1D_denoise(firstFilter,secondFilter,nBins,10);
      
      
      
      return secondFilter;
  }
  
  //####################################################################
  //                          Denoise Function
  //####################################################################
  
  void PDSNoiseFilter::TV1D_denoise(vector<float> input, vector<float>& secondFilter, const int width, const float lambda) {
    
    if (width>0) { /*to avoid invalid memory access to input[0]*/
        int k=0, k0=0;  /*k: current sample location, k0: beginning of current segment*/
        float umin=lambda, umax=-lambda;/*u is the dual variable*/
        float vmin=input[0]-lambda, vmax=input[0]+lambda;/*bounds for the segment's value*/
        int kplus=0, kminus=0; /*last positions where umax=-lambda, umin=lambda, respectively*/
        const float twolambda=2.0*lambda;/*auxiliary variable*/
        const float minlambda=-lambda;/*auxiliary variable*/
        for (;;) {/*simple loop, the exit test is inside*/
            while (k==width-1) {/*we use the right boundary condition*/
                if (umin<0.0) { /*vmin is too high -> negative jump necessary*/
                    do secondFilter[k0++]=vmin; while (k0<=kminus);
                    umax=(vmin=input[kminus=k=k0])+(umin=lambda)-vmax;
                } else if (umax>0.0) {	/*vmax is too low -> positive jump necessary*/
                    do secondFilter[k0++]=vmax; while (k0<=kplus);
                    umin=(vmax=input[kplus=k=k0])+(umax=minlambda)-vmin;
                } else {
                    vmin+=umin/(k-k0+1); 
                    do secondFilter[k0++]=vmin; while(k0<=k);
                    return;
                }
            }
            if ((umin+=input[k+1]-vmin)<minlambda) {/*negative jump necessary*/
                do secondFilter[k0++]=vmin; while (k0<=kminus);
                vmax=(vmin=input[kplus=kminus=k=k0])+twolambda;
                umin=lambda; umax=minlambda;
            } else if ((umax+=input[k+1]-vmax)>lambda) {/*positive jump necessary*/
                do secondFilter[k0++]=vmax; while (k0<=kplus);
                vmin=(vmax=input[kplus=kminus=k=k0])-twolambda;
                umin=lambda; umax=minlambda;
            } else { 	/*no jump necessary, we continue*/
                k++;
                if (umin>=lambda) {/*update of vmin*/
                    vmin+=(umin-lambda)/((kminus=k)-k0+1);
                    umin=lambda;
                } 
                if (umax<=minlambda) {/*update of vmax*/
                    vmax+=(umax+lambda)/((kplus=k)-k0+1);
                    umax=minlambda;
                }
            }

        }

    }

  }
  DEFINE_ART_MODULE(PDSNoiseFilter)
}
