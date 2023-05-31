// ===========================================================================
// Deconvolution_module.cc
// This module produces deconvoluted signals output
// Using the signal digitized with the template in digitizer module as input.
// Creating OpWaveform object
// Filter wiener/gauss -FFT
// @authors     : Daniele Guffanti, Maritza Delgado, Sergio Manthey Corchado
// @created     : Jan 26, 2022 
//============================================================================= 
 
#ifndef Deconvolution_h
#define Deconvolution_h
 
// Framework includes

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ART extensions
#include "nurandom/RandomUtils/NuRandomService.h"

// LArSoft includes
#include "larcore/CoreUtils/ServiceUtil.h" 
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpWaveform.h"
#include "lardata/Utilities/LArFFT.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/Simulation/BeamGateInfo.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"

 // CLHEP includes
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandExponential.h"
#include "CLHEP/Random/RandFlat.h"

// C++ Includes
#include <map>
#include <memory>
#include <string>
#include <stdio.h>
#include <iostream>
#include <vector>

// ROOT includes
#include "TFile.h"
#include "TMath.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "RooGaussian.h"
#include "TVirtualFFT.h"
#include "TStyle.h"
#include "TLine.h"
#include "TF1.h"
#include "TComplex.h"
#include "TFFTComplexReal.h"
#include <Rtypes.h> 

using std::string;
using std::vector;


namespace opdet {
  class Deconvolution : public art::EDProducer{
    public:
      struct Config {
        public:
          fhicl::Atom<std::string> module_type{ fhicl::Name("module_type") }; 
          fhicl::Atom<std::string> InputModule{ fhicl::Name("InputModule") }; 
          fhicl::Atom<std::string> InstanceName{ fhicl::Name("InstanceName") }; 
          fhicl::Atom<Double_t>    LineNoiseRMS{ fhicl::Name("LineNoiseRMS"), 1.0 };
          fhicl::Atom<size_t>      PreTrigger{ fhicl::Name("PreTrigger"), 0};
          fhicl::Atom<short>       Pedestal{ fhicl::Name("Pedestal"), 1500}; 
          fhicl::Atom<std::string> DigiDataFile{ fhicl::Name("DigiDataFile") }; 
          fhicl::Atom<size_t>      DigiDataColumn{ fhicl::Name("DigiDataColumn"), 1 }; 
          fhicl::Atom<Int_t>       Samples{ fhicl::Name("Samples"), 1000 };
          fhicl::Atom<Int_t>       PedestalBuffer{ fhicl::Name("PedestalBuffer"), 10 }; 
          fhicl::Atom<Double_t>    Scale{ fhicl::Name("Scale"), 1 }; 
          fhicl::Atom<bool>        ApplyPostfilter{ fhicl::Name("ApplyPostfilter"), false };
          fhicl::Atom<bool>        ApplyPostBLCorrection{ fhicl::Name("ApplyPostBLCorrection") }; 
          fhicl::Atom<bool>        AutoScale{ fhicl::Name("AutoScale"), false };
          fhicl::Atom<std::string> OutputProduct{ fhicl::Name("OutputProduct"), "decowave"};

          struct Filter {
            fhicl::Atom<std::string> Name{fhicl::Name("Name"), "Gauss"};
            fhicl::Atom<Double_t>    Cutoff{fhicl::Name("Cutoff"), 2};  
          };
          fhicl::Table<Config::Filter> Filter{ fhicl::Name("WfmFilter") };

          struct ExtraFilter {
            fhicl::Atom<std::string> Name{fhicl::Name("Name"), "Gauss"}; 
            fhicl::Atom<Double_t>    Cutoff{fhicl::Name("Cutoff"), 2}; 
          }; 
          fhicl::Table<Config::ExtraFilter> Postfilter{ fhicl::Name("WfmPostfilter") };
      };

      using Parameters = art::EDProducer::Table<Config>;

      //! Input signal shape model,for Wiener filter could use delta or scintillation light signal.
      enum EInputShape {kDelta = 0, kScint = 1};
      //! Waveform filter type
      enum EFilterType {kOther = 0, kWiener = 1, kGauss = 2, kNone = 3}; 

      struct WfmExtraFilter_t {
        TString fName; 
        Double_t fCutoff; 

        WfmExtraFilter_t() {fName = ""; fCutoff = 1.0;}

        WfmExtraFilter_t(const char* name, Double_t cutoff) {
          fName = name; fCutoff = cutoff;
        }

        WfmExtraFilter_t(const struct Config::ExtraFilter& config) {
          fName = config.Name();
          fCutoff = config.Cutoff(); 
        }
      }; 

      struct WfmFilter_t {
        TString fName; 
        EFilterType fType; 
        Double_t fCutoff; 

        WfmFilter_t() : fName("Unknown_filter"), fType(kOther), fCutoff(32) {}
        WfmFilter_t(const char* name) : fName(name), fType(kOther), fCutoff(32) {}
        WfmFilter_t(const struct Config::Filter& config) {
          fName = config.Name(); 
          fCutoff = config.Cutoff(); 
          if (fName == "Wiener") fType = kWiener; 
          else if (fName == "Gauss") fType = kGauss; 
          else fType = kOther; 
        }
      };

      /**
       * @brief Helper struct to handle waveforms' FT
       *
       * Helper data struct to store and handle waveforms' Fourier Transform. 
       */
      struct CmplxWaveform_t {
        //! vector of Fourier coefficients
        std::vector<TComplex> fCmplx; 
        //! vector of the real part of Fourier coefficients
        std::vector<Double_t> fRe; 
        //! vector of the imaginary part of Fourier coefficients
        std::vector<Double_t> fIm; 

        //! Empty constructor
        inline CmplxWaveform_t() {}
        //! Constructor
        inline CmplxWaveform_t(int size) {
          fCmplx = std::vector(size, TComplex(0., 0)); 
          fRe = std::vector(size, 0.); 
          fIm = std::vector(size, 0.); 
        }

        //! Copy constructor
        inline CmplxWaveform_t(const CmplxWaveform_t& cwf) {
          fCmplx = std::vector( cwf.fCmplx ); 
          fRe = std::vector( cwf.fRe ); 
          fIm = std::vector( cwf.fIm ); 
        }

        //! Destructor
        inline ~CmplxWaveform_t() {
          fCmplx.clear(); fRe.clear(); fIm.clear(); 
        }

        /**
         * @brief Point-side product of complex waveforms
         *
         * Compute the product of two complex waveforms point 
         * by point
         */
        inline CmplxWaveform_t operator*(const CmplxWaveform_t& cwf) const {
          CmplxWaveform_t result(fCmplx.size()); 
          if (cwf.fCmplx.size() != fCmplx.size()) {
            printf("opdet::Deconvolution_module::CmplxWaveform_t::operator* ERROR"); 
            printf(" waveforms size does not match.\n"); 
            return result;
          }

          for (size_t i=0; i<fCmplx.size(); i++) {
            result.fCmplx.at(i) = fCmplx.at(i) * cwf.fCmplx.at(i); 
            result.fRe.at(i) = result.fCmplx.at(i).Re(); 
            result.fIm.at(i) = result.fCmplx.at(i).Im(); 
          }

          return result;
        }

        //! Set the complex coefficient `i` given its real and imaginary part
        inline void MakeCmplx(size_t const i) { 
          fCmplx.at(i) = TComplex(fRe.at(i), fIm.at(i)); 
        }

        //! Set the complex coefficients given real and imaginary parts
        //!
        //! Note that the second half of the waveform is set to zero
        inline void MakeCmplx() {
          for (size_t i=0; i<fCmplx.size(); i++) {
            MakeCmplx(i);
          }
          return;
        }

        //! Set the real and imaginary part of the coefficient `i` 
        inline void MakeReAndIm(size_t const i) {
          fRe.at(i) = fCmplx.at(i).Re(); 
          fIm.at(i) = fCmplx.at(i).Im(); 
        }

        //! Set real and imaginary parts from the complex coefficients
        //!
        //! Note that the second half of the waveform is set to zero
        inline void MakeReAndIm() {
          for (size_t i=0; i<fCmplx.size(); i++) {
            MakeReAndIm(i); 
          } 
          return;
        }
      };


      explicit Deconvolution(Parameters const&);
      virtual ~Deconvolution();
      void produce(art::Event& evt);

      // Parameters we'll read from the fcl-file
      std::string fInputModule;                 //!< Module used to create OpDetWaveforms
      std::string fInstanceName;                //!< Input tag for OpDetWaveforms collection
      double fSampleFreq;                       //!< Sampling frequency in MHz 
      short  fPedestal;                         //!< In ADC counts
      double  fLineNoiseRMS;                    //!< Pedestal RMS in ADC counts
      size_t fPreTrigger;                       //!< In ticks
      std::string fDigiDataFile;                //!< single p.e. template source file
      size_t fDigiDataColumn;                   //!< single p.e. template source file column    
      double fScale;                            //!< Scaling of resulting wvfs 
      int fSamples;                             //!< Same value as ReadoutWindow in digitizer
      int fPedestalBuffer;                      //!< Used to calculate pedestal which is definded as PreTrigger - PedestalBuffer in [ticks]
                                                //!< default is 10 ticks -> should be adapted to resulting peak width (after deconvolution).
      bool fApplyPostfilter; 
      bool fApplyPostBLCorr;
      bool fAutoScale;
      
      // Additional parameters here
      
      //--- Single photoelectron variables
      std::vector<double> fSinglePEWaveform;    //!< Vector that stores the single PE template in ADC*us
      double fSinglePEAmplitude;                //!< single PE amplitude for found maximum peak in the template.
      unsigned int WfDeco;                      //!< Number of waveform processed
      
      //--------Filter Variables
      std::string fOutputProduct; 
      WfmExtraFilter_t fPostfilterConfig;
      WfmFilter_t fFilterConfig; 
      
      //Input signal shape model for Wiener filter could use delta "kDelta" or scintillation light signals "kScint"
      EInputShape fInputShape = kDelta; 
      // EInputShape fInputShape = kScint; //uncomment if set to kScint

    private:
      int  CountFileColumns(const char* file_path);
      void SourceSPEDigiDataFile(); 
      void BuildExtraFilter(CmplxWaveform_t& xF0, const WfmExtraFilter_t config); 
      void ComputeExpectedInput(std::vector<double>& s, double nmax);
      void CopyToOutput(const std::vector<float>& v, std::vector<float>& target); 
      void CopyToOutput(const CmplxWaveform_t& v, std::vector<float>& target); 
      Double_t ComputeAutoNormalization(CmplxWaveform_t& xGH, const float thrs=0); 
      TVirtualFFT* fft_r2c; 
      TVirtualFFT* fft_c2r; 
      CmplxWaveform_t fxG0; 
      CmplxWaveform_t fxG1; 
  };
}
#endif

namespace opdet{
  DEFINE_ART_MODULE(Deconvolution)
}

namespace opdet {
  //---------------------------------------------------------------------------
  // Constructor
  Deconvolution::Deconvolution(const Parameters& pars)
    : EDProducer{pars}, 
    fInputModule{ pars().InputModule()}, 
    fInstanceName{ pars().InstanceName()},
    
    fPedestal{ pars().Pedestal()},
    fLineNoiseRMS{ pars().LineNoiseRMS() },
    fPreTrigger{ pars().PreTrigger()},

    fDigiDataFile{ pars().DigiDataFile()},
    fDigiDataColumn{ pars().DigiDataColumn()},
    fScale{ pars().Scale()},
    fSamples{ pars().Samples()},
    fPedestalBuffer{ pars().PedestalBuffer()},
    fApplyPostfilter{ pars().ApplyPostfilter()},
    fApplyPostBLCorr{ pars().ApplyPostBLCorrection()},
    fAutoScale{ pars().AutoScale()},
    fOutputProduct{ pars().OutputProduct() },
    fPostfilterConfig{ WfmExtraFilter_t( pars().Postfilter()) }, 
    fFilterConfig{ WfmFilter_t( pars().Filter() ) }, 
    fxG0(fSamples), 
    fxG1(fSamples)
  {  
    // Declare that we'll produce a vector of OpDetWaveforms 
    WfDeco=0;  

    // This module produces
    produces< std::vector< recob::OpWaveform> > ();   


    // Obtain parameters from TimeService
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
    fSampleFreq = clockData.OpticalClock().Frequency();

    fft_r2c = TVirtualFFT::FFT(1, &fSamples, "M R2C K");
    fft_c2r = TVirtualFFT::FFT(1, &fSamples, "M C2R K");
    
    SourceSPEDigiDataFile(); 

    // build post filter (if required)
    if (fApplyPostfilter) BuildExtraFilter(fxG1, fPostfilterConfig); 

  }

  //---------------------------------------------------------------------------
  // Destructor
  Deconvolution::~Deconvolution(){
    delete fft_r2c; 
    delete fft_c2r;
    
  }    

  //-------------------------------------------------------------------------
  void Deconvolution::produce(art::Event& evt){

    art::Handle< std::vector< raw::OpDetWaveform > > wfHandle;
    evt.getByLabel(fInputModule, fInstanceName, wfHandle);
   
    //******************************
    //-- Read Waveform----
    //****************************** 

    std::vector<raw::OpDetWaveform> digi_wave = *wfHandle; 
    int NOpDetWaveform = digi_wave.size(); //Number of waveforms in OpDetWaveform Object
    std::cout << NOpDetWaveform << std::endl;

    //Pointer that will store produced DecoWaveform
    auto out_recob = std::make_unique< std::vector< recob::OpWaveform > >();

    std::vector<short unsigned int > out_digiwave(fSamples); //vector in which the waveform will be saved
    std::vector<float> out_recob_float(fSamples);        //vector in which the decowaveform will be saved, using float
    std::vector<double> xv(fSamples, 0.); 
    
    //******************************
    //--- Setup filter's components 
    //******************************

    // xH: FFT of spe response
    std::vector<double> xh(fSinglePEWaveform);
    CmplxWaveform_t xH(fSamples);  
    fft_r2c->SetPoints(&xh[0]);  
    fft_r2c->Transform();
    fft_r2c->GetPointsComplex(&xH.fRe[0], &xH.fIm[0]);
    xH.MakeCmplx(); 
    
    //******************************
    //--- Process waveforms
    //******************************
    for (auto const& wf: digi_wave) {
      CmplxWaveform_t xV(fSamples); 
      CmplxWaveform_t xS(fSamples); 
      CmplxWaveform_t xN(fSamples); 
      CmplxWaveform_t xG(fSamples); 
      CmplxWaveform_t xY(fSamples); 
      CmplxWaveform_t xGH(fSamples); 
      std::vector<float> xSNR(fSamples, 0.); 
      int OriginalWaveformSize = wf.Waveform().size();
      //----------------------Resize deconvoluted signals (using floats) to original waveform size
      if (static_cast<int>(OriginalWaveformSize) <= fSamples) { 
        out_recob_float.resize(fSamples,0); 
      }

      else {
         out_recob_float.resize(fSamples);
      }

      for (Int_t i= 0; i < fSamples; i++){
          // Remove baseline 
          if (i < static_cast<int>(wf.Waveform().size())) xv[i] = (wf[i]-fPedestal);
          // if waveform is shorter than fSamples fill the rest with noise
          else xv[i] = CLHEP::RandGauss::shoot(0, fLineNoiseRMS);
       }
       
      //---------------------------------------------------- Guess input signal
      // Found maximum peak in the Waveform
      Double_t SPE_Max = 0;
      double maxADC=*max_element(xv.begin(),xv.end());
      double maxAmplit= maxADC; // Pedestal already subtracted
      SPE_Max = maxAmplit/fSinglePEAmplitude;
      
      std::vector<double>xs(fSamples,0.);
      // Compute expected input (using a delta or the scint tile profile as a model)
         ComputeExpectedInput(xs, SPE_Max); 
         
      //-------------------------------------------------- Compute waveform FFT
      fft_r2c->SetPoints(&xv[0]);
      fft_r2c->Transform();
      fft_r2c->GetPointsComplex(&xV.fRe[0], &xV.fIm[0]);
      xV.MakeCmplx(); 

      //----------------------------------------------------- Compute input FFT
      fft_r2c->SetPoints(&xs[0]);
      fft_r2c->Transform();
      fft_r2c->GetPointsComplex(&xS.fRe[0], &xS.fIm[0]);
      xS.MakeCmplx(); 

      //******************************
      // Compute filters.
      //******************************      

      Double_t fFrequencyCutOff = fFilterConfig.fCutoff; 
      Double_t fTickCutOff = fSamples*fFrequencyCutOff/fSampleFreq;                       
            
      for (int i=0; i<fSamples*0.5+1; i++) {
        // Compute spectral density
        double H2 = xH.fCmplx.at(i).Rho2();
        double S2 = xS.fCmplx.at(i).Rho2();
        double N2 = fLineNoiseRMS * fLineNoiseRMS * fSamples ;
        
        xN.fCmplx.at(i) = TComplex(sqrt(N2), 0.); 

        if (fFilterConfig.fType == Deconvolution::kWiener){
          // Compute Wiener filter
          xSNR.at(i) = S2 / N2; 
          xG.fCmplx.at(i) = TComplex::Conjugate(xH.fCmplx.at(i))*S2 / (H2*S2 + N2);
        }

        else if (fFilterConfig.fType == Deconvolution::kGauss){
          // Compute gauss filter
          xG.fCmplx[0] = TComplex(0,0);
          xG.fCmplx.at(i) = TComplex::Exp(
            -0.5*TMath::Power(i*fSampleFreq/(fSamples*fTickCutOff),2))
            /xH.fCmplx.at(i);
        }

        else if (fFilterConfig.fType == Deconvolution::kOther){
          // Compute dec signal
          xG.fCmplx.at(i) = TComplex(1,0)/xH.fCmplx.at(i); // Standard deconv is just the division of signal and SPE template in Fourier space
        }

         xG.MakeReAndIm(i); 
      }

      // Apply filter to the waveform
      xY  = xG * xV;
      xY.MakeReAndIm();

      //Transform of the filtered signal
      fft_c2r->SetPointsComplex(&xY.fRe[0], &xY.fIm[0]);
      fft_c2r->Transform();
      double *xy = fft_c2r->GetPointsReal();
      std::vector<double> xvdec(xy, xy+fSamples); 
     
      Double_t scale = 1.0 / fSamples;
      if (!fAutoScale) {scale = fScale;} 
      else {
        // Apply filter to the detector response (for normalization)
        xGH = xG * xH; 
        xGH.MakeReAndIm(); 
        Double_t filter_norm = ComputeAutoNormalization(xGH); 
        scale = filter_norm / (Double_t)fSamples;
      }

      if (fApplyPostfilter) {
        CmplxWaveform_t xxY(fSamples); 
        std::vector<double> ytmp(xvdec.begin(), xvdec.end()); 
        fft_r2c->SetPoints(&ytmp[0]); 
        fft_r2c->Transform(); 
        fft_r2c->GetPointsComplex(&xxY.fRe[0], &xxY.fIm[0]); 
        xxY.MakeCmplx(); 
        xxY = xxY * fxG1; 
        xxY.MakeReAndIm(); 
        fft_c2r->SetPointsComplex(&xxY.fRe[0], &xxY.fIm[0]); 
        fft_c2r->Transform(); 
        double* xxy = fft_c2r->GetPointsReal(); 
        double g1_scale = 1.0 / fSamples; 
        for (int i=0; i<fSamples; i++) {
          xvdec[i] = (xxy[i] * g1_scale); 
        }
      }
      //
       // Correct baseline after deconvolution
      double decPedestal = 0;
      if (fApplyPostBLCorr) {
        for (size_t i=0; i<fPreTrigger-fPedestalBuffer; i++){
          decPedestal = decPedestal + xvdec[i];
        }
        decPedestal = decPedestal/int(fPreTrigger-fPedestalBuffer);
        
      }

      for (int i=0; i<fSamples; i++){ 
        out_recob_float[i] = (xvdec[i]-decPedestal)*scale;
      }
      
      //---------------------Resize deconvoluted signals (using floats) to original waveform size
      if (int(out_recob_float.size()) <= OriginalWaveformSize) { 
         out_recob_float.resize(OriginalWaveformSize,0); 
      }

      else { 
        out_recob_float.resize(OriginalWaveformSize);
      }
      
   
      if        (strcmp(fOutputProduct.c_str(), "H" ) == 0) {
        CopyToOutput(xH, out_recob_float); 
      } else if (strcmp(fOutputProduct.c_str(), "S" ) == 0) {
        CopyToOutput(xS, out_recob_float); 
      } else if (strcmp(fOutputProduct.c_str(), "N" ) == 0) {
        CopyToOutput(xN, out_recob_float); 
      } else if (strcmp(fOutputProduct.c_str(), "G0") == 0) {
        CopyToOutput(fxG0, out_recob_float); 
      } else if (strcmp(fOutputProduct.c_str(), "G1") == 0) {
        CopyToOutput(fxG1, out_recob_float); 
      } else if (strcmp(fOutputProduct.c_str(), "V" ) == 0) {
        CopyToOutput(xV, out_recob_float); 
      } else if (strcmp(fOutputProduct.c_str(), "v" ) == 0) {
        fft_c2r->SetPointsComplex(&xV.fRe[0], &xV.fIm[0]); 
        fft_c2r->Transform(); 
        Double_t* vv = fft_c2r->GetPointsReal(); 
        for (int i=0; i<fSamples; i++) out_recob_float.at(i) = scale* vv[i]; 
      } else if (strcmp(fOutputProduct.c_str(), "SNR")== 0) {
        CopyToOutput(xSNR, out_recob_float); 
      } else if (strcmp(fOutputProduct.c_str(), "G" ) == 0) {
        CopyToOutput(xG, out_recob_float);
      }

      recob::OpWaveform decwav(wf.TimeStamp(), wf.ChannelNumber(), out_recob_float );
      out_recob->emplace_back(std::move(decwav));   
    }//waveforms loop
         
    //-------------------------------------Print fType Filter
       if (fFilterConfig.fType == Deconvolution::kWiener){printf("***Wiener Filter****");}
       if (fFilterConfig.fType == Deconvolution::kGauss){printf("***Gauss Filter***");}
       if (fFilterConfig.fType == Deconvolution::kOther){printf("***Standart dec***");}
       if (fApplyPostfilter){printf("***ApplyPostfilter***");}
       
      //------------------------------------------------
         
    // Push the OpDetWaveforms and OpWaveform into the event
    evt.put(std::move(out_recob));
    WfDeco++;
  }


  /**
   * @brief Build a filter to be applied after the deconvolution
   *
   * Construct an extra filter to be applied after 
   * the deconvolution process.
   * Different filters can be implemented by switching the flag 
   * `fPostfilterConfig.fName`
   * via the Config::ExtraFilter::name parameter. 
   *
   * @param xF
   */
  void Deconvolution::BuildExtraFilter(CmplxWaveform_t& xF, const WfmExtraFilter_t config) {
    if (config.fName != "Gauss") {
      printf("Deconvolution::BuildExtraFilter WARNING: Unknown filter model %s. Skip.\n", 
          config.fName.Data()); 
      return;
    } 

    // Compute sigma corresponding to the given cutoff frequency
    const Double_t df       = fSampleFreq / (Double_t)fSamples; 
    const Double_t cutoff   = config.fCutoff / df; 
    const Double_t k_cutoff = sqrt(log(2)); 
    const double sigma = fSamples * k_cutoff / (TMath::TwoPi() * cutoff);
    const int    mu    = 0.5*fSamples; 

    printf("Deconvolution::BuildExtraFilter sigma is %g\n", sigma); 

    std::vector<Double_t> xf(fSamples, 0.);
    for (int i=0; i<fSamples; i++) {
      xf.at(i) = TMath::Gaus(i, mu, sigma, kTRUE); 
    }

    std::vector<Double_t> re_(fSamples, 0.); 
    std::vector<Double_t> im_(fSamples, 0.);

    fft_r2c->SetPoints(&xf[0]); 
    fft_r2c->Transform(); 
    fft_r2c->GetPointsComplex(&re_[0], &im_[0]);

    for (int i=0; i<0.5*fSamples+1; i++) {
      TComplex F(re_.at(i), im_.at(i));  
      TComplex phase = TComplex(0., -TMath::TwoPi()*i*mu/(fSamples));  
      xF.fCmplx.at(i) = F*TComplex::Exp(phase); 
      xF.MakeReAndIm(i); 
    }

    return; 
  }


  /**
   * @brief Compute expected input signal
   *
   * Produce a waveform representing a guess of the input signal based on the
   * estimated nr of p.e. in the waveform peak. Based on the `fInputShape`
   * variable, the input shape can be assumed as a δ-function scaled for 
   * the nr of p.e. or as the LAr scintillation time profile
   *
   * @param s input signal
   * @param nmax estimated nr of p.e. at the waveform max
   */
  void Deconvolution::ComputeExpectedInput(std::vector<double>& s, double nmax) {
    if (fInputShape == kScint) {
      //ST profile 
      auto const *Larprop = lar::providerFrom<detinfo::LArPropertiesService>();
      std::vector<double>SignalTime={(Larprop->ScintFastTimeConst()* 0.001),(Larprop->ScintSlowTimeConst()* 0.001)};
      std::vector<double>SignalScint={Larprop->ScintYieldRatio(),1.-Larprop->ScintYieldRatio()};

      double dt = 1/fSampleFreq; 
      double t = 0;

      for (size_t i=0; i<s.size(); i++) {
        double lightsignal=0;
        for (size_t j=0; j<SignalTime.size();j++){
          lightsignal+=SignalScint[j]*exp(-t/SignalTime[j]);
        }
        // s.at(i) = nmax*lightsignal; 
        s.at(i) = lightsignal; 
        t+=dt; 
      }
    }
    else if (fInputShape == kDelta) {
      s.at(fSamples*0.5) = nmax; 
    }
    else {
      printf("Deconvolution::ComputeExpectedInput WARNING\n"); 
      printf("Unknown input shape: assuming Deconvolution::kDelta\n"); 
      s.at(1) = nmax;
    }
    return; 
  }

  /**
   * @brief Source the single p.e. response from file
   *
   * Source the single p.e. response template from the txt file set by 
   * `fDigiDataFile` and set the variable `fSinglePEAmplitude` with the 
   * amplitude of the single p.e. response. In case of a multi-column
   * template file, the relevant column can be selected by setting the
   * varible `fDigiDataColumn`.
   */
  void Deconvolution::SourceSPEDigiDataFile() {
    size_t n_columns = CountFileColumns(fDigiDataFile.c_str()); 
    if (fDigiDataColumn >= n_columns) {
      printf("Deconvolution::SourceSPETemplate ERROR: "); 
      printf("The module is supposed to select column %lu, but only %lu columns are present.\n", 
          fDigiDataColumn, n_columns); 
      throw art::Exception(art::errors::InvalidNumber); 
    }

    std::ifstream SPEData; 
    SPEData.open(fDigiDataFile); 
    Double_t buff[100] = {0};  

    std::string temp_str; 
    if (SPEData.is_open()) {
      while (std::getline(SPEData, temp_str)) {
        std::stringstream ss; ss << temp_str; 
        int  icol = 0; 
        while (ss) {ss >> buff[icol]; ++icol;}

        fSinglePEWaveform.push_back(buff[fDigiDataColumn]);
      } 
    } else {
      printf("Deconvolution::produce ERROR "); 
      printf("Cannot open SPE template file.\n"); 

      throw art::Exception(art::errors::FileOpenError); 
    }

    fSinglePEWaveform.resize(fSamples, 0);

    SPEData.close(); 

    // Set single p.e. maximum value
    fSinglePEAmplitude = TMath::Max(1.0, 
      *(std::max_element(fSinglePEWaveform.begin(), fSinglePEWaveform.end()))); 
    std::cout << "SPE Amplitude: " << fSinglePEAmplitude << std::endl;
    return;
  }

  /**
   * @brief Count the nr of column in a txt file
   *
   * @param file_path
   *
   * @return nr of columns
   */
  int Deconvolution::CountFileColumns(const char* file_path) {
    std::ifstream file_; 
    file_.open(file_path); 

    if (!file_.is_open()) {
      printf("Deconvolution::CountFileColumns(%s) ERROR:\n", 
          file_path); 
      printf("Unable to open file."); 
      throw art::Exception(art::errors::FileOpenError); 
    }

    int N_COLUMNS = 0; 
    std::string line; 
    int iline = 0; 
    while ( std::getline(file_, line) ) {
      std::stringstream sstream; 
      sstream << line;
      std::string sub;
      int n_columns = 0; 

      while (sstream) {
        sstream >> sub; 
        if (sub.length()) ++n_columns; 
      }

      if (iline == 1) {N_COLUMNS = n_columns;}
      else if (iline > 1) {
        if (n_columns != N_COLUMNS) {
          printf("Deconvolution::CountFileColumns(%s): WARNING ", 
              file_path); 
          printf("Nr of columns change along the file!\n"); 
          N_COLUMNS = n_columns; 
        }
      }
      iline++; 
    }
    file_.close();  
    return N_COLUMNS; 
  }

  /**
   * @brief Compute normalization factor for a given filter
   *
   * The filter normalization factor is obtained by applying the 
   * filter to the single p.e. response template (i.e., to a 1 p.e. noiseless signal).
   * the product of the filter should be a the best approximation of the 
   * input function achievable with the signal to noise ratio in such waveform, 
   * thus, the normalization consists in the integral of the product 
   * around the peak in a region defined as the one where the signal is positive. 
   * This factor is supposed to tend to 1 for high SNR (>= 10).
   *
   * @param xGH
   *
   * @return filter normalization
   */
  Double_t Deconvolution::ComputeAutoNormalization(CmplxWaveform_t& xGH, const float thrs) {


    // Apply a linear phase shift to make the "pulse" in the middle of the
    // time window
    const int shift = 0.5*fSamples;
    for (int i=0; i<fSamples*0.5+1; i++) {
     TComplex phase = TComplex(0., TMath::TwoPi()*i*shift/(fSamples)); 
     xGH.fCmplx.at(i) = xGH.fCmplx.at(i)*TComplex::Exp(phase);
    }
    xGH.MakeReAndIm(); 

    fft_c2r->SetPointsComplex(&xGH.fRe[0], &xGH.fIm[0]); 
    fft_c2r->Transform(); 
    Double_t* x_ =  fft_c2r->GetPointsReal(); 
    std::vector<Double_t> x(x_, x_ + fSamples);

    Int_t imax = std::distance(x.begin(), std::max_element(x.begin(), x.end())); 
    Int_t ileft = imax; 
    Int_t iright = imax; 

    while (x[ileft]  > thrs && ileft > 0) ileft--; 
    while (x[iright] > thrs && iright < fSamples) iright++; 

    double norm = 0; 
    for (Int_t k=ileft; k<=iright; k++) norm += x[k]; 
    norm /= (Double_t)fSamples; 

    if (norm > 1.0) norm = 1.0; 
    else if (norm <= 0.){
      printf("Deconvolution::ComputeAutoNormalization() WARNING: "); 
      printf(" bad normalization (%g), force to 1.0\n", norm);  
      norm = 1.0; 
    }
    
    return 1.0 / norm; 
  }

  void Deconvolution::CopyToOutput(const std::vector<float>& v, std::vector<float>& target) {
    target.assign(v.begin(), v.end());
    return;
  }
  void Deconvolution::CopyToOutput(const CmplxWaveform_t& v, std::vector<float>& target) {
    for (size_t i=0; i<v.fCmplx.size(); i++) {
      target.at(i) = v.fCmplx.at(i).Rho2(); 
    }
    return; 
  }
} // namespace opdet
