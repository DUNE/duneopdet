// ===============================================================================
// OpHitFinderDeco_module.cc
// This module is based on the larana/OpHitFinder module. It has been updated
// to deal with deconvolved waveforms.The module takes either raw::OpDetWaveforms
// (raw signals) or recob:OpWaveforms (deconvolved signals) as input,  
// and generates OpHits as output.
// The HitFinder algorithm can be chosen by the user.
// Added the scaling factor inside the new RunHitFinder_deco function.  
// It scales the values of the deconvolved signals before the hit finder.
//
// @authors     : Daniele Guffanti, Maritza Delgado, Sergio Manthey Corchado
// @created     : Oct, 2022 
//================================================================================
 
 // LArSoft includes
#include "larcore/CoreUtils/ServiceUtil.h" 
#include "larcore/Geometry/Geometry.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoCFD.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoFixedWindow.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoSiPM.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoSlidingWindow.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoThreshold.h"
#include "OpHitAlg_deco.h"
#include "larana/OpticalDetector/OpHitFinder/PMTPulseRecoBase.h"
#include "larana/OpticalDetector/OpHitFinder/PedAlgoEdges.h"
#include "larana/OpticalDetector/OpHitFinder/PedAlgoRollingMean.h"
#include "larana/OpticalDetector/OpHitFinder/PedAlgoUB.h"
#include "larana/OpticalDetector/OpHitFinder/PulseRecoManager.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpWaveform.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/Simulation/BeamGateInfo.h"
#include "larreco/Calibrator/IPhotonCalibrator.h"
#include "larreco/Calibrator/IPhotonCalibratorService.h"
#include "larreco/Calibrator/PhotonCalibratorStandard.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Utilities/Exception.h"
#include "fhiclcpp/ParameterSet.h"

// ROOT includes

// C++ Includes
#include <map>
#include <memory>
#include <string>

namespace {
  template <typename T>
  pmtana::PMTPulseRecoBase* thresholdAlgorithm(fhicl::ParameterSet const& hit_alg_pset,std::optional<fhicl::ParameterSet> const& rise_alg_pset){
    if (rise_alg_pset){return new T(hit_alg_pset, art::make_tool<pmtana::RiseTimeCalculatorBase>(*rise_alg_pset) );}
    else{return new T(hit_alg_pset, nullptr);}
  }
}

namespace opdet {
  class OpHitFinderDeco : public art::EDProducer {
    public:
      // Standard constructor and destructor for an ART module.
      explicit OpHitFinderDeco(const fhicl::ParameterSet&);
      virtual ~OpHitFinderDeco();

      // The producer routine, called once per event.
      void produce(art::Event&);

    private:
      std::map<int, int> GetChannelMap();
      std::vector<double> GetSPEScales();
      std::vector<double> GetSPEShifts();

      // The parameters we'll read from the .fcl file.
      std::string fInputModule; // Input tag for OpDetWaveform collection
      std::string fInputModuledigi; 
      std::string fGenModule;
      std::string fInputDigiType;
      std::vector<std::string> fInputLabels;
      std::set<unsigned int> fChannelMasks;

      pmtana::PulseRecoManager fPulseRecoMgr;
      pmtana::PMTPulseRecoBase* fThreshAlg;
      pmtana::PMTPedestalBase* fPedAlg;

      Float_t fHitThreshold;
      Float_t fScale;
      unsigned int fMaxOpChannel;
      bool fUseStartTime;

      calib::IPhotonCalibrator const* fCalib = nullptr;
  };

}

namespace opdet {
  DEFINE_ART_MODULE(OpHitFinderDeco)
}

namespace opdet {
  //----------------------------------------------------------------------------
  // Constructor
  OpHitFinderDeco::OpHitFinderDeco(const fhicl::ParameterSet& pset) : EDProducer{pset}, fPulseRecoMgr(){
    // Indicate that the Input Module comes from .fcl
    fInputModule        = pset.get<std::string>("InputModule");
    fInputModuledigi    = pset.get<std::string>("InputModuledigi");
    fGenModule          = pset.get<std::string>("GenModule");
    fInputDigiType      = pset.get<std::string>("InputDigiType");
    fInputLabels        = pset.get<std::vector<std::string>>("InputLabels");
    fUseStartTime       = pset.get<bool>("UseStartTime", false);
    fHitThreshold       = pset.get<float>("HitThreshold");
    fScale              = pset.get<float>("ScalingFactor");
    bool useCalibrator  = pset.get<bool>("UseCalibrator");

    for (auto const& ch : pset.get<std::vector<unsigned int>>("ChannelMasks", std::vector<unsigned int>()))fChannelMasks.insert(ch);

    auto const& geometry(*lar::providerFrom<geo::Geometry>());
    fMaxOpChannel = geometry.MaxOpChannel();
    
    // If useCalibrator, get it from ART
    if (useCalibrator) {fCalib = lar::providerFrom<calib::IPhotonCalibratorService>();}
    // If not useCalibrator, make an internal one based on fhicl settings to hit finder.
    else {
      bool areaToPE   = pset.get<bool>("AreaToPE");
      float SPEArea   = pset.get<float>("SPEArea");
      float SPEShift  = pset.get<float>("SPEShift", 0.);
      // Reproduce behavior from GetSPEScales()
      if (!areaToPE) SPEArea = 20;
      // Delete and replace if we are reconfiguring
      if (fCalib) { delete fCalib; }
      fCalib = new calib::PhotonCalibratorStandard(SPEArea, SPEShift, areaToPE);
    }

    // Initialize the rise time calculator tool
    auto const rise_alg_pset  = pset.get_if_present<fhicl::ParameterSet>("RiseTimeCalculator");

    // Initialize the hit finder algorithm
    auto const hit_alg_pset   = pset.get<fhicl::ParameterSet>("HitAlgoPset");
    std::string threshAlgName = hit_alg_pset.get<std::string>("Name");
    if (threshAlgName == "Threshold")
      fThreshAlg = thresholdAlgorithm<pmtana::AlgoThreshold>(hit_alg_pset, rise_alg_pset);
    else if (threshAlgName == "SiPM")
      fThreshAlg = thresholdAlgorithm<pmtana::AlgoSiPM>(hit_alg_pset, rise_alg_pset);
    else if (threshAlgName == "SlidingWindow")
      fThreshAlg = thresholdAlgorithm<pmtana::AlgoSlidingWindow>(hit_alg_pset, rise_alg_pset);
    else if (threshAlgName == "FixedWindow")
      fThreshAlg = thresholdAlgorithm<pmtana::AlgoFixedWindow>(hit_alg_pset, rise_alg_pset);
    else if (threshAlgName == "CFD")
      fThreshAlg = thresholdAlgorithm<pmtana::AlgoCFD>(hit_alg_pset, rise_alg_pset);
    else
      throw art::Exception(art::errors::UnimplementedFeature)
        << "Cannot find implementation for " << threshAlgName << " algorithm.\n";

    // Initialize the pedestal estimation algorithm
    auto const ped_alg_pset = pset.get<fhicl::ParameterSet>("PedAlgoPset");
    std::string pedAlgName  = ped_alg_pset.get<std::string>("Name");
    if (pedAlgName == "Edges")
      fPedAlg = new pmtana::PedAlgoEdges(ped_alg_pset);
    else if (pedAlgName == "RollingMean")
      fPedAlg = new pmtana::PedAlgoRollingMean(ped_alg_pset);
    else if (pedAlgName == "UB")
      fPedAlg = new pmtana::PedAlgoUB(ped_alg_pset);
    else
      throw art::Exception(art::errors::UnimplementedFeature)
        << "Cannot find implementation for " << pedAlgName << " algorithm.\n";

    produces<std::vector<recob::OpHit>>();

    fPulseRecoMgr.AddRecoAlgo(fThreshAlg);
    fPulseRecoMgr.SetDefaultPedAlgo(fPedAlg);
  }

  //----------------------------------------------------------------------------
  // Destructor
  OpHitFinderDeco::~OpHitFinderDeco()
  {
    delete fThreshAlg;
    delete fPedAlg;
  }

  //----------------------------------------------------------------------------
  void OpHitFinderDeco::produce(art::Event& evt)
  {
    // These is the storage pointer we will put in the event
    std::unique_ptr<std::vector<recob::OpHit>> HitPtr(new std::vector<recob::OpHit>);
    std::vector<const sim::BeamGateInfo*> beamGateArray;
    try {
      evt.getView(fGenModule, beamGateArray);
    }

    catch (art::Exception const& err) {if (err.categoryCode() != art::errors::ProductNotFound) throw;}

    auto const& geometry(*lar::providerFrom<geo::Geometry>());
    auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    auto const& calibrator(*fCalib);

    if (fInputDigiType == "recob"){
      // Load pulses into WaveformVector
      art::Handle<std::vector<recob::OpWaveform>> decoHandle;
      art::Handle<std::vector<raw::OpDetWaveform>> rawHandle;
      
      evt.getByLabel(fInputModule,   decoHandle);
      evt.getByLabel(fInputModuledigi,rawHandle);
      
      assert(decoHandle.isValid());
      assert(rawHandle.isValid());
   
      RunHitFinder_deco(*decoHandle,
                        *rawHandle,
                        *HitPtr,
                        fPulseRecoMgr,
                        *fThreshAlg,
                        geometry,
                        fHitThreshold,
                        fScale,
                        clock_data,
                        calibrator,
                        fUseStartTime);
      }

    if (fInputDigiType == "raw"){
      // Load pulses into WaveformVector
      art::Handle<std::vector<raw::OpDetWaveform>> rawHandle;
      evt.getByLabel(fInputModuledigi,rawHandle);
      assert(rawHandle.isValid());
       if (fChannelMasks.empty() && fInputLabels.size() < 2) {
      art::Handle<std::vector<raw::OpDetWaveform>> rawHandle;
      if (fInputLabels.empty())
        evt.getByLabel(fInputModuledigi, rawHandle);
      else
        evt.getByLabel(fInputModuledigi, fInputLabels.front(), rawHandle);
      assert(rawHandle.isValid());
      RunHitFinder(*rawHandle,
                   *HitPtr,
                   fPulseRecoMgr,
                   *fThreshAlg,
                   geometry,
                   fHitThreshold,
                   clock_data,
                   calibrator,
                   fUseStartTime); 
    }
    else {

      // Reserve a large enough array
      int totalsize = 0;
      for (auto label : fInputLabels) {
        art::Handle<std::vector<raw::OpDetWaveform>> rawHandle;
        evt.getByLabel(fInputModuledigi, label, rawHandle);
        if (!rawHandle.isValid()) continue; // Skip non-existent collections
        totalsize += rawHandle->size();
      }

      std::vector<raw::OpDetWaveform> WaveformVector;
      WaveformVector.reserve(totalsize);

      for (auto label : fInputLabels) {
        art::Handle<std::vector<raw::OpDetWaveform>> rawHandle;
        evt.getByLabel(fInputModuledigi, label, rawHandle);
        if (!rawHandle.isValid()) continue; // Skip non-existent collections

        for (auto const& wf : *rawHandle) {
          if (fChannelMasks.find(wf.ChannelNumber()) != fChannelMasks.end()) continue;
          WaveformVector.push_back(wf);
        }
      }

      RunHitFinder(WaveformVector,
                   *HitPtr,
                   fPulseRecoMgr,
                   *fThreshAlg,
                   geometry,
                   fHitThreshold,
                   clock_data,
                   calibrator,
                   fUseStartTime);
      }
    }

    else{std::cout<<"InputDigiType not recognized"<<std::endl;}
    
    // Store results into the event
    std:: cout << "hit: " << HitPtr->size() << std::endl;
    evt.put(std::move(HitPtr));
  }   
} // namespace opdet
