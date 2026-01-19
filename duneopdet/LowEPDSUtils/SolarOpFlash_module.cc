////////////////////////////////////////////////////////////////////////////////////
// Class:       SolarOpFlash                                                      //
// Module Type: producer                                                          //
// Module Label solarflash                                                        //
// File:        SolarOpFlash_module.cc                                            //
//                                                                                //
// OpHit clusterer, based on time and proximity information.                      //
// Optimized for low energy point-like interactions.                              //
// Using Dan Pershey's OpSlicer_module.cc as template.                            //
////////////////////////////////////////////////////////////////////////////////////

#ifndef SolarOpFlash_H
#define SolarOpFlash_H 1

#include "duneopdet/LowEPDSUtils/AdjOpHitsUtils.h"
#include "dunecore/ProducerUtils/ProducerUtils.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorClocks.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

// C++ Includes
#include <vector>
#include <string>
#include <memory>
#include <limits>

using namespace producer;

namespace solar
{

  class SolarOpFlash : public art::EDProducer
  {
  public:
    // Standard constructor and destructor for an ART module.
    explicit SolarOpFlash(const fhicl::ParameterSet &);
    virtual ~SolarOpFlash();

    void beginJob();
    void endJob();
    void reconfigure(fhicl::ParameterSet const &pset);

    // The producer routine, called once per event.
    void produce(art::Event &);

    void ProduceOpFlash(const std::vector<AdjOpHitsUtils::FlashInfo> &oflashinfo,
                        std::vector<recob::OpFlash> &oflashes,
                        detinfo::DetectorClocksData const &ts) const;
  private:
    // The parameters we'll read from the .fcl file.
    art::ServiceHandle<geo::Geometry> geo;
    std::string fOpWaveformLabel; // Input tag for OpDetWaveform collection
    std::string fOpHitLabel; // Input tag for OpHit collection
    std::string fOpHitTimeVariable;
    int fOpFlashAlgoNHit;
    float fOpFlashAlgoMinTime;
    float fOpFlashAlgoMaxTime;
    bool fOpFlashAlgoWeightedTime;
    float fOpFlashAlgoRad;
    float fOpFlashAlgoPE;
    float fOpFlashAlgoTriggerPE;
    float fOpFlashTimeOffset;
    std::unique_ptr<producer::ProducerUtils> producer;
    std::unique_ptr<solar::AdjOpHitsUtils> adjophits;
  };

}

namespace solar
{
  DEFINE_ART_MODULE(SolarOpFlash)
}

#endif

namespace solar
{

  //--------------------------------------------------------------------------
  // Constructor
  SolarOpFlash::SolarOpFlash(const fhicl::ParameterSet &p)
      : EDProducer{p},
        fOpHitLabel(p.get<std::string>("OpHitLabel", "ophitspe")),
        fOpHitTimeVariable(p.get<std::string>("OpHitTimeVariable", "PeakTime")),
        fOpFlashAlgoNHit(p.get<int>("OpFlashAlgoNHit", 3)),
        fOpFlashAlgoMinTime(p.get<double>("OpFlashAlgoMinTime", 0.32)), //  [20 PDS tick]
        fOpFlashAlgoMaxTime(p.get<double>("OpFlashAlgoMaxTime", 0.96)), //  [60 PDS tick]
        fOpFlashAlgoRad(p.get<float>("OpFlashAlgoRad", 500.0)),
        fOpFlashAlgoPE(p.get<float>("OpFlashAlgoPE", 1.5)),
        fOpFlashAlgoTriggerPE(p.get<float>("OpFlashAlgoTriggerPE", 1.5)),
        fOpFlashTimeOffset(p.get<float>("OpFlashTimeOffset", 0.0)),
        producer(new producer::ProducerUtils(p)),
        adjophits(new solar::AdjOpHitsUtils(p))
  {
    reconfigure(p);
    produces<std::vector<recob::OpFlash>>();
    produces<art::Assns<recob::OpFlash, recob::OpHit>>();
  }

  //--------------------------------------------------------------------------
  void SolarOpFlash::reconfigure(fhicl::ParameterSet const &p)
  {
  }

  //--------------------------------------------------------------------------
  // Destructor
  SolarOpFlash::~SolarOpFlash()
  {
  }
  //--------------------------------------------------------------------------
  void SolarOpFlash::beginJob()
  {
  }

  //--------------------------------------------------------------------------
  void SolarOpFlash::endJob()
  {
  }

  //--------------------------------------------------------------------------
  void SolarOpFlash::produce(art::Event &evt)
  {
    ProducerUtils::PrintInColor("SolarOpFlash::produce called", ProducerUtils::GetColor("green"), "Debug");
    // These are the storage pointers we will put in the event
    std::unique_ptr<std::vector<recob::OpFlash>>
        flashPtr(new std::vector<recob::OpFlash>);
    std::unique_ptr<art::Assns<recob::OpFlash, recob::OpHit>>
        assnPtr(new art::Assns<recob::OpFlash, recob::OpHit>);

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);

    // Get OpHits from the event
    auto opHitHandle = evt.getHandle<std::vector<recob::OpHit>>(fOpHitLabel);

    std::vector<art::Ptr<recob::OpHit>> OpHitList;
    std::vector<std::vector<art::Ptr<recob::OpHit>>> OpHitVec;
    std::vector<std::vector<int>> OpHitIdx;
    std::vector<AdjOpHitsUtils::FlashInfo> FlashVec;
    for (int i = 0; i < int(opHitHandle->size()); i++)
    {
      art::Ptr<recob::OpHit> opHitPtr(opHitHandle, i);
      OpHitList.push_back(opHitPtr);
    }

    // Run the clustering
    adjophits->CalcAdjOpHits(OpHitList, OpHitVec, OpHitIdx, evt);
    adjophits->MakeFlashVector(FlashVec, OpHitVec, evt);
    ProduceOpFlash(FlashVec, *flashPtr, clockData);

    // Make the associations which we noted we need
    for (size_t i = 0; i != OpHitIdx.size(); ++i)
    {
      art::PtrVector<recob::OpHit> opHitPtrVector;
      for (int const &hitIndex : OpHitIdx.at(i))
      {
        // art::Ptr<recob::OpHit> opHitPtr(opHitHandle, hitIndex);
        opHitPtrVector.push_back(OpHitList.at(hitIndex));
      }
      if (i < 10)
      {
        ProducerUtils::PrintInColor("Generating OpFlash " + ProducerUtils::str(i) + " with " + ProducerUtils::str(opHitPtrVector.size()) + " hits", ProducerUtils::GetColor("green"), "Debug");
      }
      // Create the association between the flash and the OpHits
      util::CreateAssn(*this, evt, *flashPtr, opHitPtrVector,
                       *(assnPtr.get()), i);
      if (i == OpHitIdx.size() - 1)
      {
        ProducerUtils::PrintInColor("Generated " + ProducerUtils::str(i + 1) + " OpFlashes", ProducerUtils::GetColor("green"), "Debug");
      }
    }
    // Store results into the event
    evt.put(std::move(flashPtr));
    evt.put(std::move(assnPtr));
  }

  //--------------------------------------------------------------------------

  void SolarOpFlash::ProduceOpFlash(const std::vector<AdjOpHitsUtils::FlashInfo> &OpFlashInfo,
                                    std::vector<recob::OpFlash> &oflashes,
                                    detinfo::DetectorClocksData const &ts) const
  {
    // Loop over the flashes with TheFlash being the flash
    for (int i = 0; i < int(OpFlashInfo.size()); i++)
    {
      AdjOpHitsUtils::FlashInfo OpFlash = OpFlashInfo[i];

      // Time to make the OpFlash collect the info from the opflashinfo struct
      double OpFlashFast2Tot = OpFlash.FastToTotal;
      double OpFlashX = OpFlash.X;
      double OpFlashdX = OpFlash.XWidth;
      double OpFlashY = OpFlash.Y;
      double OpFlashdY = OpFlash.YWidth;
      double OpFlashZ = OpFlash.Z;
      double OpFlashdZ = OpFlash.ZWidth;
      double OpFlashT = -1e6;
      if (fOpFlashAlgoWeightedTime) {
        OpFlashT = OpFlash.TimeWeighted - fOpFlashTimeOffset; // Convert to time to us happens in MakeFlashVector
      }
      else {
        OpFlashT = OpFlash.MainOpHitTime - fOpFlashTimeOffset; // Convert to time to us happens in MakeFlashVector
      }
      double OpFlashdT = OpFlash.TimeWidth; // Convert to time to us happens in MakeFlashVector
      std::vector<double> OpFlashPEs = OpFlash.PEperOpDet;

      // From OpFlashAlg
      // Use Frame to save the OpFlash Plane -> more usefull for analysis
      // Use InBeamFrame to save if the OpFlash is in a correct Plane
      int Plane = -1;
      bool OnPlane = false;
      if (OpFlash.Plane > -1){
        Plane = OpFlash.Plane;
        OnPlane = true;
      }
      else
      {
        // If the plane is not set, we use the OpFlashT to determine the frame
        Plane = 9999; // Default unsigned value for no plane
      }

      int Frame = ts.OpticalClock().Frame(OpFlashT); // Hard coded 18.1 us offset moved to fcl
      int BeamFrame = ts.OpticalClock().Frame(ts.TriggerTime());
      bool InBeamFrame = false;
      if (!(ts.TriggerTime() < 0))
        InBeamFrame = (Frame == BeamFrame);

      int OnBeamTime = 0;
      if (InBeamFrame)
        OnBeamTime = 1;

      oflashes.emplace_back(OpFlashT, OpFlashdT, OpFlashT, Plane, // Using Frame to save the OpFlash Plane
                            OpFlashPEs, OnPlane, OnBeamTime, OpFlashFast2Tot,
                            OpFlashX, OpFlashdX, OpFlashY, OpFlashdY, OpFlashZ, OpFlashdZ);
    }
  }
} // namespace solar
