////////////////////////////////////////////////////////////////////////////////////
// Class:       SolarOpFlash                                                      //
// Module Type: producer                                                          //
// File:        SolarOpFlash_module.cc                                            //
//                                                                                //
// OpHit clusterer, based on time and proximity information.                      //
// Optimized for low energy point-like interactions.                              //
// Using Dan Pershey's OpSlicer_module.cc as template.                            //
////////////////////////////////////////////////////////////////////////////////////

#ifndef SolarOpFlash_H
#define SolarOpFlash_H 1

#include "duneana/SolarNuAna/SolarAuxUtils.h"
#include "duneana/SolarNuAna/AdjOpHitsUtils.h"

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

namespace solar
{

  class SolarOpFlash : public art::EDProducer
  {
  public:
    void ProduceOpFlash(const std::vector<std::vector<int>> &ohitsidx,
                        const std::vector<AdjOpHitsUtils::FlashInfo> &oflashinfo,
                        std::vector<recob::OpFlash> &oflashes,
                        std::vector<std::vector<int>> &assoc,
                        detinfo::DetectorClocksData const &ts) const;
    // Standard constructor and destructor for an ART module.
    explicit SolarOpFlash(const fhicl::ParameterSet &);
    virtual ~SolarOpFlash();

    void beginJob();
    void endJob();
    void reconfigure(fhicl::ParameterSet const &pset);

    // The producer routine, called once per event.
    void produce(art::Event &);

  private:
    // The parameters we'll read from the .fcl file.
    art::ServiceHandle<geo::Geometry> geo;
    std::string fGeometry;
    std::string fOpHitLabel; // Input tag for OpHit collection
    int fOpFlashAlgoNHit;
    double fDetectorSizeX;
    float fOpFlashAlgoTime;
    float fOpFlashAlgoRad;
    float fOpFlashAlgoPE;
    float fOpFlashAlgoTriggerPE;
    // bool fOpFlashAlgoCentroid;
    std::unique_ptr<solar::SolarAuxUtils> solaraux;
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
        fGeometry(p.get<std::string>("Geometry")),
        fOpHitLabel(p.get<std::string>("OpHitLabel")),
        fOpFlashAlgoNHit(p.get<int>("OpFlashAlgoNHit")),
        fDetectorSizeX(p.get<double>("DetectorSizeX")),
        fOpFlashAlgoTime(p.get<float>("OpFlashAlgoTime")),
        fOpFlashAlgoRad(p.get<float>("OpFlashAlgoRad")),
        fOpFlashAlgoPE(p.get<float>("OpFlashAlgoPE")),
        fOpFlashAlgoTriggerPE(p.get<float>("OpFlashAlgoTriggerPE")),
        // fOpFlashAlgoCentroid(p.get<bool>("OpFlashAlgoCentroid")),
        solaraux(new solar::SolarAuxUtils(p)),
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

    // These are the storage pointers we will put in the event
    std::unique_ptr<std::vector<recob::OpFlash>>
        flashPtr(new std::vector<recob::OpFlash>);
    std::unique_ptr<art::Assns<recob::OpFlash, recob::OpHit>>
        assnPtr(new art::Assns<recob::OpFlash, recob::OpHit>);

    // This will keep track of what flashes will assoc to what ophits
    // at the end of processing
    std::vector<std::vector<int>> assocList;

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
    adjophits->CalcAdjOpHits(OpHitList, OpHitVec, OpHitIdx);
    adjophits->MakeFlashVector(FlashVec, OpHitVec, evt);
    ProduceOpFlash(OpHitIdx, FlashVec, *flashPtr, assocList, clockData);

    // Make the associations which we noted we need
    for (size_t i = 0; i != assocList.size(); ++i)
    {
      art::PtrVector<recob::OpHit> opHitPtrVector;
      for (int const &hitIndex : assocList.at(i))
      {
        art::Ptr<recob::OpHit> opHitPtr(opHitHandle, hitIndex);
        opHitPtrVector.push_back(opHitPtr);
      }

      util::CreateAssn(*this, evt, *flashPtr, opHitPtrVector,
                       *(assnPtr.get()), i);
    }

    // Store results into the event
    evt.put(std::move(flashPtr));
    evt.put(std::move(assnPtr));
  }

  //--------------------------------------------------------------------------

  void SolarOpFlash::ProduceOpFlash(const std::vector<std::vector<int>> &OpHitsIdx,
                                    const std::vector<AdjOpHitsUtils::FlashInfo> &OpFlashInfo,
                                    std::vector<recob::OpFlash> &oflashes,
                                    std::vector<std::vector<int>> &assoc,
                                    detinfo::DetectorClocksData const &ts) const
  {
    // Loop over the flashes with TheFlash being the flash
    for (int i = 0; i < int(OpFlashInfo.size()); i++)
    {
      AdjOpHitsUtils::FlashInfo OpFlash = OpFlashInfo[i];

      // Time to make the OpFlash collect the info from the opflashinfo struct
      double OpFlashFast2Tot = OpFlash.FastToTotal;
      double OpFlashY = OpFlash.Y;
      double OpFlashdY = OpFlash.YWidth;
      double OpFlashZ = OpFlash.Z;
      double OpFlashdZ = OpFlash.ZWidth;
      double OpFlashT = OpFlash.Time;
      double OpFlashdT = OpFlash.TimeWidth;
      std::vector<double> OpFlashPEs = OpFlash.PEperOpDet;

      // Loop over the hits in the flash
      std::vector<int> OpHitIdx;
      for (int j = 0; j < int(OpHitsIdx[i].size()); j++)
      {
        OpHitIdx.push_back(OpHitsIdx[i][j]);
      }

      // From OpFlashAlg
      int Frame = ts.OpticalClock().Frame(OpFlashT - 18.1);
      if (Frame == 0)
        Frame = 1;

      int BeamFrame = ts.OpticalClock().Frame(ts.TriggerTime());
      bool InBeamFrame = false;
      if (!(ts.TriggerTime() < 0))
        InBeamFrame = (Frame == BeamFrame);

      int OnBeamTime = 0;
      if (InBeamFrame)
        OnBeamTime = 1;

      oflashes.emplace_back(OpFlashT, OpFlashdT, OpFlashT, Frame,
                            OpFlashPEs, InBeamFrame, OnBeamTime, OpFlashFast2Tot,
                            OpFlashY, OpFlashdY, OpFlashZ, OpFlashdZ);
      assoc.emplace_back(OpHitIdx);
    }
  }
} // namespace solar
