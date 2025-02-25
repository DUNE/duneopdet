// ===============================================================
// OpFlashAnaDeco_module.cc
// This module is based on the larana/OpFlashAna_module.cc.  
// This analyzer writes out a TTree containing the properties of
// each reconstructed ophit and opflash
// ================================================================

// -*- mode: c++; c-basic-offset: 2; -*-
// This analyzer writes out a TTree containing the properties of
// each reconstructed flash
//

#ifndef OpFlashAnaDeco_h
#define OpFlashAnaDeco_h

// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

// C++ includes
#include "math.h"
#include <cstring>

// LArSoft includes
#include "larcore/Geometry/WireReadout.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"

// ART includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "fhiclcpp/ParameterSet.h"

namespace duneopdet {

  class OpFlashAnaDeco : public art::EDAnalyzer {
  public:
    // Standard constructor and destructor for an ART module.
    OpFlashAnaDeco(const fhicl::ParameterSet&);

    // The analyzer routine, called once per event.
    void analyze(const art::Event&);

  private:
    // The stuff below is the part you'll most likely have to change to
    // go from this custom example to your own task.

    // The parameters we'll read from the .fcl file.
    std::string fOpFlashModuleLabel; // Input tag for OpFlash collection
    std::string fOpHitModuleLabel;   // Input tag for OpHit collection
    float fSampleFreq;               // in MHz
    float fTimeBegin;                // in us
    float fTimeEnd;                  // in us

    float fYMin, fYMax, fZMin, fZMax;

    int PosHistYRes, PosHistZRes;

    bool fMakeFlashTimeHist;
    bool fMakeFlashPosHist;
    bool fMakePerFlashHists;

    bool fMakePerEventFlashTree;
    bool fMakePerFlashTree;
    bool fMakePerOpHitTree;
    bool fMakeFlashBreakdownTree;
    bool fMakeFlashHitMatchTree;

    TTree* fPerEventFlashTree;
    TTree* fPerFlashTree;
    TTree* fPerOpHitTree;
    TTree* fFlashBreakdownTree;
    TTree* fFlashHitMatchTree;

    Int_t fEventID;
    Int_t fFlashID;
    Int_t fHitID;
    Double_t fFlashTime;
    Double_t fAbsTime;
    bool fInBeamFrame;
    int fOnBeamTime;
    Float_t fTotalPE;
    Int_t fFlashFrame;

    Float_t fNPe;
    Float_t fYCenter;
    Float_t fYWidth;
    Float_t fZCenter;
    Float_t fZWidth;

    Int_t fOpChannel;
    Double_t fPeakTimeAbs;
    Double_t fPeakTime;
    Int_t fFrame;
    Float_t fWidth;
    Float_t fArea;
    Float_t fAmplitude;
    Float_t fPE;
    Float_t fFastToTotal;

    int fNFlashes;
    std::vector<int> fFlashIDVector;
    std::vector<float> fYCenterVector;
    std::vector<float> fZCenterVector;
    std::vector<float> fYWidthVector;
    std::vector<float> fZWidthVector;
    std::vector<double> fFlashTimeVector;
    std::vector<double> fAbsTimeVector;
    std::vector<int> fFlashFrameVector;
    std::vector<bool> fInBeamFrameVector;
    std::vector<int> fOnBeamTimeVector;
    std::vector<float> fTotalPEVector;
    int fNChannels;
    std::vector<float> fPEsPerFlashPerChannelVector;
    std::string fInputrecobOpflash;
  };

}

#endif

namespace duneopdet {

  //-----------------------------------------------------------------------
  // Constructor
  OpFlashAnaDeco::OpFlashAnaDeco(fhicl::ParameterSet const& pset) : EDAnalyzer(pset)
  {

    // Indicate that the Input Module comes from .fcl
    fOpFlashModuleLabel = pset.get<std::string>("OpFlashModuleLabel");
    fOpHitModuleLabel = pset.get<std::string>("OpHitModuleLabel");
    fInputrecobOpflash = pset.get<std::string>("InputrecobOpflash");
    
    auto const clock_data =
      art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
    fTimeBegin = clock_data.OpticalClock().Time();
    fTimeEnd = clock_data.OpticalClock().FramePeriod();
    fSampleFreq = clock_data.OpticalClock().Frequency();

    fYMin = pset.get<float>("YMin");
    fYMax = pset.get<float>("YMax");
    fZMin = pset.get<float>("ZMin");
    fZMax = pset.get<float>("ZMax");

    fMakeFlashTimeHist = pset.get<bool>("MakeFlashTimeHist");
    fMakeFlashPosHist = pset.get<bool>("MakeFlashPosHist");
    fMakePerFlashHists = pset.get<bool>("MakePerFlashHists");

    fMakePerEventFlashTree = pset.get<bool>("MakePerEventFlashTree");
    fMakePerFlashTree = pset.get<bool>("MakePerFlashTree");
    fMakePerOpHitTree = pset.get<bool>("MakePerOpHitTree");
    fMakeFlashBreakdownTree = pset.get<bool>("MakeFlashBreakdownTree");
    fMakeFlashHitMatchTree = pset.get<bool>("MakeFlashHitMatchTree");

    PosHistYRes = 100;
    PosHistZRes = 100;

    art::ServiceHandle<art::TFileService const> tfs;
    
    if (fMakePerOpHitTree) {
      fPerOpHitTree = tfs->make<TTree>("PerOpHitTree", "PerOpHitTree");
      fPerOpHitTree->Branch("EventID", &fEventID, "EventID/I");
      fPerOpHitTree->Branch("HitID", &fHitID, "HitID/I");
      fPerOpHitTree->Branch("OpChannel", &fOpChannel, "OpChannel/I");
      fPerOpHitTree->Branch("PeakTimeAbs", &fPeakTimeAbs, "PeakTimeAbs/D");
      fPerOpHitTree->Branch("PeakTime", &fPeakTime, "PeakTime/D");
      fPerOpHitTree->Branch("Frame", &fFrame, "Frame/I");
      fPerOpHitTree->Branch("Width", &fWidth, "Width/F");
      fPerOpHitTree->Branch("Area", &fArea, "Area/F");
      fPerOpHitTree->Branch("Amplitude", &fAmplitude, "Amplitude/F");
      fPerOpHitTree->Branch("PE", &fPE, "PE/F");
      fPerOpHitTree->Branch("FastToTotal", &fFastToTotal, "FastToTotal/F");
    }
    
    if (fInputrecobOpflash == "true"){

      if (fMakeFlashBreakdownTree) {
        fFlashBreakdownTree = tfs->make<TTree>("FlashBreakdownTree", "FlashBreakdownTree");
        fFlashBreakdownTree->Branch("EventID", &fEventID, "EventID/I");
        fFlashBreakdownTree->Branch("FlashID", &fFlashID, "FlashID/I");
        fFlashBreakdownTree->Branch("OpChannel", &fOpChannel, "OpChannel/I");
        fFlashBreakdownTree->Branch("FlashTime", &fFlashTime, "FlashTime/D");
        fFlashBreakdownTree->Branch("NPe", &fNPe, "NPe/F");
        fFlashBreakdownTree->Branch("AbsTime", &fAbsTime, "AbsTime/D");
        fFlashBreakdownTree->Branch("FlashFrame", &fFlashFrame, "FlashFrame/I");
        fFlashBreakdownTree->Branch("InBeamFrame", &fInBeamFrame, "InBeamFrame/B");
        fFlashBreakdownTree->Branch("OnBeamTime", &fOnBeamTime, "OnBeamTime/I");
        fFlashBreakdownTree->Branch("YCenter", &fYCenter, "YCenter/F");
        fFlashBreakdownTree->Branch("ZCenter", &fZCenter, "ZCenter/F");
        fFlashBreakdownTree->Branch("YWidth", &fYWidth, "YWidth/F");
        fFlashBreakdownTree->Branch("ZWidth", &fZWidth, "ZWidth/F");
        fFlashBreakdownTree->Branch("TotalPE", &fTotalPE, "TotalPE/F");
     }

      if (fMakePerFlashTree) {
        fPerFlashTree = tfs->make<TTree>("PerFlashTree", "PerFlashTree");
        fPerFlashTree->Branch("EventID", &fEventID, "EventID/I");
        fPerFlashTree->Branch("FlashID", &fFlashID, "FlashID/I");
        fPerFlashTree->Branch("YCenter", &fYCenter, "YCenter/F");
        fPerFlashTree->Branch("ZCenter", &fZCenter, "ZCenter/F");
        fPerFlashTree->Branch("YWidth", &fYWidth, "YWidth/F");
        fPerFlashTree->Branch("ZWidth", &fZWidth, "ZWidth/F");
        fPerFlashTree->Branch("FlashTime", &fFlashTime, "FlashTime/D");
        fPerFlashTree->Branch("AbsTime", &fAbsTime, "AbsTime/D");
        fPerFlashTree->Branch("FlashFrame", &fFlashFrame, "FlashFrame/I");
        fPerFlashTree->Branch("InBeamFrame", &fInBeamFrame, "InBeamFrame/B");
        fPerFlashTree->Branch("OnBeamTime", &fOnBeamTime, "OnBeamTime/I");
        fPerFlashTree->Branch("TotalPE", &fTotalPE, "TotalPE/F");
      }

      if (fMakePerEventFlashTree) {
        fPerEventFlashTree = tfs->make<TTree>("PerEventFlashTree", "PerEventFlashTree");
        fPerEventFlashTree->Branch("EventID", &fEventID, "EventID/I");
        fPerEventFlashTree->Branch("NFlashes", &fNFlashes, "NFlashes/I");
        fPerEventFlashTree->Branch("FlashIDVector", &fFlashIDVector);
        fPerEventFlashTree->Branch("YCenterVector", &fYCenterVector);
        fPerEventFlashTree->Branch("ZCenterVector", &fZCenterVector);
        fPerEventFlashTree->Branch("YWidthVector", &fYWidthVector);
        fPerEventFlashTree->Branch("ZWidthVector", &fZWidthVector);
        fPerEventFlashTree->Branch("FlashTimeVector", &fFlashTimeVector);
        fPerEventFlashTree->Branch("AbsTimeVector", &fAbsTimeVector);
        fPerEventFlashTree->Branch("FlashFrameVector", &fFlashFrameVector);
        fPerEventFlashTree->Branch("InBeamFrameVector", &fInBeamFrameVector);
        fPerEventFlashTree->Branch("OnBeamTimeVector", &fOnBeamTimeVector);
        fPerEventFlashTree->Branch("TotalPEVector", &fTotalPEVector);
        fPerEventFlashTree->Branch("NChannels", &fNChannels, "NChannels/I");
        // The only way I can think of to record a two-dimensional variable-size array in a TTree
        // is by flattening it into a one-dimension variable-size array
        fPerEventFlashTree->Branch("PEsPerFlashPerChannelVector", &fPEsPerFlashPerChannelVector);
      }

      if (fMakeFlashHitMatchTree) {
        fFlashHitMatchTree = tfs->make<TTree>("FlashHitMatchTree", "FlashHitMatchTree");
        fFlashHitMatchTree->Branch("EventID", &fEventID, "EventID/I");
        fFlashHitMatchTree->Branch("FlashID", &fFlashID, "FlashID/I");
        fFlashHitMatchTree->Branch("HitID", &fHitID, "HitID/I");
        fFlashHitMatchTree->Branch("OpChannel", &fOpChannel, "OpChannel/I");
        fFlashHitMatchTree->Branch("HitPeakTimeAbs", &fPeakTimeAbs, "HitPeakTimeAbs/F");
        fFlashHitMatchTree->Branch("HitPeakTime", &fPeakTime, "HitPeakTime/F");
        fFlashHitMatchTree->Branch("HitPE", &fPE, "HitPE/F");
        fFlashHitMatchTree->Branch("FlashPE", &fTotalPE, "FlashPE/F");
        fFlashHitMatchTree->Branch("FlashTimeAbs", &fAbsTime, "FlashTimeAbs/D");
        fFlashHitMatchTree->Branch("FlashTime", &fFlashTime, "FlashTime/D");
        fFlashHitMatchTree->Branch("HitFrame", &fFrame, "HitFrame/I");
        fFlashHitMatchTree->Branch("FlashFrame", &fFlashFrame, "FlashFrame/I");
      }

      fFlashID = 0;
    }  
  }

  //-----------------------------------------------------------------------
 void OpFlashAnaDeco::analyze(const art::Event& evt)
  {
  
  // Access ART's TFileService, which will handle creating and writing
  // histograms for us.
  art::ServiceHandle<art::TFileService const> tfs;
  fEventID = evt.id().event();

  if (fInputrecobOpflash == "true"){

    // Get flashes from event
    art::Handle<std::vector<recob::OpFlash>> FlashHandle;
    evt.getByLabel(fOpFlashModuleLabel, FlashHandle);

    // Get assosciations between flashes and hits
    art::FindManyP<recob::OpHit> Assns(FlashHandle, evt, fOpFlashModuleLabel);

    // Create string for histogram name
    char HistName[50];

    fFlashID = 0;
       
    std::vector<TH1D*> FlashHist;

    sprintf(HistName, "Event %d Flash Times", evt.id().event());
    TH1D* FlashTimes = nullptr;
    if (fMakeFlashTimeHist) {
      FlashTimes = tfs->make<TH1D>(HistName,
                                   ";t (ns);",
                                   int((fTimeEnd - fTimeBegin) * fSampleFreq),
                                   fTimeBegin * 1000.,
                                   fTimeEnd * 1000.);
    }

    TH2D* FlashPositions = nullptr;
    if (fMakeFlashPosHist) {
      sprintf(HistName, "Event %d All Flashes YZ", evt.id().event());

      FlashPositions =
        tfs->make<TH2D>(HistName, ";y ;z ", PosHistYRes, fYMin, fYMax, PosHistZRes, fZMin, fZMax);
    }

    unsigned int NOpChannels = art::ServiceHandle<geo::WireReadout const>()->Get().NOpChannels();

    if (fMakePerEventFlashTree) {
      fNFlashes = FlashHandle->size();
      fNChannels = NOpChannels;
    }

    // For every OpFlash in the vector
    for (unsigned int i = 0; i < FlashHandle->size(); ++i) {

      // Get OpFlash
      art::Ptr<recob::OpFlash> TheFlashPtr(FlashHandle, i);
      recob::OpFlash TheFlash = *TheFlashPtr;

      fFlashTime = TheFlash.Time();
      fFlashID = i; //++;

      TH2D* ThisFlashPosition = nullptr;
      if (fMakePerFlashHists) {
        sprintf(HistName, "Event %d t = %f", evt.id().event(), fFlashTime);
        FlashHist.push_back(
          tfs->make<TH1D>(HistName, ";OpChannel;PE", NOpChannels, 0, NOpChannels));

        sprintf(HistName, "Event %d Flash %f YZ", evt.id().event(), fFlashTime);

        ThisFlashPosition =
          tfs->make<TH2D>(HistName, ";y ;z ", PosHistYRes, fYMin, fYMax, PosHistZRes, fZMin, fZMax);
      }
      fYCenter = TheFlash.YCenter();
      fZCenter = TheFlash.ZCenter();
      fYWidth = TheFlash.YWidth();
      fZWidth = TheFlash.ZWidth();
      fInBeamFrame = TheFlash.InBeamFrame();
      fOnBeamTime = TheFlash.OnBeamTime();
      fAbsTime = TheFlash.AbsTime();
      fFlashFrame = TheFlash.Frame();
      fTotalPE = TheFlash.TotalPE();

      if (fMakePerEventFlashTree) {
        fFlashIDVector.emplace_back(i);
        fYCenterVector.emplace_back(TheFlash.YCenter());
        fZCenterVector.emplace_back(TheFlash.ZCenter());
        fYWidthVector.emplace_back(TheFlash.YWidth());
        fZWidthVector.emplace_back(TheFlash.ZWidth());
        fFlashTimeVector.emplace_back(TheFlash.Time());
        fAbsTimeVector.emplace_back(TheFlash.AbsTime());
        fFlashFrameVector.emplace_back(TheFlash.Frame());
        fInBeamFrameVector.emplace_back(TheFlash.InBeamFrame());
        fOnBeamTimeVector.emplace_back(TheFlash.OnBeamTime());
        fTotalPEVector.emplace_back(TheFlash.TotalPE());
      }

      for (unsigned int j = 0; j != NOpChannels; ++j) {
        if (fMakePerFlashHists) FlashHist.at(FlashHist.size() - 1)->Fill(j, TheFlash.PE(j));
        fNPe = TheFlash.PE(j);
        fOpChannel = j;

        if (fMakePerEventFlashTree) fPEsPerFlashPerChannelVector.emplace_back(TheFlash.PE(j));

        if ((fMakeFlashBreakdownTree) && (fNPe > 0)) fFlashBreakdownTree->Fill();
      }

      for (int y = 0; y != PosHistYRes; ++y)
        for (int z = 0; z != PosHistZRes; ++z) {
          float ThisY = fYMin + (fYMax - fYMin) * float(y) / PosHistYRes + 0.0001;
          float ThisZ = fZMin + (fZMax - fZMin) * float(z) / PosHistZRes + 0.0001;
          if (fMakePerFlashHists)
            ThisFlashPosition->Fill(ThisY,
                                    ThisZ,
                                    fTotalPE * exp(-pow((ThisY - fYCenter) / fYWidth, 2) / 2. -
                                                   pow((ThisZ - fZCenter) / fZWidth, 2) / 2.));
          if (fMakeFlashPosHist)
            FlashPositions->Fill(ThisY,
                                 ThisZ,
                                 fTotalPE * exp(-pow((ThisY - fYCenter) / fYWidth, 2) -
                                                pow((ThisZ - fZCenter) / fZWidth, 2)));
        }

      if (fMakeFlashTimeHist) FlashTimes->Fill(fFlashTime, fTotalPE);

      if (fMakePerFlashTree) fPerFlashTree->Fill();

      if (fMakeFlashHitMatchTree) {
        // Extract the assosciations
        fHitID = 0;
        std::vector<art::Ptr<recob::OpHit>> matchedHits = Assns.at(i);
        for (auto ophit : matchedHits) {
          fOpChannel = ophit->OpChannel();
          fPeakTimeAbs = ophit->PeakTimeAbs();
          fPeakTime = ophit->PeakTime();
          fFrame = ophit->Frame();
          fWidth = ophit->Width();
          fArea = ophit->Area();
          fAmplitude = ophit->Amplitude();
          fPE = ophit->PE();
          fFastToTotal = ophit->FastToTotal();
          fFlashHitMatchTree->Fill();
          fHitID++;
        }
      }
    }

    
    if (fMakePerEventFlashTree) {
      fPerEventFlashTree->Fill();
      fFlashIDVector.clear();
      fYCenterVector.clear();
      fZCenterVector.clear();
      fYWidthVector.clear();
      fZWidthVector.clear();
      fFlashTimeVector.clear();
      fAbsTimeVector.clear();
      fFlashFrameVector.clear();
      fInBeamFrameVector.clear();
      fOnBeamTimeVector.clear();
      fTotalPEVector.clear();
      fPEsPerFlashPerChannelVector.clear();
    }
  }  
    if (fMakePerOpHitTree) {
      art::Handle<std::vector<recob::OpHit>> OpHitHandle;
      evt.getByLabel(fOpHitModuleLabel, OpHitHandle);

      for (size_t i = 0; i != OpHitHandle->size(); ++i) {
        fOpChannel = OpHitHandle->at(i).OpChannel();
        fPeakTimeAbs = OpHitHandle->at(i).PeakTimeAbs();
        fPeakTime = OpHitHandle->at(i).PeakTime();
        fFrame = OpHitHandle->at(i).Frame();
        fWidth = OpHitHandle->at(i).Width();
        fArea = OpHitHandle->at(i).Area();
        fAmplitude = OpHitHandle->at(i).Amplitude();
        fPE = OpHitHandle->at(i).PE();
        fFastToTotal = OpHitHandle->at(i).FastToTotal();
        fHitID = i;
        fPerOpHitTree->Fill();
      }
    }
  }

} // namespace duneopdet

namespace duneopdet {
  DEFINE_ART_MODULE(OpFlashAnaDeco)
}