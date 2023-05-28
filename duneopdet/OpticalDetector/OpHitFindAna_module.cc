// ===============================================================
// OpHitFindAna_module.cc
// This module is based on the larana/OpFlashAna_module.cc.  
// This analyzer writes out a TTree containing the properties of
// each reconstructed ophit
// ================================================================

#ifndef OpHitFindAna_h
#define OpHitFindAna_h

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpWaveform.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

// ART includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"

// ROOT includes
#include "TH1.h"
#include "THStack.h"
#include "TF1.h"
#include "TVector3.h"
#include "TTree.h"

// C++ Includes
#include <map>
#include <vector>
#include <cstring>
#include <sstream>
#include "math.h"


namespace opdet {

  class OpHitFindAna : public art::EDAnalyzer {
  public:
    // Standard constructor and destructor for an ART module.
    OpHitFindAna(const fhicl::ParameterSet&);
    virtual ~OpHitFindAna();
    // The analyzer routine, called once per event.
    void analyze(const art::Event&);

  private:
    
    // The parameters we'll read from the .fcl file.
    std::string fInputModule;    // Input tag for OpDet collection
    std::string fInstanceName;   // Input tag for OpDet collection
    double fSampleFreq;          // in MHz
    float fTimeBegin;            // in us
    float fTimeEnd;              // in us
     
    // Flags to enable or disable output of debugging 
    bool fMakePerOpHitTree;
         
    TTree* fPerOpHitTree;
    TTree* fPerChannelTree;

    // Output TTree and its branch variables
    
    Int_t fEventID;
    Int_t fHitID;
    Int_t fOpChannelID;
    Int_t fOpChannel;
    Double_t fPeakTimeAbs;
    Double_t fPeakTime;
    Int_t fFrame;
    Float_t fWidth;
    Float_t fArea;
    Float_t fAmplitude;
    Int_t fHit;
    Float_t fPE;
    Float_t fFastToTotal;
       
  };

}

#endif

 namespace opdet {

  DEFINE_ART_MODULE(OpHitFindAna)

 }


 namespace opdet {

  //-----------------------------------------------------------------------
  // Constructor
  OpHitFindAna::OpHitFindAna(fhicl::ParameterSet const& pset) : EDAnalyzer(pset)
  {

    // Indicate that the Input Module comes from .fcl
    fInputModule  =   pset.get< std::string >("InputModule");
    fInstanceName =   pset.get< std::string >("InstanceName");
    
    art::ServiceHandle<art::TFileService const> tfs;
    
    // Obtaining parameters from the DetectorClocksService
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
    fTimeBegin = clockData.OpticalClock().Time();
    fTimeEnd = clockData.OpticalClock().FramePeriod();
    fSampleFreq = clockData.OpticalClock().Frequency();
        
    fMakePerOpHitTree   = pset.get<bool>("MakePerOpHitTree");
   
    // TTree for output

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
        
  }
 
  //-----------------------------------------------------------------------
  // Destructor 
  OpHitFindAna::~OpHitFindAna() 
  {  
  }    
  
  //-----------------------------------------------------------------------
  void OpHitFindAna::analyze(const art::Event& evt)
  {

    // Create a handle for our vector of pulses
    art::Handle<std::vector<recob::OpHit>> OpHitHandle;
    
    // Read in HitHandle
    evt.getByLabel(fInputModule, OpHitHandle);

    // Access ART's TFileService, which will handle creating and writing
    art::ServiceHandle<art::TFileService const> tfs;

    fEventID = evt.id().event();
    
    art::ServiceHandle<geo::Geometry const> geom;
     
    
    // For every OpHit in the vector
    
    if (fMakePerOpHitTree) {
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
    
} // namespace opdet