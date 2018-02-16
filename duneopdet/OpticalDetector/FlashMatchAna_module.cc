// -*- mode: c++; c-basic-offset: 2; -*-
// This analyzer writes out a TTree for studying the matching
// between flashes and events
//

#ifndef FlashMatchAna_H
#define FlashMatchAna_H 1

// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TTree.h"

// C++ includes
#include <map>
#include <vector>
#include <iostream>
#include <cstring>
#include <sstream>
#include "math.h"
#include <climits>

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larsim/MCCheater/PhotonBackTracker.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "dune/OpticalDetector/OpFlashSort.h"

// ART includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"

namespace opdet {
 
  class FlashMatchAna : public art::EDAnalyzer{
  public:
 
    // Standard constructor and destructor for an ART module.
    FlashMatchAna(const fhicl::ParameterSet&);
    virtual ~FlashMatchAna();

    // This method is called once, at the start of the job. In this
    // example, it will define the histogram we'll write.
    void beginJob();

    // The analyzer routine, called once per event. 
    void analyze (const art::Event&); 

  private:

    // The stuff below is the part you'll most likely have to change to
    // go from this custom example to your own task.

    // The parameters we'll read from the .fcl file.
    std::string fOpFlashModuleLabel;       // Input tag for OpFlash collection
    std::string fOpHitModuleLabel;         // Input tag for OpHit collection
    std::string fSignalLabel;              // Input tag for the signal generator label
    std::string fGeantLabel;               // Input tag for GEANT
    
    TTree * fFlashMatchTree;

    Int_t fEventID;

    Float_t fTrueX;
    Float_t fTrueY;
    Float_t fTrueZ;
    Float_t fTrueT;
    Float_t fDetectedT;
    Float_t fTrueE;
    Int_t   fTruePDG;
    
    Int_t fNFlashes;    
    std::vector< Int_t >   fFlashIDVector;
    std::vector< Float_t > fYCenterVector;
    std::vector< Float_t > fZCenterVector;
    std::vector< Float_t > fYWidthVector;
    std::vector< Float_t > fZWidthVector;
    std::vector< Float_t > fTimeVector;
    std::vector< Float_t > fTimeWidthVector;
    std::vector< Float_t > fTimeDiffVector;
    std::vector< Int_t >   fFlashFrameVector;
    std::vector< Bool_t >  fInBeamFrameVector;
    std::vector< Int_t >   fOnBeamTimeVector;
    std::vector< Float_t > fTotalPEVector;
    std::vector< Bool_t >  fSignalVector;
    std::vector< Float_t > fPurityVector;
    Int_t fNOpDets;
    std::vector<Int_t> fNHitOpDetVector;
    //std::vector< std::vector<Float_t> > fPEsPerFlashPerOpDetVector;
  };

} 

#endif // FlashMatchAna_H

namespace opdet {

  //-----------------------------------------------------------------------
  // Constructor
  FlashMatchAna::FlashMatchAna(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
  {

    // Indicate that the Input Module comes from .fcl
    fOpFlashModuleLabel = pset.get<std::string>("OpFlashModuleLabel");
    fOpHitModuleLabel   = pset.get<std::string>("OpHitModuleLabel");
    fSignalLabel        = pset.get<std::string>("SignalLabel");
    fGeantLabel         = pset.get<std::string>("GeantLabel");

    art::ServiceHandle< art::TFileService > tfs;

    fFlashMatchTree = tfs->make<TTree>("FlashMatchTree","FlashMatchTree");
    fFlashMatchTree->Branch("EventID",                     &fEventID,   "EventID/I");
    fFlashMatchTree->Branch("TrueX",                       &fTrueX,     "TrueX/F");
    fFlashMatchTree->Branch("TrueY",                       &fTrueY,     "TrueY/F");
    fFlashMatchTree->Branch("TrueZ",                       &fTrueZ,     "TrueZ/F");
    fFlashMatchTree->Branch("TrueT",                       &fTrueT,     "TrueT/F");
    fFlashMatchTree->Branch("DetectedT",                   &fDetectedT, "DetectedT/F");
    fFlashMatchTree->Branch("TrueE",                       &fTrueE,     "TrueE/F");
    fFlashMatchTree->Branch("TruePDG",                     &fTruePDG,   "TruePDG/I");
    fFlashMatchTree->Branch("NFlashes",                    &fNFlashes,  "NFlashes/I");
    fFlashMatchTree->Branch("FlashIDVector",               &fFlashIDVector);
    fFlashMatchTree->Branch("YCenterVector",               &fYCenterVector);
    fFlashMatchTree->Branch("ZCenterVector",               &fZCenterVector);
    fFlashMatchTree->Branch("YWidthVector",                &fYWidthVector);
    fFlashMatchTree->Branch("ZWidthVector",                &fZWidthVector);
    fFlashMatchTree->Branch("TimeVector",                  &fTimeVector);
    fFlashMatchTree->Branch("TimeWidthVector",             &fTimeWidthVector);
    fFlashMatchTree->Branch("TimeDiffVector",              &fTimeDiffVector);
    fFlashMatchTree->Branch("TotalPEVector",               &fTotalPEVector);
    fFlashMatchTree->Branch("NOpDets",                     &fNOpDets,    "NOpDets/I");
    fFlashMatchTree->Branch("NHitOpDetVector",             &fNHitOpDetVector);
    //fFlashMatchTree->Branch("PEsPerFlashPerOpDetVector",   &fPEsPerFlashPerOpDetVector);
    fFlashMatchTree->Branch("Signal",                      &fSignalVector);
    fFlashMatchTree->Branch("Purity",                      &fPurityVector);
  }

  //-----------------------------------------------------------------------
  // Destructor
  FlashMatchAna::~FlashMatchAna() 
  {}
   
  //-----------------------------------------------------------------------
  void FlashMatchAna::beginJob()
  {}

  //-----------------------------------------------------------------------
  void FlashMatchAna::analyze(const art::Event& evt) 
  {
    // Get the required services
    auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    auto const* timeService = lar::providerFrom<detinfo::DetectorClocksService>();
    art::ServiceHandle< geo::Geometry > geom;
    art::ServiceHandle< cheat::PhotonBackTracker > pbt;
    art::ServiceHandle<cheat::ParticleInventoryService> pinv;
    art::ServiceHandle< art::TFileService > tfs;
    //pbt->Rebuild(evt);
    

    // Record the event ID
    fEventID = evt.id().event();


    //////////////////////////////////////
    // Access all the truth information //
    //////////////////////////////////////
    
    auto MClistHandle = evt.getValidHandle<std::vector<simb::MCTruth> >(fSignalLabel);

    art::Ptr<simb::MCTruth> mctruth(MClistHandle, 0);
    if (mctruth->NParticles() == 0) {
      mf::LogError("FlashMatchAna") << "No MCTruth Particles";
    }

    // Get all the track ids associated with the signal event.
    std::set<int> signal_trackids;
    art::FindManyP<simb::MCParticle> SignalGeantAssns(MClistHandle,evt,fGeantLabel);
    for ( size_t i = 0; i < SignalGeantAssns.size(); i++) {
      auto parts = SignalGeantAssns.at(i);
      for (auto part = parts.begin(); part != parts.end(); part++) {
        signal_trackids.emplace((*part)->TrackId());
      }
    }

    // Get just the neutrino, entry 0 from the list, and record its properties
    const simb::MCParticle& part(mctruth->GetParticle(0));
    fTrueX     = part.Vx();
    fTrueY     = part.Vy();
    fTrueZ     = part.Vz();
    fTrueT     = part.T()*1000; // ns -> us
    fTrueE     = part.E();
    fTruePDG   = part.PdgCode();

    // Get the PlaneID which describes the location of the true vertex
    int plane = 0;
    double loc[] = {part.Vx(), part.Vy(), part.Vz()};
    geo::TPCID tpc = geom->FindTPCAtPosition(loc);
    geo::PlaneID planeid(tpc, plane);

    // Convert true X to would-be charge arrival time, and convert from ticks to us, add to MC time
    double deltaTicks = detprop->ConvertXToTicks(part.Vx(), planeid);
    double deltaT = timeService->TPCTick2Time(deltaTicks);
    fDetectedT = fTrueT + deltaT;

    // Get the maximum possible time difference by getting number of ticks corresponding to
    // one full drift distance, and converting to time. 
    double maxT = timeService->TPCTick2Time(detprop->NumberTimeSamples());
    

    
    //////////////////////////////////////
    // Access all the Flash Information //
    //////////////////////////////////////

    // Get flashes from event
    art::Handle< std::vector< recob::OpFlash > > FlashHandle;
    std::vector<art::Ptr<recob::OpFlash> > flashlist;
    if (evt.getByLabel(fOpFlashModuleLabel, FlashHandle)) {
      art::fill_ptr_vector(flashlist, FlashHandle);
      std::sort(flashlist.begin(), flashlist.end(), recob::OpFlashPtrSortByPE);
    }
    
    // Get assosciations between flashes and hits
    art::FindManyP< recob::OpHit > Assns(FlashHandle, evt, fOpFlashModuleLabel);


    fNOpDets   = geom->NOpDets();
    fNFlashes  = FlashHandle->size();

    double maxPurity = 0;
    int maxPurityID = 0;
    

    
    // For every OpFlash in the vector
    for(unsigned int i = 0; i < FlashHandle->size(); ++i)
    {
      // Get OpFlash and associated hits
      art::Ptr< recob::OpFlash > TheFlashPtr(FlashHandle, i);
      recob::OpFlash TheFlash = *TheFlashPtr;
      std::vector< art::Ptr<recob::OpHit> > matchedHits = Assns.at(i);

      // Calculate the flash purity
      double purity = pbt->OpHitCollectionPurity(signal_trackids, matchedHits);

      if (purity > maxPurity) {
        maxPurity = purity;
        maxPurityID = i;
      }
      
      // Calcuate relative detection time
      double timeDiff = fDetectedT - TheFlash.Time();
      
      // Check if this is a possible flash (w/in 1 drift window)
      // Otherwise, don't store it
      if (timeDiff < -10 || timeDiff > maxT)
        continue;
      
      fFlashIDVector    .emplace_back(i);
      fYCenterVector    .emplace_back(TheFlash.YCenter());
      fZCenterVector    .emplace_back(TheFlash.ZCenter());
      fYWidthVector     .emplace_back(TheFlash.YWidth());
      fZWidthVector     .emplace_back(TheFlash.ZWidth());
      fTimeVector       .emplace_back(TheFlash.Time());
      fTimeWidthVector  .emplace_back(TheFlash.TimeWidth());
      fTimeDiffVector   .emplace_back(timeDiff);
      fTotalPEVector    .emplace_back(TheFlash.TotalPE());
      fSignalVector     .emplace_back(false);
      fPurityVector     .emplace_back(purity);

      std::vector<Float_t> PEPerOpDet;
      for(unsigned int iOD = 0; iOD < geom->NOpDets(); ++iOD){
        PEPerOpDet.emplace_back(0);
      }
      for(unsigned int iC=0; iC < geom->NOpChannels(); ++iC)
      {
        unsigned int iOD = geom->OpDetFromOpChannel(iC);
        PEPerOpDet[iOD] += TheFlash.PE(iC);
      }
      
      int NHitOpDets = 0;
      for(unsigned int iOD = 0; iOD < geom->NOpDets(); ++iOD){
        if (PEPerOpDet[iOD] > 0) ++NHitOpDets;
      }
      fNHitOpDetVector.emplace_back(NHitOpDets);
      //fPEsPerFlashPerOpDetVector.push_back(PEPerOpDet);
    }

    // Mark the flash with the highest purity as signal
    auto it = std::find(fFlashIDVector.begin(), fFlashIDVector.end(), maxPurityID);
    if (it != fFlashIDVector.end()) {
      auto index = std::distance(fFlashIDVector.begin(), it);
      fSignalVector[index] = true;
    }
    

    ////////////////////////////
    // Write out and clean up //
    ////////////////////////////
    
    fFlashMatchTree->Fill();
    fFlashIDVector              .clear();
    fYCenterVector              .clear();
    fZCenterVector              .clear();
    fYWidthVector               .clear();
    fZWidthVector               .clear();
    fTimeVector                 .clear();
    fTimeWidthVector            .clear();
    fTimeDiffVector             .clear();
    fTotalPEVector              .clear();
    fNHitOpDetVector            .clear();
    //fPEsPerFlashPerOpDetVector  .clear();
    fSignalVector               .clear();
    fPurityVector               .clear();
  }
} // namespace opdet

namespace opdet {
  DEFINE_ART_MODULE(FlashMatchAna)
}
