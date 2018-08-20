// -*- mode: c++; c-basic-offset: 2; -*-
// This analyzer writes out a TTree for studying the matching
// between flashes and events
//

#ifndef FlashMatchAna_H
#define FlashMatchAna_H 1

// ROOT includes
#include "TH1.h"
#include "TEfficiency.h"
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
#include "larsim/MCCheater/PhotonBackTrackerService.h"
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
    TTree * fLargestFlashTree;
    TTree * fSelectedFlashTree;

    TEfficiency * fRecoEfficiencyVsE;
    TEfficiency * fRecoEfficiencyVsX;
    TEfficiency * fRecoEfficiencyVsXandE;
    TEfficiency * fLargestEfficiencyVsE;
    TEfficiency * fLargestEfficiencyVsX;
    TEfficiency * fLargestEfficiencyVsXandE;
    TEfficiency * fSelectedEfficiencyVsE;
    TEfficiency * fSelectedEfficiencyVsX;
    TEfficiency * fSelectedEfficiencyVsXandE;

    // Parameters from the fhicl
    int   fNBinsE;
    float fLowE;
    float fHighE;
    int   fNBinsX;
    float fLowX;
    float fHighX;
    float fDistanceCut;

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
    std::vector< Float_t > fPurityVector;
    std::vector< Float_t > fDistanceVector;
    Int_t fNOpDets;
    std::vector<Int_t> fNHitOpDetVector;

    Int_t    fFlashID;
    Float_t  fYCenter;
    Float_t  fZCenter;
    Float_t  fYWidth;
    Float_t  fZWidth;
    Float_t  fTime;
    Float_t  fTimeWidth;
    Float_t  fTimeDiff;
    //Int_t    fFlashFrame; // unused
    //Bool_t   fInBeamFrame; // unused
    //Int_t    fOnBeamTime; // unused
    Float_t  fTotalPE;
    Float_t  fPurity;
    Float_t  fDistance;
    Int_t    fNHitOpDets;
    std::vector< Float_t > fPEsPerOpDetVector;

    
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
    fNBinsE             = pset.get<int>("NBinsE");
    fLowE               = pset.get<float>("LowE");
    fHighE              = pset.get<float>("HighE");
    fNBinsX             = pset.get<int>("NBinsX");
    fLowX               = pset.get<float>("LowX");
    fHighX              = pset.get<float>("HighX");
    fDistanceCut        = pset.get<float>("DistanceCut");

    

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
    fFlashMatchTree->Branch("Purity",                      &fPurityVector);
    fFlashMatchTree->Branch("Distance",                    &fDistanceVector);

    fLargestFlashTree = tfs->make<TTree>("LargestFlashTree","LargestFlashTree");
    fLargestFlashTree->Branch("EventID",                     &fEventID,   "EventID/I");
    fLargestFlashTree->Branch("TrueX",                       &fTrueX,     "TrueX/F");
    fLargestFlashTree->Branch("TrueY",                       &fTrueY,     "TrueY/F");
    fLargestFlashTree->Branch("TrueZ",                       &fTrueZ,     "TrueZ/F");
    fLargestFlashTree->Branch("TrueT",                       &fTrueT,     "TrueT/F");
    fLargestFlashTree->Branch("DetectedT",                   &fDetectedT, "DetectedT/F");
    fLargestFlashTree->Branch("TrueE",                       &fTrueE,     "TrueE/F");
    fLargestFlashTree->Branch("TruePDG",                     &fTruePDG,   "TruePDG/I");
    fLargestFlashTree->Branch("NFlashes",                    &fNFlashes,  "NFlashes/I");
    fLargestFlashTree->Branch("FlashID",                     &fFlashID,   "FlashID/I");
    fLargestFlashTree->Branch("YCenter",                     &fYCenter,   "YCenter/F");
    fLargestFlashTree->Branch("ZCenter",                     &fZCenter,   "ZCenter/F");
    fLargestFlashTree->Branch("YWidth",                      &fYWidth,    "YWidth/F");
    fLargestFlashTree->Branch("ZWidth",                      &fZWidth,    "ZWidth/F");
    fLargestFlashTree->Branch("Time",                        &fTime,      "Time/F");
    fLargestFlashTree->Branch("TimeWidth",                   &fTimeWidth, "TimeWidth/F");
    fLargestFlashTree->Branch("TimeDiff",                    &fTimeDiff,  "TimeDiff/F");
    fLargestFlashTree->Branch("TotalPE",                     &fTotalPE,   "TotalPE/F");
    fLargestFlashTree->Branch("NOpDets",                     &fNOpDets,   "NOpDets/I");
    fLargestFlashTree->Branch("NHitOpDets",                  &fNHitOpDets,"NHitOpDets/I");
    fLargestFlashTree->Branch("PEsPerOpDetVector",           &fPEsPerOpDetVector);
    fLargestFlashTree->Branch("Purity",                      &fPurity,    "Purity/F");
    fLargestFlashTree->Branch("Distance",                    &fDistance,  "Distance/F");


    fSelectedFlashTree = tfs->make<TTree>("SelectedFlashTree","SelectedFlashTree");
    fSelectedFlashTree->Branch("EventID",                     &fEventID,   "EventID/I");
    fSelectedFlashTree->Branch("TrueX",                       &fTrueX,     "TrueX/F");
    fSelectedFlashTree->Branch("TrueY",                       &fTrueY,     "TrueY/F");
    fSelectedFlashTree->Branch("TrueZ",                       &fTrueZ,     "TrueZ/F");
    fSelectedFlashTree->Branch("TrueT",                       &fTrueT,     "TrueT/F");
    fSelectedFlashTree->Branch("DetectedT",                   &fDetectedT, "DetectedT/F");
    fSelectedFlashTree->Branch("TrueE",                       &fTrueE,     "TrueE/F");
    fSelectedFlashTree->Branch("TruePDG",                     &fTruePDG,   "TruePDG/I");
    fSelectedFlashTree->Branch("NFlashes",                    &fNFlashes,  "NFlashes/I");
    fSelectedFlashTree->Branch("FlashID",                     &fFlashID,   "FlashID/I");
    fSelectedFlashTree->Branch("YCenter",                     &fYCenter,   "YCenter/F");
    fSelectedFlashTree->Branch("ZCenter",                     &fZCenter,   "ZCenter/F");
    fSelectedFlashTree->Branch("YWidth",                      &fYWidth,    "YWidth/F");
    fSelectedFlashTree->Branch("ZWidth",                      &fZWidth,    "ZWidth/F");
    fSelectedFlashTree->Branch("Time",                        &fTime,      "Time/F");
    fSelectedFlashTree->Branch("TimeWidth",                   &fTimeWidth, "TimeWidth/F");
    fSelectedFlashTree->Branch("TimeDiff",                    &fTimeDiff,  "TimeDiff/F");
    fSelectedFlashTree->Branch("TotalPE",                     &fTotalPE,   "TotalPE/F");
    fSelectedFlashTree->Branch("NOpDets",                     &fNOpDets,   "NOpDets/I");
    fSelectedFlashTree->Branch("NHitOpDets",                  &fNHitOpDets,"NHitOpDets/I");
    fSelectedFlashTree->Branch("PEsPerOpDetVector",           &fPEsPerOpDetVector);
    fSelectedFlashTree->Branch("Purity",                      &fPurity,    "Purity/F");
    fSelectedFlashTree->Branch("Distance",                    &fDistance,    "Distance/F");


    fRecoEfficiencyVsE         = tfs->make<TEfficiency>("recoEfficiencyVsE",         ";Energy (GeV);Efficiency",  fNBinsE, fLowE, fHighE);
    fRecoEfficiencyVsX         = tfs->make<TEfficiency>("recoEfficiencyVsX",         ";Position (cm);Efficiency", fNBinsX, fLowX, fHighX);
    fRecoEfficiencyVsXandE     = tfs->make<TEfficiency>("recoEfficiencyVsXandE",     ";Position (cm);Energy (GeV);Efficiency", fNBinsX, fLowX, fHighX, fNBinsE, fLowE, fHighE);
    fLargestEfficiencyVsE      = tfs->make<TEfficiency>("largestEfficiencyVsE",      ";Energy (GeV);Efficiency",  fNBinsE, fLowE, fHighE);
    fLargestEfficiencyVsX      = tfs->make<TEfficiency>("largestEfficiencyVsX",      ";Position (cm);Efficiency", fNBinsX, fLowX, fHighX);
    fLargestEfficiencyVsXandE  = tfs->make<TEfficiency>("largestEfficiencyVsXandE",  ";Position (cm);Energy (GeV);Efficiency", fNBinsX, fLowX, fHighX, fNBinsE, fLowE, fHighE);
    fSelectedEfficiencyVsE     = tfs->make<TEfficiency>("selectedEfficiencyVsE",     ";Energy (GeV);Efficiency",  fNBinsE, fLowE, fHighE);
    fSelectedEfficiencyVsX     = tfs->make<TEfficiency>("selectedEfficiencyVsX",     ";Position (cm);Efficiency", fNBinsX, fLowX, fHighX);
    fSelectedEfficiencyVsXandE = tfs->make<TEfficiency>("selectedEfficiencyVsXandE", ";Position (cm);Energy (GeV);Efficiency", fNBinsX, fLowX, fHighX, fNBinsE, fLowE, fHighE);
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
    art::ServiceHandle< cheat::PhotonBackTrackerService > pbt;
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
    if (! geom->HasTPC(tpc) ) {
      mf::LogInfo("FlashMatchAna") << "No valid TPC for " << tpc;
      return;
    }
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
    else {
      mf::LogError("FlashMatchAna") << "Cannot load any flashes. Failing";
      abort();
    }
    
    // Get assosciations between flashes and hits
    art::FindManyP< recob::OpHit > Assns(flashlist, evt, fOpFlashModuleLabel);


    // Set up some flags to fill as we loop
    // through flashes. These will control
    // filling of efficiency plots after the loop.
    bool AnyReconstructed = false;
    bool LargestFound     = false;
    bool LargestRight     = false;
    bool SelectedFound    = false;
    bool SelectedRight    = false;
    
    
    // For every OpFlash in the vector 
    fNOpDets   = geom->NOpDets();
    fNFlashes  = flashlist.size();
    for(unsigned int i = 0; i < flashlist.size(); ++i)
    {
      // Get OpFlash and associated hits
      recob::OpFlash TheFlash = *flashlist[i];
      std::vector< art::Ptr<recob::OpHit> > matchedHits = Assns.at(i);

      // Calculate the flash purity
      double purity = pbt->OpHitCollectionPurity(signal_trackids, matchedHits);
      
      // Calcuate relative detection time
      double timeDiff = fDetectedT - TheFlash.Time();
      
      // Check if this is a possible flash (w/in 1 drift window)
      // Otherwise, skip it
      if (timeDiff < -10 || timeDiff > maxT)
        continue;

      // Put flash info into variables
      fFlashID     = i;
      fYCenter     = TheFlash.YCenter();
      fZCenter     = TheFlash.ZCenter();
      fYWidth      = TheFlash.YWidth();
      fZWidth      = TheFlash.ZWidth();
      fTime        = TheFlash.Time();
      fTimeWidth   = TheFlash.TimeWidth();
      fTimeDiff    = timeDiff;
      fTotalPE     = TheFlash.TotalPE();
      fPurity      = purity;

      // Calculate distance from MC truth vertex in the Y-Z plane
      fDistance = sqrt( pow(fTrueY-fYCenter,2) +  pow(fTrueZ-fZCenter,2) );


      // Loop through all the opdets with hits in this flash
      fPEsPerOpDetVector.clear();
      for(unsigned int iOD = 0; iOD < geom->NOpDets(); ++iOD){
        fPEsPerOpDetVector.emplace_back(0);
      }
      for(unsigned int iC=0; iC < geom->NOpChannels(); ++iC)
      {
        unsigned int iOD = geom->OpDetFromOpChannel(iC);
        fPEsPerOpDetVector[iOD] += TheFlash.PE(iC);
      }
      
      fNHitOpDets = 0;
      for(unsigned int iOD = 0; iOD < geom->NOpDets(); ++iOD){
        if (fPEsPerOpDetVector[iOD] > 0) ++fNHitOpDets;
      }
      fNHitOpDetVector.emplace_back(fNHitOpDets);


      // Add flash info to the tree of all possible flashes
      fFlashIDVector    .emplace_back(fFlashID);
      fYCenterVector    .emplace_back(fYCenter);
      fZCenterVector    .emplace_back(fZCenter);
      fYWidthVector     .emplace_back(fYWidth);
      fZWidthVector     .emplace_back(fZWidth);
      fTimeVector       .emplace_back(fTime);
      fTimeWidthVector  .emplace_back(fTimeWidth);
      fTimeDiffVector   .emplace_back(fTimeDiff);
      fTotalPEVector    .emplace_back(fTotalPE);
      fPurityVector     .emplace_back(fPurity);
      fDistanceVector   .emplace_back(fDistance);


      // Did we reconstruct any flashes with signal in them?
      if (fPurity > 0) AnyReconstructed = true;

        
      // First == Largest, so if this is the first flash it is also the largest.
      // So, fill the LargestFlash tree and the LargestFlash efficiency plots
      if (!LargestFound) {

        // Write out the info into the tree for the largest flash
        fLargestFlashTree->Fill();

        // Record that we found the largest flash
        // and if we got it right
        LargestFound = true;
        if (fPurity > 0) LargestRight = true;
      }


      // The first time we get into here we have the largest flash that is
      // within the distance cut. So, fill the SelectedFlash tree and the
      // selected flash efficiency plots
      if (!SelectedFound && fDistance < fDistanceCut) {

        // Write out the info for the selected flash
        fSelectedFlashTree->Fill();

        // Record that we found the selected flash
        // and if we got it right
        SelectedFound = true;
        if (fPurity > 0) SelectedRight = true;
      }

    }

    // Fill these TEfficiencies once for every event
    // but use the booleans to decide if it was
    // "selected" or not.
    fRecoEfficiencyVsE->Fill(AnyReconstructed, fTrueE);
    fRecoEfficiencyVsX->Fill(AnyReconstructed, fTrueX);
    fRecoEfficiencyVsXandE->Fill(AnyReconstructed, fTrueX, fTrueE);

    fLargestEfficiencyVsE->Fill(LargestRight, fTrueE);
    fLargestEfficiencyVsX->Fill(LargestRight, fTrueX);
    fLargestEfficiencyVsXandE->Fill(LargestRight, fTrueX, fTrueE);
    
    fSelectedEfficiencyVsE->Fill(SelectedRight, fTrueE);
    fSelectedEfficiencyVsX->Fill(SelectedRight, fTrueX);
    fSelectedEfficiencyVsXandE->Fill(SelectedRight, fTrueX, fTrueE);



    ///////////////////////////////////////////////
    // Write out the FlashMatchTree and clean up //
    ///////////////////////////////////////////////
    
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
    fPurityVector               .clear();
    fDistanceVector             .clear();
  }
} // namespace opdet

namespace opdet {
  DEFINE_ART_MODULE(FlashMatchAna)
}
