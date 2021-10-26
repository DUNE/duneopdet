// -*- mode: c++; c-basic-offset: 2; -*-
// This analyzer writes out a folder containing the 
// PDS event display. For every OpFLsh, the OpHits are 
// plotted. 
//
//


#ifndef ProtoDUNE_opdet_eventdisplay_H
#define ProtoDUNE_opdet_eventdisplay_H 

// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TTree.h"
#include "TLine.h"
#include "TLegend.h"

// C++ includes
#include <map>
#include <vector>
#include <iostream>
#include <cstring>
#include <sstream>
#include "math.h"
#include <climits>
#include <functional>
#include <algorithm>
#include <map>


// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "lardataobj/RawData/OpDetPulse.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "lardataalg/DetectorInfo/DetectorClocks.h"
//#include "OpticalDetector/OpDigiProperties.h"

#include "larana/OpticalDetector/OpFlashAlg.h"
#include <numeric>

// ART includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"

namespace opdet {
 
  class ProtoDUNE_opdet_eventdisplay : public art::EDAnalyzer{
  public:
 
    // Standard constructor and destructor for an ART module.
    ProtoDUNE_opdet_eventdisplay(const fhicl::ParameterSet&);
    virtual ~ProtoDUNE_opdet_eventdisplay();

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
    float fSampleFreq;                     // in MHz
    float fTimeBegin;                      // in us
    float fTimeEnd;                        // in us
    
    float fYMin, fYMax, fZMin, fZMax;

    int PosHistYRes, PosHistZRes;

//    bool fMakeFlashTimeHist;
//    bool fMakeFlashPosHist;
//    bool fMakePerFlashHists;
//    bool fMakePerOpHitTree;
          
//    TTree * fPerOpHitTree;
   
    Int_t   fEventID;
    Int_t   fFlashID;
    Int_t   fHitID;
    Float_t fFlashTime; 
//    Float_t fAbsTime;
//    bool    fInBeamFrame;
//    int     fOnBeamTime;
//    Float_t fTotalPE;
//    Int_t   fFlashFrame;
    
//    Float_t fNPe;
//    Float_t fYCenter;
//    Float_t fYWidth;
//    Float_t fZCenter;
//    Float_t fZWidth;

    Int_t   fOpChannel;
    Float_t fPeakTimeAbs;
    Float_t fPeakTime;
//    Int_t   fFrame;
//    Float_t fWidth;
//    Float_t fArea;
//    Float_t fAmplitude;
    Float_t fPE;
//    Float_t fFastToTotal;
    Float_t fYOpCenter;
    Float_t fZOpCenter;
    Float_t fXOpCenter;

//    int fNFlashes;
    std::vector< int >   fFlashIDVector;
    std::vector< float > fYCenterVector;
    std::vector< float > fZCenterVector;
    std::vector< float > fYWidthVector;
    std::vector< float > fZWidthVector;
    std::vector< float > fFlashTimeVector;
    std::vector< float > fAbsTimeVector;
    std::vector< int >   fFlashFrameVector;
    std::vector< bool >  fInBeamFrameVector;
    std::vector< int >   fOnBeamTimeVector;
    std::vector< float > fTotalPEVector;
//    int fNChannels;
    std::vector< float > fPEsPerFlashPerChannelVector;
  };

} 

#endif // ProtoDUNE_opdet_eventdisplay_H

namespace opdet {

  //-----------------------------------------------------------------------
  // Constructor
  ProtoDUNE_opdet_eventdisplay::ProtoDUNE_opdet_eventdisplay(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
  {

    // Indicate that the Input Module comes from .fcl
    fOpFlashModuleLabel = pset.get<std::string>("OpFlashModuleLabel");
    fOpHitModuleLabel   = pset.get<std::string>("OpHitModuleLabel");

//    art::ServiceHandle<OpDigiProperties> odp;
//    fTimeBegin  = odp->TimeBegin();
//    fTimeEnd    = odp->TimeEnd();
//    fSampleFreq = odp->SampleFreq();

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
    fTimeBegin  = clockData.OpticalClock().Time();
    fTimeEnd    = clockData.OpticalClock().FramePeriod();
    fSampleFreq = clockData.OpticalClock().Frequency();
   
    fYMin = pset.get<float>("YMin");
    fYMax = pset.get<float>("YMax");
    fZMin = pset.get<float>("ZMin");
    fZMax = pset.get<float>("ZMax");
    
    PosHistYRes = 20;
    PosHistZRes = 42;

    art::ServiceHandle< art::TFileService > tfs;
    
    fFlashID = 0;
  }

  //-----------------------------------------------------------------------
  // Destructor
  ProtoDUNE_opdet_eventdisplay::~ProtoDUNE_opdet_eventdisplay() 
  {}
   
  //-----------------------------------------------------------------------
  void ProtoDUNE_opdet_eventdisplay::beginJob()
  {}

  //-----------------------------------------------------------------------
  void ProtoDUNE_opdet_eventdisplay::analyze(const art::Event& evt) 
  {

    float displace[288]={0};

     for(int j=0;j<287;j++) displace[j]=j;
     for(int j=0;j<284;j+=4){
       displace[j]=-100;
       displace[j+1]=-40;
       displace[j+2]=5;
       displace[j+3]=60;
       if((j>=132 && j<=143) || (j>=264 && j<=275)){
	 displace[j]=-100;
      	 displace[j+1]=-75;
       	 displace[j+2]=-60;
       	 displace[j+3]=-40;
       	 displace[j+4]=-20;
       	 displace[j+5]=-5;
       	 displace[j+6]=5;
       	 displace[j+7]=20;
	 displace[j+8]=40;
	 displace[j+9]=60;
	 displace[j+10]=75;
	 displace[j+11]=100;
	 j+=8;
       }
       
    }

      Float_t binsz[] = {0,10, 29, 48, 66, 84, 102, 120, 138, 156, 174, 192, 211, 230, 240, 257, 278, 296, 314, 332, 350, 368, 386, 404, 422, 441, 460, 470, 489, 508, 526, 544, 562, 580, 598, 616, 634, 650, 671, 690, 700};
      Int_t  binnumz = sizeof(binsz)/sizeof(Float_t) - 1;

      Float_t binsy[]={-600,-540, -480, -420, -360, -300, -240, -180, -120, -60, 0, 60, 120, 180, 240, 300, 360, 420, 480, 540, 600};
      Int_t  binnumy = sizeof(binsy)/sizeof(Float_t) - 1;
    

    
    // Get flashes from event
      auto FlashHandle = evt.getHandle< std::vector< recob::OpFlash > >(fOpFlashModuleLabel);

    // Get assosciations between flashes and hits
    art::FindManyP< recob::OpHit > Assns(FlashHandle, evt, fOpFlashModuleLabel);

    // Create string for histogram name
    char HistName[50];
    
    fFlashID = 0;

    // Access ART's TFileService, which will handle creating and writing
    // histograms for us.
        art::ServiceHandle< art::TFileService > tfs;
    
    //std::vector<TH1D*> FlashHist;
    std::vector<TH2D*> DispHist;
    
    fEventID = evt.id().event();

    sprintf(HistName, "Event %d Disps", evt.id().event()); 
    
    art::ServiceHandle< geo::Geometry > geom;
    //unsigned int NOpChannels = geom->NOpChannels(); 

    // For every OpFlash in the vector
    for(unsigned int i = 0; i < FlashHandle->size(); ++i)
    {

      // Get OpFlash
      art::Ptr< recob::OpFlash > TheFlashPtr(FlashHandle, i);
      recob::OpFlash TheFlash = *TheFlashPtr;

      fFlashTime = TheFlash.Time();
      fFlashID   = i; //++;
    

      TH2D * ThisEventDisp = nullptr;    
      ThisEventDisp = tfs->make<TH2D>(HistName,"PE from OpHits in OpFlash;z ;DaS             RaS (y)",
				      binnumz, binsz,
				      binnumy, binsy);
    

      double xyz[3];
      auto OpHitHandle = evt.getHandle< std::vector< recob::OpHit > >(fOpHitModuleLabel);
      
      for(size_t i=0; i!=OpHitHandle->size(); ++i)
	{
	  fOpChannel   = OpHitHandle->at(i).OpChannel();
	  geom->OpDetGeoFromOpChannel(fOpChannel).GetCenter(xyz);
	  fXOpCenter   = xyz[0];
	  fYOpCenter   = xyz[1];
	      fZOpCenter   = xyz[2];
	      fPeakTimeAbs = OpHitHandle->at(i).PeakTimeAbs();
	      fPeakTime    = OpHitHandle->at(i).PeakTime();
	      fPE          = OpHitHandle->at(i).PE();
	      fHitID       = i;
	      //fPerOpHitTree->Fill();
	      //ThisEventDisp->Fill(fZOpCenter, fYOpCenter,fPE);
	      if(fXOpCenter >0){
		if(fOpChannel>=132 && fOpChannel<=143)  
		  ThisEventDisp->Fill(fZOpCenter+displace[fOpChannel], fYOpCenter,fPE);
		else{
		  ThisEventDisp->Fill(fZOpCenter+displace[fOpChannel], fYOpCenter,fPE);
		  ThisEventDisp->Fill(fZOpCenter+displace[fOpChannel]+18, fYOpCenter,fPE);
		  ThisEventDisp->Fill(fZOpCenter+displace[fOpChannel]+36, fYOpCenter,fPE);
		}
	      }
	      else if(fXOpCenter<0){
		if(fOpChannel>=264 && fOpChannel<=275)  
		  ThisEventDisp->Fill(fZOpCenter+displace[fOpChannel], fYOpCenter-600,fPE);
		else{
		  ThisEventDisp->Fill(fZOpCenter+displace[fOpChannel], fYOpCenter-600,fPE);
		  ThisEventDisp->Fill(fZOpCenter+displace[fOpChannel]+18, fYOpCenter-600,fPE);
		  ThisEventDisp->Fill(fZOpCenter+displace[fOpChannel]+38, fYOpCenter-600,fPE);
		}
	      }
	      TLine *line = new TLine(0,0,600,0);
	      line->SetLineColor(kBlack);
	      line->Draw();
	    }
    




    }
  

  }

} // namespace opdet

namespace opdet {
  DEFINE_ART_MODULE(ProtoDUNE_opdet_eventdisplay)
}
