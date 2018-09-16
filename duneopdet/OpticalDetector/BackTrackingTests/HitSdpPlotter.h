/*
 * \file: HitSdpPlotter.h
 * \author: Jason Stock (jason.stock@mines.sdsmt.edu
 * \brief: This is a small Tree made for use in calibration analysis.
 *
 */

#ifndef TOY_H
#define TOY_H

//Includes
#include "dune/DuneObj/CalibTreeRecord.h"

//LArSoft Includes
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/Simulation/OpDetBacktrackerRecord.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/PhotonBackTrackerService.h"

//FrameworkIncludes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Utilities/Exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TH2D.h"

//CPP includes
#include <vector>
#include <map>


namespace{}//

namespace HitSdpPlotter {

  class HitSdpPlotter : public art::EDAnalyzer
  {
    public:

      //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      HitSdpPlotter(fhicl::ParameterSet const& pSet);

      virtual void beginJob();
      virtual void analyze(const art::Event& evt) override;
      virtual void endJob();


    private:

      art::InputTag private_OpHitLabel;
      art::InputTag private_BtrLabel;
      art::ServiceHandle<art::TFileService> private_service_tfs;
      art::ServiceHandle<geo::Geometry> private_service_geom;
      art::ServiceHandle<cheat::PhotonBackTrackerService> private_service_pbt;
      std::unique_ptr<art::TFileDirectory> mDir;
      std::map<UInt_t, TH1D*> sdp_time_hists;
      std::map<UInt_t, TH1D*> detTimeHist;
      TH1D* widths;


      //Also need buffers

  };// end class HitSdpPlotter

}

#endif//endif HitSdpPlotter
