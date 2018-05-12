/*
 * \file: HitSdpPlotter.h
 * \author: Jason Stock (jason.stock@mines.sdsmt.edu
 * \brief: This is a small Tree made for use in calibration analysis.
 *
 */

#ifndef TOY_H
#define TOY_H

//Includes
#include "dune/DuneObjBase/EventRecord.h"

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

      art::InputTag private_Label;
      art::ServiceHandle<art::TFileService> private_service_tfs;
      TH2D* private_Hist;


      //Also need buffers

  };// end class HitSdpPlotter

}

#endif//endif HitSdpPlotter
