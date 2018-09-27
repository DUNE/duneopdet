/*
 * \file: HitSdpPlotter_module.cc
 * \author: JStock (jason.stock@mines.sdsmt.edu)
 * \brief: This is a small analysis Tree made for use in the Calibration group to investigate radiological backgrounds and calibration sources.
 *
 */

#include "HitSdpPlotter.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include <vector>
#include <TObject.h>

namespace{}


namespace HitSdpPlotter {

  //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  HitSdpPlotter::HitSdpPlotter(fhicl::ParameterSet const& pSet)
    :EDAnalyzer(pSet),
    private_OpHitLabel(pSet.get<art::InputTag>("OpHitLabel", "ophit")),
    private_BtrLabel(pSet.get<art::InputTag>("BTRLabel", "largeant"))
  {  }

  //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  /*
     HitSdpPlotter::CalibrationTree(fhiclConfig const& config)
     :EDAnalyzer(config),
     private_HitLabel(config.HitLabel()),
     private_OpHitLabel(config.OpHitLabel())
     {  }
     */ //Commented out until I have a subtable for EDAnalyzer

  //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  void HitSdpPlotter::beginJob(){
//int nOpDets=private_geom->NOpDets();
    mDir = std::make_unique<art::TFileDirectory>(private_service_tfs->mkdir("hists","hists"));
    widths = private_service_tfs->make<TH1D>("ophit_widths", "width of ophits", 1000, -1, 2001);

    //PEs On X Axis
  }

  //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  void HitSdpPlotter::analyze(const art::Event& evt) {

    art::Handle<std::vector<recob::OpHit>> ophitHandle;
    std::vector<art::Ptr<recob::OpHit>> ophitList;
    if(evt.getByLabel(private_OpHitLabel, ophitHandle) )
      art::fill_ptr_vector(ophitList, ophitHandle);

    art::Handle<std::vector<sim::OpDetBacktrackerRecord>> btrHandle;
    std::vector<art::Ptr<sim::OpDetBacktrackerRecord>> btrList;
    if(evt.getByLabel(private_BtrLabel, btrHandle) )
      art::fill_ptr_vector(btrList, btrHandle);

    for( auto& btr : btrList ){
      TH1D* dummy = 0;
      UInt_t opdet = btr->OpDetNum();
      auto empIt = sdp_time_hists.emplace(std::make_pair(opdet, dummy));
      if(empIt.second==true){
        TString histName = TString("SDPTime")+TString(std::to_string(opdet))+TString("Hist");
        TString histTitle = TString("PD ")+TString(std::to_string(opdet))+TString(" SDP Times.");
        TH1D* tmphist1 = mDir->make<TH1D>(histName, histTitle, 30000, -1500000, 1500000);
        empIt.first->second = std::move(tmphist1);
      }
      auto time_sdps = btr->timePDclockSDPsMap();
      for( auto& time_sdp : time_sdps ){
        double time = time_sdp.first;
        empIt.first->second->Fill(time);
      }
    }

    for ( auto& ophit : ophitList ){
      TH1D* dummy = 0;
      UInt_t opdet = private_service_geom->OpDetFromOpChannel(ophit->OpChannel());
      auto empIt = detTimeHist.emplace(std::make_pair(opdet, dummy));
      //auto empIt = detTimeHist.emplace(std::make_pair(opdet, TH1D*(0)));
//      auto empLowIt = detLowerTimeHist.emplace(std::make_pair(opdet, dummy));
      //auto empLowIt = detLowerTimeHist.emplace(std::make_pair(opdet, TH1D*(0)));
//      auto empUpIt = detUpperTimeHist.emplace(std::make_pair(opdet, dummy));
      //auto empUpIt = detUpperTimeHist.emplace(std::make_pair(opdet, TH1D*(0)));
      if(empIt.second==true){
        TString histName = TString("PD")+TString(std::to_string(opdet))+TString("Hist");
        TString histTitle = TString("PD ")+TString(std::to_string(opdet))+TString(" OpHit Times.");
        /*
        TString histLowName = TString("PD")+TString(std::to_string(opdet))+TString("LowHist");
        TString histLowTitle = TString("PD ")+TString(std::to_string(opdet))+TString(" OpHit Lower Times.");
        TString histUpName = TString("PD")+TString(std::to_string(opdet))+TString("UpHist");
        TString histUpTitle = TString("PD ")+TString(std::to_string(opdet))+TString(" OpHit Upper Times.");
        */
        TH1D* tmphist1 = mDir->make<TH1D>(histName, histTitle, 30000, -1500000, 1500000);
        empIt.first->second = std::move(tmphist1);
        //TH1D* tmphist2 = mDir->make<TH1D>(histLowName, histLowTitle, 3000000, -1500000, 1500000);
        //empLowIt.first->second = std::move(tmphist2);
        //TH1D* tmphist3 = mDir->make<TH1D>(histUpName, histUpTitle, 3000000, -1500000, 1500000);
        //empUpIt.first->second = std::move(tmphist3);
      }
      if((private_service_pbt->OpHitToSimSDPs_Ps(ophit)).size()>1){
        double time = ophit->PeakTime() * 1000;
        double width = ophit->Width() * 1000;
        widths->Fill(width);
        empIt.first->second->Fill(time);
        //empLowIt.first->second->Fill(time - width);
        //empUpIt.first->second->Fill(time + width);
      }
    }


  } //end analyze

  //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  void HitSdpPlotter::endJob(){
  }

}//end namespace


DEFINE_ART_MODULE(HitSdpPlotter::HitSdpPlotter)
