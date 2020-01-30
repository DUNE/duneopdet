///////////////////////////////////////////////////////////////
// This script generate TDR plot for calorimetry study.
//////////////////////////////////////////////////////////////

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "THStack.h"
#include "TLegend.h"

void plot_tdr_sn_calE(){

  TFile* in_true_21pe = new TFile("/dune/app/users/bbehera/develop/2019/PD/tdr/fout_resoln_sn_true_21pe.root", "read");
  TFile* in_file_phy  = new TFile("/dune/app/users/bbehera/develop/2018/PD/tdr/fout_physics.root", "read");
  TFile* in_file_tpc  = new TFile("/dune/app/users/bbehera/develop/2019/PD/tdr/RMSvsEnergyPlots_std.root", "read");

  TH1F*   hpds      = (TH1F*)in_true_21pe->Get("hresn_x_pecut_21pe");
  TH1D*   hphysics  = (TH1D*)in_file_phy->Get("hphysics");
  TGraph* htpc      = (TGraph*)in_file_tpc->Get("StdDev_ChargeDCR");

  gStyle->SetOptStat(0);  

  TCanvas* c1 = new TCanvas("c1","c1");
  
  hphysics->GetXaxis()->SetRangeUser(0, 30);
  htpc->GetXaxis()->SetRangeUser(0, 30);
  hpds->GetXaxis()->SetRangeUser(0, 30);
  hpds->GetYaxis()->SetRangeUser(0, 0.4);
  hpds->SetLineWidth(3);
  hphysics->SetLineColor(kRed);
  hphysics->SetLineWidth(3); 
  htpc->SetLineWidth(3);
  hpds->SetLineColor(kBlue);

  hpds->Draw("e hist ");
  htpc->Draw("CL same");
  hphysics->Draw("hist ][ same");

  TLegend *leg = AutoPlaceLegend(0.3, 0.2, -1);
  leg->SetTextSize(0.04);
  leg->AddEntry(htpc,"TPC resolution","l");
  leg->AddEntry(hpds,"PDS resolution","l");
  leg->AddEntry(hphysics,"Physics limited resolution","l");
  leg->Draw(); 

  c1->SaveAs("plots/pds-snb-res-vs-truex_optimized.png");
  c1->SaveAs("plots/pds-snb-res-vs-truex_optimized.pdf");
}
