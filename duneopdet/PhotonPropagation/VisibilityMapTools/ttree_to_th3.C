/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : ttree_to_th3.C
 * @created     : 2024-10-18 08:13
 */

#include <cstdio>
#include <iostream>
#include <fstream>
#include <string.h>
#include <sstream>

#include "RtypesCore.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include "TDirectoryFile.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TClonesArray.h"
#include "TH1D.h"
#include "TH3F.h"
#include <THn.h>
#include <THnSparse.h>
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TROOT.h"

void ttree_to_th3(const TString input_file, const TString vis_dir, const Double_t vis_threshold = 1e-4) 
{
  gStyle->SetTitleSize(0.06, "xyz");
  gStyle->SetLabelSize(0.06, "xyz");

  TFile* file = new TFile(input_file); 
  TTree* tOpDet  = nullptr; 
  TTree* tVisTPC = nullptr; 
  TTree* tVisBuf = nullptr; 
  TDirectoryFile* dunevis_dir = file->Get<TDirectoryFile>(vis_dir); 
  dunevis_dir->cd(); 

  TH1D* h_tmp[3] = {0}; 
  h_tmp[0] = dunevis_dir->Get<TH1D>("hgrid0"); 
  h_tmp[1] = dunevis_dir->Get<TH1D>("hgrid1"); 
  h_tmp[2] = dunevis_dir->Get<TH1D>("hgrid2"); 

  const int nbins[3] = {h_tmp[0]->GetNbinsX(), h_tmp[1]->GetNbinsX(), h_tmp[2]->GetNbinsX()};
  const double xmin[3] ={
    h_tmp[0]->GetXaxis()->GetXmin(), 
    h_tmp[1]->GetXaxis()->GetXmin(), 
    h_tmp[2]->GetXaxis()->GetXmin()};
  const double xmax[3] ={
    h_tmp[0]->GetXaxis()->GetXmax(), 
    h_tmp[1]->GetXaxis()->GetXmax(), 
    h_tmp[2]->GetXaxis()->GetXmax()}; 

  tOpDet = dunevis_dir->Get<TTree>("opDetMap"); 
  tVisTPC = dunevis_dir->Get<TTree>("photoVisMap"); 
  tVisBuf = dunevis_dir->Get<TTree>("photoVisMapBuffer");

  TTreeReader reader(tVisTPC); 
  TTreeReaderArray<double> _xyz(reader, "coords"); 
  TTreeReaderArray<double> _opdet_visd(reader, "opDet_visDirect"); 
  //TTreeReaderArray<double> _opdet_visr(reader, "opDet_visReflct"); 
  TTreeReaderValue<double> _vis_d(reader, "total_visDirect"); 
  //TTreeReaderValue<double> _vis_r(reader, "total_visReflct"); 

  const Long64_t N_OPDET = tOpDet->GetEntries();

  THnF* h3total_vis = new THnF("h3VisMap", "Visibility Map", 3, nbins, xmin, xmax);
  std::vector<THnSparseF*> h3opDet; 
  h3opDet.reserve(N_OPDET);
  for (size_t iopdet = 0; iopdet < N_OPDET; iopdet++) {
    TString name = Form("h3VisMap_opDet%ld", iopdet); 
    TString titl = Form("Visibility Map - OpDet %ld", iopdet);
    h3opDet.emplace_back(new THnSparseF(name, titl, 3, nbins, xmin, xmax));
  }

  TH1D* h_opdet_vis_tpc = new TH1D("h_opdet_vis_tpc", "opdet visibility (TPC)", 100, -15, 1); 
  TH1D* h_opdet_vis_buf = new TH1D("h_opdet_vis_buf", "opdet visibility (buffer)", 100, -15, 1); 

  while( reader.Next() ) {
    double xyz[3] = {_xyz[0], _xyz[1], _xyz[2]};
    h3total_vis->Fill(xyz, *(_vis_d)); 
    for (Long64_t i = 0; i < N_OPDET; i++) {
      h_opdet_vis_tpc->Fill( log10(_opdet_visd[i]) ); 
      if (_opdet_visd[i] > vis_threshold) {
        h3opDet[i]->Fill(xyz, _opdet_visd[i]); 
      }
    }
  }

  if (tVisBuf) {
    if (tVisBuf->GetEntries() > 0) {
      TTreeReader readerBuff(tVisBuf); 
      TTreeReaderArray<double> _buf_xyz(readerBuff, "coords"); 
      TTreeReaderArray<double> _buf_opdet_visd(readerBuff, "opDet_visDirectBuff"); 
      //TTreeReaderArray<double> _buf_opdet_visr(readerBuff, "opDet_visReflctBuff"); 
      TTreeReaderValue<double> _buf_vis_d(readerBuff, "total_visDirectBuff"); 
      //TTreeReaderValue<double> _buf_vis_r(readerBuff, "total_visReflctBuff"); 

      while( readerBuff.Next() ) {
        double xyz[3] = {_buf_xyz[0], _buf_xyz[1], _buf_xyz[2]};
        h3total_vis->Fill(xyz, *(_buf_vis_d)); 
        for (Long64_t i = 0; i < N_OPDET; i++) {
          h_opdet_vis_buf->Fill( log10(_buf_opdet_visd[i]) ); 
          if ( _buf_opdet_visd[i] > vis_threshold ) {
            h3opDet[i]->Fill(xyz, _buf_opdet_visd[i]); 
          }
        }
      }
    }
  }

  TString output_name = input_file; 
  output_name.Resize( output_name.Index(".root") ); 
  output_name += "_"+vis_dir+".root"; 

  TFile* output_h3 = new TFile(output_name, "recreate"); 
  h3total_vis->Write(); 
  for (auto& h3 : h3opDet) {
    h3->Write(); 
  }
  output_h3->Close(); 

  return;
}

TH3F* rebin_visibility_map(const TH3F* h3, const int rbx, const int rby, const int rbz)
{
  TH3F* h3rb = (TH3F*)h3->Clone(); 
  h3rb = (TH3F*)h3rb->Rebin3D(rbx, rby, rbz, Form("%s_rb_%i_%i_%i", h3->GetName(), rbx, rby, rbz)); 

  float scale_factor = 1.0 / (rbx * rby *rbz); 
  h3rb->Scale( scale_factor ); 
  
  return h3rb; 
}


