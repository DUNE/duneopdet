/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : ttree_to_th3.C
 * @created     : 2024-10-18 08:13
 * @description : Fill 3D visibility maps from TTree produced by the PhotonVisibilityExport
 */

#include <cstdio>
#include <iostream>
#include <string.h>

#include <RtypesCore.h>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <TDirectoryFile.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TClonesArray.h>
#include <THn.h>
#include <THnSparse.h>
#include <TH1D.h>
#include <TROOT.h>

/** 
 * Create THn visibility maps from the PhotonVisibilityExport_module's output
 *
 * Create visibility maps for the entire detector volume and for each individual OpDet.
 * To reduce the resources needed for handling the visibility of the full set of 
 * OpDet, we use three-dimensional THnSparse to store the visiblity map, rounding 
 * the visibility to 0.0 when below the threshold set in vis_threshold.
 *
 * @param input_file input file path
 * @param visexport_label label used for the PhotonVisibilityExport module 
 * @param vis_threshold Set visibility to 0 if below this value. (Only for single OpDet maps)
 * @param split_fill Split the filling of OpDet visivility in N runs to reduce memory impact (default 5)
 */
void ttree_to_th3(const TString input_file, 
    const TString visexport_label, 
    const Double_t vis_threshold = 1e-4, 
    const UInt_t split_fill = 5) 
{
  gStyle->SetTitleSize(0.06, "xyz");
  gStyle->SetLabelSize(0.06, "xyz");

  TFile* file = new TFile(input_file); 
  TTree* tOpDet  = nullptr; 
  TTree* tVisTPC = nullptr; 
  TTree* tVisBuf = nullptr; 
  TDirectoryFile* dunevis_dir = file->Get<TDirectoryFile>(visexport_label); 
  dunevis_dir->cd(); 

  const TString axis_label[3] = {"x [cm]", "y [cm]", "z [cm]"}; 
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
  TTreeReaderArray<double> _opdet_visr(reader, "opDet_visReflct"); 
  TTreeReaderValue<double> _vis_d(reader, "total_visDirect"); 
  TTreeReaderValue<double> _vis_r(reader, "total_visReflct"); 

  const size_t N_OPDET = tOpDet->GetEntries();
  size_t nOpDetPerRun = N_OPDET;
  UInt_t nRuns = 1;
  if (split_fill > 1) {
    nRuns = split_fill;
    nOpDetPerRun = std::ceil( N_OPDET / (double)nRuns); 
  }
  std::vector<std::vector<size_t>> OpDetChunks(nRuns);
  size_t iopdet = 0; 
  for (size_t irun = 0; irun < nRuns; irun++) {
    std::vector<size_t>& opdets_ = OpDetChunks.at(irun); 
    opdets_.reserve(nOpDetPerRun); 
    while (iopdet < N_OPDET && opdets_.size() <= nOpDetPerRun) {
      opdets_.emplace_back(iopdet); 
      iopdet++; 
    }
  }

  TString output_name = input_file; 
  output_name.Resize( output_name.Index(".root") ); 
  output_name += "_"+visexport_label+".root"; 
  TFile* output_h3 = new TFile(output_name, "recreate"); 

  size_t irun = 0; 
  for (const auto& opdets : OpDetChunks) {
    printf("Processing OpDet chunk %ld - [%ld - %ld]\n", irun, opdets.front(), opdets.back()); 

    TClonesArray h3opDet("THnSparseF", opdets.size());
    size_t iopdet_idx = 0;
    for (const auto& iopdet : opdets) {
      TString name = Form("h3VisMap_opDet%ld", iopdet); 
      TString titl = Form("Visibility Map - OpDet %ld", iopdet);
      new (h3opDet[iopdet_idx]) THnSparseF(name, titl, 3, nbins, xmin, xmax);
      for (size_t idim = 0; idim < 3; idim++) {
        ((THnSparseF*)h3opDet[iopdet_idx])->GetAxis(idim)->SetTitle(axis_label[idim]); 
      }
      iopdet_idx++; 
    }

    THnF* hnTotal = nullptr; 
    if (irun == 0) {
      hnTotal = new THnF("h3VisMap", "Detector visibility map", 3, nbins, xmin, xmax);
      for (size_t idim = 0; idim < 3; idim++) {
        hnTotal->GetAxis(idim)->SetTitle(axis_label[idim]); 
      }
    }

    reader.Restart();

    while (reader.Next()) {
      double xyz[3] = {_xyz[0], _xyz[1], _xyz[2]};
      iopdet_idx = 0; 
      for (const auto& iopdet : opdets) {
        double& vis = _opdet_visd[iopdet]; 
        if (vis > vis_threshold) {
          ((THnSparseF*)h3opDet[iopdet_idx])->Fill(xyz, vis); 
        }
        iopdet_idx++; 
      }

      if (irun == 0) {
        hnTotal->Fill( xyz, *(_vis_d) ); 
      }
    }

    if (tVisBuf) {
      if (tVisBuf->GetEntries() > 0) {
        TTreeReader readerBuff(tVisBuf); 
        TTreeReaderArray<double> _buf_xyz(readerBuff, "coords"); 
        TTreeReaderArray<double> _buf_opdet_visd(readerBuff, "opDet_visDirectBuff"); 
        TTreeReaderArray<double> _buf_opdet_visr(readerBuff, "opDet_visReflctBuff"); 
        TTreeReaderValue<double> _buf_vis_d(readerBuff, "total_visDirectBuff"); 
        TTreeReaderValue<double> _buf_vis_r(readerBuff, "total_visReflctBuff"); 

        while( readerBuff.Next() ) {
          double xyz[3] = {_buf_xyz[0], _buf_xyz[1], _buf_xyz[2]};
          iopdet_idx = 0; 
          for (const auto& iopdet : opdets) {
            double& vis = _opdet_visd[iopdet]; 
            if ( vis > vis_threshold ) {
              ((THnSparseF*)h3opDet[iopdet_idx])->Fill(xyz, vis); 
            }
            iopdet_idx++; 
          }

          if (irun == 0) {
            hnTotal->Fill( xyz, *(_buf_vis_d) ); 
          }

        }
      }
    }

    output_h3->cd(); 
    if (irun == 0) {
      hnTotal->Write(); 
      delete hnTotal;
    }
    for (const auto& h3 : h3opDet) {
      h3->Write(); 
    }
    h3opDet.Clear();

    irun++; 
  }

  file->Close();

  return;
}

