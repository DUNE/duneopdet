/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : DUNEVisUtils.hpp
 * @created     : Friday Oct 18, 2024 17:19:58 CEST
 */

#ifndef DUNEVISUTILS_HPP

#define DUNEVISUTILS_HPP

#include <TFile.h>
#include <THn.h>
#include <THnSparse.h>
#include <TH3D.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TProfile2D.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TROOT.h>

inline TH3* rebin_visibility_map(const TH3* h3, const int rbx, const int rby, const int rbz)
{
  TH3* h3rb = (TH3*)h3rb->Rebin3D(rbx, rby, rbz, Form("%s_rb_%i_%i_%i", h3->GetName(), rbx, rby, rbz)); 

  float scale_factor = 1.0 / (rbx * rby *rbz); 
  h3rb->Scale( scale_factor ); 
  return h3rb; 
}

inline THnSparse* rebin_visibility_map(const THnSparse* hn, const int rbx, const int rby, const int rbz) 
{
  const Int_t rbidx[3] = {rbx, rby, rbz}; 
  THnSparse* hnrb = hn->Rebin(rbidx); 
  
  float scale_factor = 1.0 / (rbx * rby * rbz); 
  hnrb->Scale( scale_factor ); 
  return hnrb; 
}

inline THn* rebin_visibility_map(const THn* hn, const int rbx, const int rby, const int rbz) {
  const Int_t rbidx[3] = {rbx, rby, rbz}; 
  THn* hnrb = hn->Rebin(rbidx); 
  
  float scale_factor = 1.0 / (rbx * rby * rbz); 
  hnrb->Scale( scale_factor ); 
  return hnrb; 
}

inline const int get_projection_axis(const int& axis0, const int& axis1) {
  assert(axis0 != axis1); 
  assert(axis1 < 2); 
  assert(axis0 < 2); 

  return (3-axis0-axis1);
}

inline TH2D* draw_projection(const THnBase* hn, 
    const int& axis0,  
    const int& axis1, 
    const double& x0, 
    const double& x1)
{
  const int proj_axis = get_projection_axis(axis0, axis1);
  const int ibin0 = hn->GetAxis(proj_axis)->FindBin(x0); 
  const int ibin1 = hn->GetAxis(proj_axis)->FindBin(x1); 
  const int proj_nbin = ibin1 - ibin0 + 1; 
  hn->GetAxis(proj_axis)->SetRange(ibin0, ibin1); 

  TH2D* h_proj = hn->Projection(axis0, axis1); 
  h_proj->Scale( 1.0 / proj_nbin ); 

  hn->GetAxis(proj_axis)->SetRange(); 
  return h_proj;
}

inline void CanvasPartition(TCanvas *C,const Int_t Nx,const Int_t Ny,
    Float_t lMargin, Float_t rMargin,
    Float_t bMargin, Float_t tMargin)
{
  if (!C) return;
  // Setup Pad layout:
  Float_t vSpacing = 0.0;
  Float_t vStep  = (1.- bMargin - tMargin - (Ny-1) * vSpacing) / Ny;
  Float_t hSpacing = 0.0;
  Float_t hStep  = (1.- lMargin - rMargin - (Nx-1) * hSpacing) / Nx;
  Float_t vposd,vposu,vmard,vmaru,vfactor;
  Float_t hposl,hposr,hmarl,hmarr,hfactor;
  for (Int_t i=0;i<Nx;i++) {
    if (i==0) {
      hposl = 0.0;
      hposr = lMargin + hStep;
      hfactor = hposr-hposl;
      hmarl = lMargin / hfactor;
      hmarr = 0.0;
    } else if (i == Nx-1) {
      hposl = hposr + hSpacing;
      hposr = hposl + hStep + rMargin;
      hfactor = hposr-hposl;
      hmarl = 0.0;
      hmarr = rMargin / (hposr-hposl);
    } else {
      hposl = hposr + hSpacing;
      hposr = hposl + hStep;
      hfactor = hposr-hposl;
      hmarl = 0.0;
      hmarr = 0.0;
    }
    for (Int_t j=0;j<Ny;j++) {
      if (j==0) {
        vposd = 0.0;
        vposu = bMargin + vStep;
        vfactor = vposu-vposd;
        vmard = bMargin / vfactor;
        vmaru = 0.0;
      } else if (j == Ny-1) {
        vposd = vposu + vSpacing;
        vposu = vposd + vStep + tMargin;
        vfactor = vposu-vposd;
        vmard = 0.0;
        vmaru = tMargin / (vposu-vposd);
      } else {
        vposd = vposu + vSpacing;
        vposu = vposd + vStep;
        vfactor = vposu-vposd;
        vmard = 0.0;
        vmaru = 0.0;
      }
      C->cd(0);
      char name[16];
      sprintf(name,"pad_%i_%i",i,j);
      TPad *pad = (TPad*) gROOT->FindObject(name);
      if (pad) delete pad;
      pad = new TPad(name,"",hposl,vposd,hposr,vposu);
      pad->SetLeftMargin(hmarl);
      pad->SetRightMargin(hmarr);
      pad->SetBottomMargin(vmard);
      pad->SetTopMargin(vmaru);
      pad->SetFrameBorderMode(0);
      pad->SetBorderMode(0);
      pad->SetBorderSize(0);
      pad->Draw();
    }
  }
}


#endif /* end of include guard DUNEVISUTILS_HPP */

