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


/**
 * @brief Rebin the visibility map averaging the bin content after rebin
 *
 * @param hn visibility map
 * @param rb0 bin group for axis 0
 * @param rb1 bin group for axis 1
 * @param rb2 bin group for axis 2
 * @return rebinned visibility map
 */
inline THnSparse* rebin_visibility_map(
    const THnSparse* hn, 
    const int rb0, 
    const int rb1, 
    const int rb2) 
{
  const Int_t rbidx[3] = {rb0, rb1, rb2}; 
  THnSparse* hnrb = hn->Rebin(rbidx); 
  
  float scale_factor = 1.0 / (rb0 * rb1 * rb2); 
  hnrb->Scale( scale_factor ); 
  return hnrb; 
}

/**
 * @brief Rebin the visibility map averaging the bin content after rebin
 *
 * @param hn visibility map
 * @param rb0 bin group for axis 0
 * @param rb1 bin group for axis 1
 * @param rb2 bin group for axis 2
 * @return rebinned visibility map
 */
inline THn* rebin_visibility_map(const THn* hn, const int rb0, const int rb1, const int rb2) {
  const Int_t rbidx[3] = {rb0, rb1, rb2}; 
  THn* hnrb = hn->Rebin(rbidx); 
  
  float scale_factor = 1.0 / (rb0 * rb1 * rb2); 
  hnrb->Scale( scale_factor ); 
  return hnrb; 
}

/**
 * @brief Get axis along which averaging the visibility (internal use)
 */
inline const int get_projection_axis(const int& axis0, const int& axis1) {
  assert(axis0 != axis1); 
  assert(axis1 < 2); 
  assert(axis0 < 2); 

  return (3-axis0-axis1);
}

/**
 * @brief Compute a 2D average visibility profile
 *
 * @param hn visibility map
 * @param axis0 axis index on the projection y-axis
 * @param axis1 axis index on the projection x-axis
 * @param x0 lower edge of the projected axis range
 * @param x1 upper edge of the projected axis range
 * @return projected visibility
 */
inline TH2D* get_projection(const THnBase* hn, 
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

#endif /* end of include guard DUNEVISUTILS_HPP */

