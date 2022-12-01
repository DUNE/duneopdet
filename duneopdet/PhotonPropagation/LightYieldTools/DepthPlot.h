// Tyler LaBree
// Northern Illinois University

/* Draw average and minimum light yield vs detector depth on a given photon
 * library and on the range defined.
 */
void DrawDepthPlot(PhotonLibrary* lib, Int_t minDepth, Int_t maxDepth) {
  Int_t dist = maxDepth-minDepth;
  Double_t depths[dist];
  Double_t averageLYs[dist];
  Double_t minimumLYs[dist];
  for (int i=0; i<dist; i++) {
    depths[i] = lib->GetPosInM(2,Double_t(i+minDepth)+0.5);
  }
  for (int i=0; i<dist; i++) {
    Int_t slice[2] = {i+minDepth,i+minDepth+1};
    LightYieldMap* depthLYMap = new LightYieldMap(lib, 2, slice, "");
    averageLYs[i] = depthLYMap->average;
    minimumLYs[i] = depthLYMap->minimum;
  }

  TCanvas* myC2 = new TCanvas("myC2", "Light Yield vs Depth", 20, 52, 1250, 650);
  TGraph *lyAvgDepth = new TGraph(dist,depths, averageLYs);
  TGraph *lyMinDepth = new TGraph(dist,depths, minimumLYs);
  lyAvgDepth->SetMarkerStyle(20);
  lyAvgDepth->SetMarkerColor(kBlue);
  lyMinDepth->SetMarkerStyle(20);
  lyMinDepth->SetMarkerColor(kRed);
  lyAvgDepth->SetTitle("Avg and Min Light Yield vs Detector Depth;z [m];LY [PE/MeV]");
  lyAvgDepth->Draw("AP");
  lyAvgDepth->GetYaxis()->SetRangeUser(0.,32.);
  lyMinDepth->Draw("SAMEP");
}
