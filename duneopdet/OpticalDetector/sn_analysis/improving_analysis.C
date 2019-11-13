/////////////////////////////////////////////////////////
// This script is in under development
// and analysis is improving 
// In this script we have removed the PE Cut based on 
// signal over background events.
// 
// To run this script:
// 1. root -l
// 2. [] .L improving_analysis.C
// 3. all_configuration()
////////////////////////////////////////////////////////

#include <fstream>
#include <vector>
#include <glob.h>

using namespace std;

float RecoE(float trueX, float pe, float  p0, float  p1, 
	    float  p2, float  p3, float  p4){
  
  float value = 0.0;
  value = p0
    + p1 * trueX
    + p2 * std::pow(trueX, 2)
    + p3 * std::pow(trueX, 3)
    + p4 * std::pow(trueX, 4);
  
  float nuE = 0.0;
  nuE = pe/value;
  return nuE;
}

void PrintPlots(TH2F*& h, bool profile,
                TString name)
{
  TCanvas *c = new TCanvas();
  gStyle->SetPalette(1);
  c->SetLeftMargin(0.15);
  c->SetRightMargin(0.15);
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  if(profile){
    TProfile* phisto1 = h->ProfileX("i");
    h->Draw("colz");
    phisto1->Draw("ep same");
    c->SaveAs( "plots/" + name);
  }else{
    h->Draw("colz");
    c->SaveAs( "plots/" + name);
  }
}

void CenterTitles(TH2F*& h, TCanvas *c)
{
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  h->GetZaxis()->CenterTitle();
  c->SetLeftMargin(0.15);
  c->SetRightMargin(0.15);
  c->SetBottomMargin(0.15);
}

void CenterTitles(TH1F*& h, TCanvas *c)
{
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  h->GetZaxis()->CenterTitle();
  c->SetLeftMargin(0.15);
  c->SetRightMargin(0.15);
  c->SetBottomMargin(0.15);
}
//----------------------------------
TString title(int i, TString s){
  TString title = s;
  title+="_";
  title+=i;
  return title;
}

//-----------------------------------------------

TH1F* ResolnBias(TH2F*& h, bool bias,
		 TString xyname)
{
  int ibinx =  h->GetNbinsX();
  TH1F* h1   = new TH1F("h1"," " + xyname + " ",30, 0, 60);
  TString title(int i, TString s);
  TH1F* hold[ibinx];
  double rms[ibinx];
  double rmserror[ibinx];
  double mean[ibinx];
  double meanerror[ibinx];
  
  for (int ibiny = 0; ibiny < ibinx; ibiny++){
    hold[ibiny]=(TH1F*)h->ProjectionY(title(ibiny,"h1pol"),ibiny+1,ibiny+1);
    
    if (bias){
      mean[ibiny]      = hold[ibiny]->GetMean();
      meanerror[ibiny] = hold[ibiny]->GetMeanError();
    }else{
      rms[ibiny]      = hold[ibiny]->GetRMS();
      rmserror[ibiny] = hold[ibiny]->GetRMSError();
    }
  }
  
  
  for (int j = 0; j< ibinx; j++){
    if (bias){
      h1->SetBinContent(j+1, mean[j]);
      h1->SetBinError(j+1, meanerror[j]);
    }else{
      h1->SetBinContent(j+1,rms[j]);
      h1->SetBinError(j+1,rmserror[j]);
    }
  }
  return h1;
}

std::vector<std::string> glob(const char *pattern) {
  glob_t g;
  glob(pattern, GLOB_TILDE, nullptr, &g); // one should ensure glob returns 0!
  std::vector<std::string> filelist;
  filelist.reserve(g.gl_pathc);
  for (size_t i = 0; i < g.gl_pathc; ++i) {
    filelist.emplace_back(g.gl_pathv[i]);
  }
  globfree(&g);
  return filelist;
}

//---------------------------------------------------------------------------------------
void improving_analysis(TString dummy = "flashmatchDEF45cm0100Hz5snrNonRefl/SelectedFlashTree", 
			int nbinx = 24, int minx = 0, int maxx = 360, 
			int nbiny = 30, int miny = 0, int maxy = 60,
			int nspy = 5, int spminy = 0, int spmaxy = 10,
			int nsx = 30, int sminx = 0, int smaxx = 60,
			int nsy = 100, int sminy = 0, int smaxy = 100,
			int nry = 100, int rminy = -1, int rmaxy = 1,
			bool eff35 = true,  bool eff15 = true,  bool eff25 = true, bool eff45 = true,
			bool opt = true,  bool pes = true, bool pe2 = true, bool pe3 = true)
{
  TFile *f = new TFile("/dune/data/users/bbehera/2019/Marley/marley_rad_ana_hist.root");
  TTree * chain = (TTree*) f->Get(dummy);
  
  // Attach variables to the branches we need
  float Purity;     chain->SetBranchAddress("Purity",      &Purity);
  float TotalPE;    chain->SetBranchAddress("TotalPE",     &TotalPE);
  float TrueE;      chain->SetBranchAddress("TrueE",       &TrueE);
  float TrueX;      chain->SetBranchAddress("TrueX",       &TrueX);
  float RecoX;      chain->SetBranchAddress("RecoX",       &RecoX);
  float TrueY;      chain->SetBranchAddress("TrueY",       &TrueY);
  float TrueZ;      chain->SetBranchAddress("TrueZ",       &TrueZ);
  float YWidth;     chain->SetBranchAddress("YWidth",       &YWidth);
  float ZWidth;     chain->SetBranchAddress("ZWidth",       &ZWidth);
  float ZCenter;    chain->SetBranchAddress("ZCenter",       &ZCenter);
  float YCenter;    chain->SetBranchAddress("YCenter",       &YCenter);
  int   NHitOpDets; chain->SetBranchAddress("NHitOpDets",    &NHitOpDets);
  
  TH2F *hnue_x_pecut;
  TH2F *hsmear_x_pecut; TH2F *hreso_x_pecut;
  TH1F *hpe;   TH1F *hpelt; 
  
  hnue_x_pecut     = new TH2F("hnue_x_pecut",     TString::Format(" ; |True X| (cm); #frac{Total PE}{True E} (#frac{PE}{MeV})"), nbinx, minx, maxx, nbiny, miny, maxy);
  hsmear_x_pecut = new TH2F("hsmear_x_pecut", TString::Format(" ; True Neutrino Energy (MeV);  Reconstructed Neutrino Energy (MeV)"), nsx, sminx, smaxx, nsy, sminy, smaxy);
  hreso_x_pecut  = new TH2F("hreso_x_pecut",  TString::Format(" ; True Neutrino Energy (MeV) ; #frac{Reco - True}{True}"), nsx, sminx, smaxx, nry, rminy, rmaxy);

  hpe   = new TH1F("hpe",  TString::Format(" Purity > 0.5; Total PE ; Events"), 100, 0, 1500);
  hpelt = new TH1F("hpelt",  TString::Format("Purity < 0.5; Total PE ; Events"), 100, 0, 1500);

  

  int ientry = 0;
  Int_t nevent = chain->GetEntries();

  float recoe_x_pecut;  
  
  for (Int_t i=0; i < nevent;i++) {
    chain->GetEntry(i);              //read complete accepted event in memory
    
    entry1++;  
    float deltaz = ZCenter-TrueZ;
    float deltay = YCenter-TrueY;

    // Avoid effects of walls
    if (TrueX < 20 ) continue;
    if (TrueY > 300  || TrueY < -300) continue;

    // new set of cuts for removing the almost all wrong flashes event.    
    if (NHitOpDets < 10) continue; 
    if (deltaz < -100 || deltaz > 100) continue;
    if (deltay < -100 || deltay > 100) continue;
    if (YWidth < 60) continue; 
    if (ZWidth < 60) continue;


    if (Purity > 0.5) {
      hpe->Fill(TotalPE);
    }else{
      hpelt->Fill(TotalPE);
    }

    // Filling TrueX vs TotalPE here.....
    //------------------------------------------------------------------
    if (pe2){
      
      hnue_x_pecut->Fill(std::abs(TrueX), TotalPE/(TrueE*pow(10,3)));
      recoe_x_pecut = RecoE(std::abs(TrueX), TotalPE,
			    69.49,  -0.5379,  0.00269, 
			    -9.027e-06,  1.26e-08);
        
      hsmear_x_pecut->Fill(TrueE*pow(10,3), recoe_x_pecut);
      
      hreso_x_pecut->Fill(TrueE*pow(10,3), (recoe_x_pecut - TrueE*(pow(10,3)))/(TrueE*pow(10,3)));
    }

  }// end of event loop

  //------------------------------------------------------------------
  TH1F* hresn_x_pecut;

  
  
  if(pe2){
    hresn_x_pecut     = ResolnBias(hreso_x_pecut,     false, "; True Neutrino Energy (MeV); Resolution");
    TFile fout("fout_resoln_sn_true_pe2_fixed.root", "RECREATE");
    TProfile* phisto_x_pecut = hnue_x_pecut->ProfileX("i");
    phisto_x_pecut->Fit("pol4");
    hnue_x_pecut->Draw("colz");
    phisto_x_pecut->Draw("ep same");
    hnue_x_pecut->Write("truex_vs_pe");
    phisto_x_pecut->Write("profile_truex_vs_pe");
    hsmear_x_pecut->Write("smear_x_pecut");
    hreso_x_pecut->Write("reso_x_pecut");
    hresn_x_pecut->Write("hresn_x_pecut");
    hpe->Write("PEgt");
    hpelt->Write("PElt");
  }
  
}// end of void

void all_configuration(){
  
  improving_analysis("flashmatchTHR35QENonRefl2PE/SelectedFlashTree",
		     24,  0,  360, 30,  0,  60, 5,  0,  10,
		     30,  0, 60, 100,  0,  100, 100,  -1,  1,
		     false, false, false, false, false, false, true, false);
} // end of void
