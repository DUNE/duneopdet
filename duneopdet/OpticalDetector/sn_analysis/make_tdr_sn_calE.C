/////////////////////////////////////////////////////////
// This script will help you on providing TDR plot for
// calorimetry study of SN neutrino.
// Note: Cut has been improved later
//
// for running this script
// root -l
// [] .L make_tdr_sn_calE.C
// [] make_tdr_sn_calE_all()  
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

void PrintPlots(TH2F*& h, bool profile, TString name)
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

TString title(int i, TString s){
  TString title = s;
  title+="_";
  title+=i;
  return title;
}

//-----------------------------------------------
TH1F* ResolnBias(TH2F*& h, bool bias, TString xyname)
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
void make_tdr_sn_calE(TString dummy = "flashmatchDEF45cm0100Hz5snrNonRefl/SelectedFlashTree", 
		      int nbinx = 24, int minx = 0, int maxx = 360, 
		      int nbiny = 30, int miny = 0, int maxy = 60,
		      int nspy = 5, int spminy = 0, int spmaxy = 10,
		      int nsx = 30, int sminx = 0, int smaxx = 60,
		      int nsy = 100, int sminy = 0, int smaxy = 100,
		      int nry = 100, int rminy = -1, int rmaxy = 1,
		      bool pe7 = true,  bool pe14 = true,  bool pe21 = true, bool pe28 = true,
		      bool opt = true,  bool pes = true)
{

  TChain *chain = new TChain(dummy);
  //  chain->Add("/pnfs/dune/scratch/users/bbehera/recox/mcc11_pd_tdr_dra*.root");
  //  for (const auto &filename : glob("/pnfs/dune/tape_backed/dunepro/mcc11/protodune/unknown/root-hist/09/99/46/*/mcc11_pd_tdr_dra_hist*.root")) {

  for (const auto &filename : glob("/pnfs/dune/tape_backed/dunepro/mcc11/protodune/unknown/root-hist/09/99/*/*/mcc11_pd_tdr_dra_hist*.root")) {
    chain->Add(filename.c_str());
  }


  // Attach variables to the branches we need
  float Purity;    chain->SetBranchAddress("Purity",      &Purity);
  float TotalPE;   chain->SetBranchAddress("TotalPE",     &TotalPE);
  float TrueE;     chain->SetBranchAddress("TrueE",       &TrueE);
  float TrueX;     chain->SetBranchAddress("TrueX",       &TrueX);
  //float RecoX;     chain->SetBranchAddress("RecoX",       &RecoX);
  float TrueY;     chain->SetBranchAddress("TrueY",       &TrueY);
  float TrueZ;     chain->SetBranchAddress("TrueZ",       &TrueZ);

  TH2F *hnue_x_pecut; 
  TH2F *hsmear_x_pecut; TH2F *hreso_x_pecut;
  TH1F *hpe;   TH1F *hpelt;
  
  hnue_x_pecut   = new TH2F("hnue_x_pecut",   TString::Format("; |True X| (cm); #frac{Total PE}{True E} (#frac{PE}{MeV})"), nbinx, minx, maxx, nbiny, miny, maxy);
  hsmear_x_pecut = new TH2F("hsmear_x_pecut", TString::Format(" ; True Neutrino Energy (MeV);  Reconstructed Neutrino Energy (MeV)"), nsx, sminx, smaxx, nsy, sminy, smaxy);
  hreso_x_pecut  = new TH2F("hreso_x_pecut",  TString::Format(" ; True Neutrino Energy (MeV) ; #frac{Reco - True}{True}"), nsx, sminx, smaxx, nry, rminy, rmaxy);
    
  hpe   = new TH1F("hpe",    TString::Format(" purity > 0.5; Total PE ; Events"), 100, 0, 1000);
  hpelt = new TH1F("hpelt",  TString::Format(" purity < 0.5; Total PE ; Events"), 100, 0, 1000);

  float recoe_x_pecut;

  int ientry = 0;
  Int_t nevent = chain->GetEntries();

  for (Int_t i=0; i < nevent;i++) {
    chain->GetEntry(i);              //read complete accepted event in memory
  
    // Avoid effects of walls  
    if (TrueX < 20 ) continue;
    if (TrueZ > 1000 || TrueZ <  300) continue;
    if (TrueY > 300  || TrueY < -300) continue;

    // Fill the signal and background event
    // based on purity distribution.
  
    if (Purity > 0.5) {
      hpe->Fill(Purity);
    }else{
      hpelt->Fill(Purity);
    }

    // PE cuts applied here.....
    //------------------------------------------------------------------

    // Light yield of 7 PE/MeV    
    if (pe7){

      // PE cut has applied to remove the background events.
      // for 7 PE/MeV, PE > 30 was decided based on looking
      // at signal(purity > 0.5) over background (purity < 0.5)

      if (TotalPE > 30){
	hnue_x_pecut->Fill(std::abs(TrueX), TotalPE/(TrueE*pow(10,3)));
	
	recoe_x_pecut = RecoE(std::abs(TrueX), TotalPE, 
			      16.48, -0.1653, 0.0008907, 
			      -2.596e-06, 2.986e-09);
	
	hsmear_x_pecut->Fill(TrueE*pow(10,3), recoe_x_pecut);
	
	hreso_x_pecut->Fill(TrueE*pow(10,3), (recoe_x_pecut - TrueE*(pow(10,3)))/(TrueE*pow(10,3)));
      }
    }
    
    // Light yield of 14 PE/MeV
    if (pe14){

      // PE cut has applied to remove the background events.
      // for 14 PE/MeV, PE > 40 was decided based on looking
      // at signal(purity > 0.5) over background (purity < 0.5)

      if (TotalPE > 40){
	hnue_x_pecut->Fill(std::abs(TrueX), TotalPE/(TrueE*pow(10,3)));
	
	recoe_x_pecut = RecoE(std::abs(TrueX), TotalPE,
			      34.49, -0.2769, 0.001023,
			      -1.986e-06, 1.601e-09);
	
	hsmear_x_pecut->Fill(TrueE*pow(10,3), recoe_x_pecut);
	
	hreso_x_pecut->Fill(TrueE*pow(10,3), (recoe_x_pecut - TrueE*(pow(10,3)))/(TrueE*pow(10,3)));
      }
    }

    // Light yield of 21 PE/MeV
    if (pe21){
      
      // PE cut has applied to remove the background events.
      // for 21 PE/MeV, PE > 80 was decided based on looking
      // at signal(purity > 0.5) over background (purity < 0.5)

      if ((TotalPE > 80)){
	hnue_x_pecut->Fill(std::abs(TrueX), TotalPE/(TrueE*pow(10,3)));
	
	recoe_x_pecut = RecoE(std::abs(TrueX), TotalPE, 
			      54.55, -0.4315, 0.001619,
			      -3.316e-06, 2.899e-09);
	
	hsmear_x_pecut->Fill(TrueE*pow(10,3), recoe_x_pecut);
	
	hreso_x_pecut->Fill(TrueE*pow(10,3), (recoe_x_pecut - TrueE*(pow(10,3)))/(TrueE*pow(10,3)));
      }
    }

    // Light yield of 28 PE/MeV
    if (pe28){

      // PE cut has applied to remove the background events.
      // for 28 PE/MeV, PE > 100 was decided based on looking
      // at signal(purity > 0.5) over background (purity < 0.5)

      if (TotalPE > 100){
	hnue_x_pecut->Fill(std::abs(TrueX), TotalPE/(TrueE*pow(10,3)));
	recoe_x_pecut = RecoE(std::abs(TrueX), TotalPE,
			      74.54, -0.5649, 0.001977, 
			      -3.733e-06, 3.006e-09);
	
	hsmear_x_pecut->Fill(TrueE*pow(10,3), recoe_x_pecut);
	
	hreso_x_pecut->Fill(TrueE*pow(10,3), (recoe_x_pecut - TrueE*(pow(10,3)))/(TrueE*pow(10,3)));
      }
    }

    // Optimistic reflection
    if (opt){

      // PE cut has applied to remove the background events.
      // PE > 100 was decided based on looking at
      // signal(purity > 0.5) over background (purity < 0.5)

      if (TotalPE > 100){
	
	hnue_x_pecut->Fill(std::abs(TrueX), TotalPE/(TrueE*pow(10,3)));
	recoe_x_pecut = RecoE(std::abs(TrueX), TotalPE,
			      54.96, -0.326, 0.0007629,
			      2.874e-07, -1.531e-09);
	
	hsmear_x_pecut->Fill(TrueE*pow(10,3), recoe_x_pecut);
	
	hreso_x_pecut->Fill(TrueE*pow(10,3), (recoe_x_pecut - TrueE*(pow(10,3)))/(TrueE*pow(10,3)));
      }
    }

    // Pessimistic reflection
    if (pes){

      // PE cut has applied to remove the background events.
      // PE > 50 was decided based on looking at 
      // signal(purity > 0.5) over background (purity < 0.5)

      if (TotalPE > 50){
	hnue_x_pecut->Fill(std::abs(TrueX), TotalPE/(TrueE*pow(10,3)));

	recoe_x_pecut = RecoE(std::abs(TrueX), TotalPE,
			      25.14, -0.1942, 0.0007413,
			      -1.328e-06, 9.895e-10);
		  
	hsmear_x_pecut->Fill(TrueE*pow(10,3), recoe_x_pecut);
	
	hreso_x_pecut->Fill(TrueE*pow(10,3), (recoe_x_pecut - TrueE*(pow(10,3)))/(TrueE*pow(10,3)));
      }
    }

  }// end of event loop
  
  
  TH1F* hresn_x_pecut;
  
  if(pe21){
    hresn_x_pecut     = ResolnBias(hreso_x_pecut,     false, "; True Neutrino Energy (MeV); Resolution");
    TFile fout("fout_resoln_sn_true_21pe.root", "RECREATE");  
    hresn_x_pecut->Write("hresn_x_pecut_21pe");
    hpe->Write("PEgt_21pe");
    hpelt->Write("PElt_21pe");
  }
  
  if(pe7){
    hresn_x_pecut     = ResolnBias(hreso_x_pecut,     false, "; True Neutrino Energy (MeV); Resolution");
    TFile fout("fout_resoln_sn_true_7pe.root", "RECREATE");
    hresn_x_pecut->Write("hresn_x_pecut_7pe");
    hpe->Write("PEgt_7pe");
    hpelt->Write("PElt_7pe");
  }
  
  if(pe14){
    hresn_x_pecut     = ResolnBias(hreso_x_pecut,     false, "; True Neutrino Energy (MeV); Resolution");
    TFile fout("fout_resoln_sn_true_14pe.root", "RECREATE");
    hresn_x_pecut->Write("hresn_x_pecut_14pe");
    hpe->Write("PEgt_21pe");
    hpelt->Write("PElt_21pe");
  }
  
  if(pe28){
    hresn_x_pecut     = ResolnBias(hreso_x_pecut,     false, "; True Neutrino Energy (MeV); Resolution");
    TFile fout("fout_resoln_sn_true_28pe.root", "RECREATE");
    hresn_x_pecut->Write("hresn_x_pecut_28pe");
    hpe->Write("PEgt_28pe");
    hpelt->Write("PElt_28pe");
  }
  
  if(opt){
    hresn_x_pecut     = ResolnBias(hreso_x_pecut,     false, "; True Neutrino Energy (MeV); Resolution");
    TFile fout("fout_resoln_sn_true_opt.root", "RECREATE");
    hresn_x_pecut->Write("hresn_x_pecut_opt");
    hpe->Write("PEgt_opt");
    hpelt->Write("PElt_opt");
  }

  if(pes){
    hresn_x_pecut     = ResolnBias(hreso_x_pecut,     false, "; True Neutrino Energy (MeV); Resolution");
    TFile fout("fout_resoln_sn_true_pes.root", "RECREATE");
    hresn_x_pecut->Write("hresn_x_pecut_pes");
    hpe->Write("PEgt_pes");
    hpelt->Write("PElt_pes");
  }
}

void make_tdr_sn_calE_all(){
  make_tdr_sn_calE("flashmatchEFF15cm0100Hz5snrNonRefl/SelectedFlashTree", 
		   24, 0, 360, 10, 0, 20, 10, 0, 10,
		   30,  0, 60, 100,  0,  100, 100,  -1,  1,
		   true, false, false, false, false, false);
  
  make_tdr_sn_calE("flashmatchEFF30cm0100Hz5snrNonRefl/SelectedFlashTree", 
		   24, 0, 360, 20, 0, 40, 10, 0, 10, 
		   30,  0, 60, 100,  0,  100, 100,  -1,  1,
		   false, true, false, false, false, false);
  
  make_tdr_sn_calE("flashmatchDEF45cm0100Hz5snrNonRefl/SelectedFlashTree", 
		   24,  0,  360, 30,  0,  60, 5,  0,  10,
		   30,  0, 60, 100,  0,  100, 100,  -1,  1,
		   false, false, true, false, false, false);
  
  make_tdr_sn_calE("flashmatchEFF60cm0100Hz5snrNonRefl/SelectedFlashTree", 
		   24, 0, 360, 40, 0, 80, 10, 0, 20,
		   30,  0, 60, 100,  0,  100, 100,  -1,  1,
		   false, false, false, true, false, false);
  
  make_tdr_sn_calE("flashmatchREF45cm0100Hz5snrOptRefl/SelectedFlashTree",
		   24,  0,  360, 30,  0,  60, 5,  0,  10,
                   30,  0, 60, 100,  0,  100, 100,  -1,  1,
		   false, false, false, false, true, false);
  
  make_tdr_sn_calE("flashmatchREF45cm0100Hz5snrPesRefl/SelectedFlashTree", 
		   24, 0, 360, 17, 0, 34, 8, 0, 15,
		   30,  0, 60, 100,  0,  100, 100,  -1,  1,
		   false, false, false, false, false, true);
  
  
} // end of void
