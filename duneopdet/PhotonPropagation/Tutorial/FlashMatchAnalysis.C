#define FlashMatchAnalysis_cxx
#include "FlashMatchAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

//This loop was created from the TTree::MakeClass() function and thus relies on the corresponding FlashMatchAnalysis.h file.

void FlashMatchAnalysis::Loop()
{
   if (fChain == 0) return;

   //Get each event
   Long64_t nentries = fChain->GetEntriesFast();

   //Initialize relevant histograms
   //In this macro we will be comparing specturms and efficiencies of the reconstruction while cheating to
   // select the flash we know is the supernova via mctruth to the spectrums and efficienies of adding various
   // cuts to come as close as possible to reproduce the former effeciency whithout mctruth

   //MCTruth Histograms
   TH1F *htrueE = new TH1F("mctruth E","energy specturm",50,0,0.06); //mctruth E
   TH1F *htrueX = new TH1F("mctruth x","position spectrum",50,-400,400); //mctruth X

   //Histograms related to the reconstructed flashes matched using mctruth
   TH1F *hpossibleselE = new TH1F("SN interactions which create flashes","SN interactions which create flashes",50,0,0.06);
   TH1F *hpossibleselX = new TH1F("SN interactions which create flashes","SN interactions which create flashes",50,-400,400);
   TH1F *hpossibleeffE = new TH1F("peff","energy efficiency",50,0,0.06);
   TH1F *hpossibleeffX = new TH1F("SN interactions which create flashes","position efficiency",50,-400,400);

   //Histograms related to selection and efficiency of a simple cut that looks for the largest pe flash in an event
   TH1F *hlargeflashselE = new TH1F("Largest flash selection","largest flash selection",50,0,0.06);
   TH1F *hlargeflashselX = new TH1F("Largest flash selection","largest flash selection",50,-400,400);
   TH1F *hlargeflasheffE = new TH1F("lfeff","largest flash selection",50,0,0.06);
   TH1F *hlargeflasheffX = new TH1F("Largest flash selection","largest flash selection",50,-400,400);

   //Histograms related to selection and efficiency of looking only in a 240cm radius around the supernova event and
   //in the yz plane and then finds the largest flash
   TH1F *hcutflashselE = new TH1F("Largest flash selection","cutst flash selection",50,0,0.06);
   TH1F *hcutflashselX = new TH1F("Largest flash selection","cutst flash selection",50,-400,400);
   TH1F *hcutflasheffE = new TH1F("ceff","cutst flash selection",50,0,0.06);
   TH1F *hcutflasheffX = new TH1F("Largest flash selection","cutst flash selection",50,-400,400);

   //Signal vs noise histograms. These are illsutratuve of why we use the cuts we use
   TH1F *hTotalPEVectorSN = new TH1F("TotalPESN","TotalPE SN",50,0,140);
   TH1F *hTotalPEVectorRAD = new TH1F("TotalPERAD","PE per flash",50,0,140);
   TH1F *hdistanceSN = new TH1F("distanceSN","distance SN",50,0,2000);
   TH1F *hdistanceRAD = new TH1F("distanceRAD","distance b/t flash and mctruth",50,0,2000);

//These two variables hold the selected flash for each event subject to the different cuts
//selflash0 corresponds to the largest pe cut while selflash1 corresponds to the largest pe + distance cut
int selflash0 = 0; 
int selflash1 = 0; 

//We'll use these to find the largest pe while looping over the flashes
float largestpe0 = 0;
float largestpe1 = 0;

//this variable is the distance between the reconstructed flash position and the mctrue position on the yz plane
//Though its true that we are cheating a little bit by using the mctrue position, its an ok approximation because
//when we are actually taking data we will have tpc information about y and z which should be quite accurate.
float distance;

//loop over every event
   for (Long64_t jentry=0; jentry<nentries;jentry++) { 
      fChain->GetEntry(jentry);
      int nflashes = (*Signal).size();
  
      htrueE->Fill(TrueE);
      htrueX->Fill(TrueX);

//loop over every flash
      for(int iflash=0; iflash<nflashes; iflash++) {

         distance = sqrt( (TrueY - (*YCenterVector)[iflash]) * (TrueY - (*YCenterVector)[iflash]) +  (TrueZ - (*ZCenterVector)[iflash]) * (TrueZ - (*ZCenterVector)[iflash]) );

//see if the supernovas in this flash has any corresponding reconstructed flashes   
         if( (*Signal)[iflash] ){
            hpossibleselE->Fill(TrueE);
            hpossibleselX->Fill(TrueX);
            hpossibleeffE->Fill(TrueE);
            hpossibleeffX->Fill(TrueX);

            hTotalPEVectorSN->Fill((*TotalPEVector)[iflash]);
            hdistanceSN->Fill(distance);
         }
         else {
            hTotalPEVectorRAD->Fill((*TotalPEVector)[iflash]);
            hdistanceRAD->Fill(distance);
           
         }

//Find the largest pe in the event.
         if( (*TotalPEVector)[iflash] > largestpe0 ) {
            largestpe0 = (*TotalPEVector)[iflash];
            selflash0 = iflash;
         }

//Find the largest pe within 240 cm of the true yz on the yz plane. 
         if( (*TotalPEVector)[iflash] > largestpe1 && distance < 240) {
            largestpe1 = (*TotalPEVector)[iflash];
            selflash1 = iflash;
         }
 

      }// done looping over flashes

//Fill the cut histograms for each event.      
      if( (*Signal)[selflash0] ) {
         hlargeflashselE->Fill(TrueE);
         hlargeflashselX->Fill(TrueX);
         hlargeflasheffE->Fill(TrueE);
         hlargeflasheffX->Fill(TrueX);

      }

      if( (*Signal)[selflash1] ) {
         hcutflashselE->Fill(TrueE);
         hcutflashselX->Fill(TrueX);
         hcutflasheffE->Fill(TrueE);
         hcutflasheffX->Fill(TrueX);

      }


      selflash0 = 0;
      selflash1 = 0;
      largestpe0 = 0;
      largestpe1 = 0;
   }// done looping over events

//The efficiency is the number of correct events divided by the number of events.
   hpossibleeffE->Divide(htrueE);
   hlargeflasheffE->Divide(htrueE);
   hpossibleeffX->Divide(htrueX);
   hlargeflasheffX->Divide(htrueX);

   hcutflasheffE->Divide(htrueE);
   hcutflasheffX->Divide(htrueX);

//The remainder of the code creates the pretty plots.
   hpossibleselE->SetLineColor(kGreen+2); 
   hlargeflashselE->SetLineColor(kRed);
   hpossibleselX->SetLineColor(kGreen+2); 
   hlargeflashselX->SetLineColor(kRed);

   hpossibleeffE->SetLineColor(kGreen+2); 
   hlargeflasheffE->SetLineColor(kRed);
   hpossibleeffX->SetLineColor(kGreen+2); 
   hlargeflasheffX->SetLineColor(kRed);

   hcutflashselE->SetLineColor(kCyan);
   hcutflashselX->SetLineColor(kCyan);

   hcutflasheffE->SetLineColor(kCyan);
   hcutflasheffX->SetLineColor(kCyan);

   TCanvas *c1 = new TCanvas("c1","c1");
   c1->cd();
   htrueE->GetXaxis()->SetTitle("energy (GeV)");
   htrueE->GetYaxis()->SetTitle("entries");
   htrueE->Draw();
   hpossibleselE->Draw("same");
   hlargeflashselE->Draw("same");
   hcutflashselE->Draw("same");


   TCanvas *c2 = new TCanvas("c2","c2");
   c2->cd();
   htrueX->Draw();
   htrueX->GetXaxis()->SetTitle("x position (cm)");
   htrueX->GetYaxis()->SetTitle("entries");
   hpossibleselX->Draw("same");
   hlargeflashselX->Draw("same");
   hcutflashselX->Draw("same");

   TCanvas *c3 = new TCanvas("c3","c3");
   c3->cd();
   hpossibleeffE->GetXaxis()->SetTitle("energy (GeV)");
   hpossibleeffE->GetYaxis()->SetTitle("efficiency");
   hpossibleeffE->Draw();
   hlargeflasheffE->Draw("same");
   hcutflasheffE->Draw("same");

   TCanvas *c4 = new TCanvas("c4","c4");
   c4->cd();
   hpossibleeffX->Draw();
   hpossibleeffX->GetYaxis()->SetTitle("efficiency");
   hpossibleeffX->GetXaxis()->SetTitle("x position (cm)");
   hlargeflasheffX->Draw("same");
   hcutflasheffX->Draw("same");

   TCanvas *c5 = new TCanvas("c5","c5");
   c5->cd();
   hTotalPEVectorRAD->GetXaxis()->SetTitle("total PE in flash");
   hTotalPEVectorRAD->GetYaxis()->SetTitle("entries");
   hTotalPEVectorRAD->Draw();
   hTotalPEVectorSN->SetLineColor(kRed);
   hTotalPEVectorSN->Draw("same");

   TCanvas *c6 = new TCanvas("c6","c6");
   c6->cd();
   hdistanceRAD->GetXaxis()->SetTitle("distance (cm)");
   hdistanceRAD->GetYaxis()->SetTitle("entries");
   hdistanceRAD->Draw();
   hdistanceSN->SetLineColor(kRed);
   hdistanceSN->Draw("same");


}
