class LightYieldMap {
  public:
    PhotonLibrary* lib;
    Int_t normDir; // Direction normal to the light yield plot
    Int_t normBoundsInVx[2];
    Int_t mapDirs[2];
    Double_t map[maxVxDim][maxVxDim];
    Double_t average;
    Double_t minimum;
    TString title;
    TString subtitle;

    LightYieldMap(PhotonLibrary* myLib, Int_t myNormDir, Int_t myNormBoundsInVx[2], TString descriptor) {
      lib = myLib;

      // Error handling
      if (!(0<= myNormDir && myNormDir<=2)) {
        cout << "ERROR: Normal direction not valid." << endl;
      } 
      if (!(0 <= myNormBoundsInVx[0] && myNormBoundsInVx[0] < lib->cryostatBoundsInVx[myNormDir][1])) {
        cout << "ERROR: Position is out of cryo range." << endl;
      } else if (!(lib->detectorBoundsInVx[myNormDir][0] <= myNormBoundsInVx[0] && myNormBoundsInVx[0] < lib->detectorBoundsInVx[myNormDir][1])) {
        cout << "WARNING: Position is out of fiducial range." << endl;
      }
      if (!(0 <= myNormBoundsInVx[1]-1 && myNormBoundsInVx[1]-1 < lib->cryostatBoundsInVx[myNormDir][1])) {
        cout << "ERROR: Position is out of cryo range." << endl;
      } else if (!(lib->detectorBoundsInVx[myNormDir][0] <= myNormBoundsInVx[1]-1 && myNormBoundsInVx[1]-1 < lib->detectorBoundsInVx[myNormDir][1])) {
        cout << "WARNING: Position is out of fiducial range." << endl;
      }

      normDir = myNormDir;
      for (int i=0; i<2; i++) normBoundsInVx[i] = myNormBoundsInVx[i];
      SetMapDirs();
      SetMap();
      SetAvg();
      SetMin();
      title = Form("Light Yield (%s%c = %.2f - %.2f m);%c [m];%c [m]; Light Yield [PE/MeV]"
          , descriptor.Data(), dirNames[normDir], lib->GetPos(normDir,Double_t(normBoundsInVx[0]))
          , lib->GetPos(normDir,Double_t(normBoundsInVx[1])), dirNames[mapDirs[1]], dirNames[mapDirs[0]]);
      subtitle = Form("LY_min = %.2f   <LY> = %.2f", minimum, average);
      return 0;
    }
    void SetMapDirs() {
      int i=0;
      for (int j=0; j<3; j++)
        if (j!=normDir)
          mapDirs[i++] = j;
    }
    void SetMap() {
      Int_t indices[3] = {0,0,0};
      for (indices[mapDirs[0]]=lib->detectorBoundsInVx[mapDirs[0]][0]; indices[mapDirs[0]]<lib->detectorBoundsInVx[mapDirs[0]][1]; indices[mapDirs[0]]++) {
        for (indices[mapDirs[1]]=lib->detectorBoundsInVx[mapDirs[1]][0]; indices[mapDirs[1]]<lib->detectorBoundsInVx[mapDirs[1]][1]; indices[mapDirs[1]]++) {
          map[indices[mapDirs[0]]][indices[mapDirs[1]]] = 0;
          for (indices[normDir]=normBoundsInVx[0]; indices[normDir]<normBoundsInVx[1]; indices[normDir]++) {
            map[indices[mapDirs[0]]][indices[mapDirs[1]]] += lib->LYPerVx[indices[0]][indices[1]][indices[2]];
          }
          map[indices[mapDirs[0]]][indices[mapDirs[1]]] /= (normBoundsInVx[1] - normBoundsInVx[0]);
        }
      }
    }
    void SetAvg() {
      Double_t numerator = 0;
      Double_t denominator =
        (lib->detectorBoundsInVx[mapDirs[0]][1] - lib->detectorBoundsInVx[mapDirs[0]][0])
        * (lib->detectorBoundsInVx[mapDirs[1]][1] - lib->detectorBoundsInVx[mapDirs[1]][0]);
      for (int i=lib->detectorBoundsInVx[mapDirs[0]][0]; i<lib->detectorBoundsInVx[mapDirs[0]][1]; i++) {
        for (int j=lib->detectorBoundsInVx[mapDirs[1]][0]; j<lib->detectorBoundsInVx[mapDirs[1]][1]; j++) {
          numerator += map[i][j];
        }
      }
      average = numerator/denominator;
    }
    void SetMin() {
      Double_t minVal = 1.79769e+308;
      Int_t minPos[2] = {-1,-1};
      for (int i=lib->detectorBoundsInVx[mapDirs[0]][0]; i<lib->detectorBoundsInVx[mapDirs[0]][1]; i++) {
        for (int j=lib->detectorBoundsInVx[mapDirs[1]][0]; j<lib->detectorBoundsInVx[mapDirs[1]][1]; j++) {
          if (map[i][j] < minVal) {
            minVal = map[i][j];
            minPos[0] = i;
            minPos[1] = j;
          }
        }
      }
      minimum = minVal;
    }
    void Draw() {
      TH2D *mapHist = new TH2D("mapHist",title.Data(),
          (lib->detectorBoundsInVx[mapDirs[1]][1] - lib->detectorBoundsInVx[mapDirs[1]][0]),
          lib->GetPos(mapDirs[1], lib->detectorBoundsInVx[mapDirs[1]][0]),
          lib->GetPos(mapDirs[1], lib->detectorBoundsInVx[mapDirs[1]][1]),
          (lib->detectorBoundsInVx[mapDirs[0]][1] - lib->detectorBoundsInVx[mapDirs[0]][0]),
          lib->GetPos(mapDirs[0], lib->detectorBoundsInVx[mapDirs[0]][0]),
          lib->GetPos(mapDirs[0], lib->detectorBoundsInVx[mapDirs[0]][1]));
      for (int i=lib->detectorBoundsInVx[mapDirs[0]][0]; i<lib->detectorBoundsInVx[mapDirs[0]][1]; i++) {
        Double_t x = lib->GetPos(mapDirs[0],Double_t(i)+0.5);
        for (int j=lib->detectorBoundsInVx[mapDirs[1]][0]; j<lib->detectorBoundsInVx[mapDirs[1]][1]; j++) {
          Double_t y = lib->GetPos(mapDirs[1],Double_t(j)+0.5);
          Int_t bin = mapHist->Fill(y,x,map[i][j]);
        }
      }

      TCanvas* myC = new TCanvas("myC", "Light Yield", 20, 52, 1250,
          Double_t(lib->detectorBoundsInVx[mapDirs[0]][1])
          /Double_t(lib->detectorBoundsInVx[mapDirs[1]][1])*1250.);
      myC->SetRightMargin(0.13);
      myC->SetTopMargin(0.13);
      mapHist->SetStats(0);
      gStyle->SetNumberContours(20);
      gStyle->SetPalette(kBird);
      mapHist->GetZaxis()->SetRangeUser(0, 110);
      mapHist->Draw("colz");
      //mapHist->GetZaxis()->SetRangeUser(0.4, 0.9);

      TPaveText *pt = new TPaveText(0.3,0.88,0.7,0.93,"blNDC");
      pt->SetName("subtitle");
      pt->SetBorderSize(2);
      pt->SetFillColor(19);
      pt->AddText(subtitle.Data());
      pt->Draw("same");
    }
};


