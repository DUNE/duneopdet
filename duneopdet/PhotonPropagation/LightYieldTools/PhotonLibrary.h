#include "TH2D.h"
#include "TH1.h"
// Detector dimensions
Double_t cryostatBoundsInM[3][2] = {
  {-4.2512,  4.252 },
  {-7.7412,  7.7412},
  {-1.0012, 21.9173}
};
Double_t detectorBoundsInM[3][2] = {
  {-3.25  ,  3.25  },
  {-6.75  ,  6.75  },
  { 0.    , 21.    }
};

// Extra constants
Char_t dirNames[3] = {'x', 'y', 'z'};
Double_t nPhotPerEDep = 25000.; // photons per MeV
Double_t arapucaEfficiency = 0.03;
const int maxVxDim = 200; // This reserves memory for the voxel data

class PhotonLibrary {
  public:
    Int_t cryostatBoundsInVx[3][2]; // Include min, but exclude max
    Int_t detectorBoundsInVx[3][2]; // Include min, but exclude max
    Double_t LYPerVx[maxVxDim][maxVxDim][maxVxDim];

    PhotonLibrary(TString myFilename, Int_t myVxDims[3]) {
      SetCryostatVol(myVxDims);
      SetDetectorVol();
      SetLYPerVx(myFilename);
    }

    void SetCryostatVol(Int_t vxDims[3]) {
      for (int i=0; i<3; i++) {
        cryostatBoundsInVx[i][0] = 0;
        cryostatBoundsInVx[i][1] = vxDims[i];
      }
    }
    void SetDetectorVol() {
      for (int i=0; i<3; i++) {
        detectorBoundsInVx[i][0] = Int_t(GetVox(i, detectorBoundsInM[i][0])+1);
        detectorBoundsInVx[i][1] = Int_t(GetVox(i, detectorBoundsInM[i][1]));
      }
    }
    void SetLYPerVx(TString filename) {
      Int_t nVx = 1;
      for (int i=0; i<3; i++)
        nVx *= (cryostatBoundsInVx[i][1] - cryostatBoundsInVx[i][0]);
      Double_t visPerVx[nVx];
      for (int i=0; i<nVx; i++)
        visPerVx[nVx] = 0;

      TFile *file = TFile::Open(filename);
      TTree *tree = (TTree*)file->Get("PhotonLibraryData");
      Int_t nEntries = tree->GetEntries();
      Int_t voxel, opChannel;
      Float_t visibility;
      tree->SetBranchAddress("Voxel",&voxel);
      tree->SetBranchAddress("OpChannel",&opChannel);
      tree->SetBranchAddress("Visibility",&visibility);
      for (Int_t i=0; i<nEntries; i++) {
        tree->GetEntry(i);
        //if (opChannel>=56) // cathode plane
        //if (opChannel<56) // field cage walls
        //if (opChannel<56 && opChannel%2==0) // right field cage wall
        //if (opChannel<56 && opChannel%2==1) // left field cage wall
        visPerVx[voxel] += visibility;
      }
      for (int i=0; i<cryostatBoundsInVx[0][1]; i++) {
        for (int j=0; j<cryostatBoundsInVx[1][1]; j++) {
          for (int k=0; k<cryostatBoundsInVx[2][1]; k++) {
            int vxIndex = i + j*cryostatBoundsInVx[0][1]
              + k*cryostatBoundsInVx[1][1]*cryostatBoundsInVx[0][1];
            LYPerVx[i][j][k] = GetLightYield(visPerVx[vxIndex]);
          }
        }
      }
    }
    Double_t GetVox(Int_t dir, Double_t pos) {
      return (pos - cryostatBoundsInM[dir][0])
        * Double_t(cryostatBoundsInVx[dir][1] - cryostatBoundsInVx[dir][0])
        / (cryostatBoundsInM[dir][1] - cryostatBoundsInM[dir][0]);
    }
    Double_t GetPos(Int_t dir, Double_t vox) {
      return vox * (cryostatBoundsInM[dir][1] - cryostatBoundsInM[dir][0])
        / Double_t(cryostatBoundsInVx[dir][1] - cryostatBoundsInVx[dir][0]) 
        + cryostatBoundsInM[dir][0];
    }
    Double_t GetLightYield(Double_t vis) {
      return vis * arapucaEfficiency * nPhotPerEDep; // PE/MeV
    }
};
