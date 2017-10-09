#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <fstream>

void AssembleSingleFile(std::string FileList, std::string BaseDirectory, std::string OutputName);
TChain * CreateChainFromList_opt(std::string ListFileName, std::string ChainName, std::string BaseDirectory, bool DoCheck=false);


void AssembleSingleFile(std::string FileList, std::string BaseDirectory, std::string OutputName)
{
  std::cout<<"Begining AssembleSingleFile.C\n";
  
  TChain * ch = CreateChainFromList_opt(FileList.c_str(),BaseDirectory.c_str(),"pmtresponse/PhotonLibraryData",false);
  
  Int_t Voxel, OpChannel;
  Float_t Visibility;
  ch->SetBranchAddress("Voxel",      &Voxel);
  ch->SetBranchAddress("OpChannel",  &OpChannel);
  ch->SetBranchAddress("Visibility", &Visibility);

  std::cout<<"Create Output File\n";
  TFile *f = new TFile(OutputName.c_str(),"RECREATE");
  TTree * tt = new TTree("PhotonLibraryData","PhotonLibraryData");
  tt->Branch("Voxel",      &Voxel,      "Voxel/I");
  tt->Branch("OpChannel",  &OpChannel,  "OpChannel/I");
  tt->Branch("Visibility", &Visibility, "Visibility/F");
  
  std::cout<<"Start Fill of output file branches\n";
  long int nEnt = ch->GetEntries();
  std::cout<<"Need to loop over "<<nEnt<<" entries.";
  for(long int i=0; i<nEnt; ++i)
    {
      if( i%1000==0) {
        std::cout<<"Now on entry "<<i<<" of "<<nEnt<<"\n";
      }else if( i>nEnt || i<0 ){
        std::cout<<"ERROR: i="<<i<<" does not fall in the range 0 to "<<nEnt<<". Aborting.\n";
//        throw std::logic_error();
        throw 20;
      }

      ch->GetEntry(i);
      tt->Fill();
    }
  
  std::cout<<"write output file\n";
  f->Write();
  f->Close();
}




TChain * CreateChainFromList_opt(std::string ListFileName, std::string BaseDirectory, std::string ChainName, bool DoCheck)
{
  ifstream InputFile(ListFileName.c_str());
  std::string FileName;

  TChain * TheChain = new TChain(ChainName.c_str());
  if(!DoCheck)
    {
      while(getline(InputFile, FileName))
        {
          FileName=BaseDirectory+FileName;
          std::cout<<"Adding file "<< FileName.c_str()<<" to the TChain"<<std::endl;
          TheChain->Add(FileName.c_str());
        }
    }
  else
    {
      while(getline(InputFile, FileName))
        {
          FileName=BaseDirectory+FileName;
          TFile*f=TFile::Open(FileName.c_str());
          if(f->Get(ChainName.c_str()))
            {
              TheChain->Add(FileName.c_str());
              std::cout<<"Adding file "<< FileName.c_str()<<" to the TChain"<<std::endl;
            }
          else std::cout<<"Chain " <<ChainName << " not found in file " << FileName<<std::endl;
        }


    }
  return TheChain;
}



