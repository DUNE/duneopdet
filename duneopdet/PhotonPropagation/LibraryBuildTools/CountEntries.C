#include "TTree.h"
#include "TFile.h"
#include "TDirectoryFile.h"

long int CountEntries(const char* inFile){
  TFile* f = TFile::Open(inFile);
  long int out = ((TTree*)(((TDirectoryFile*)(f->Get("pmtresponse")))->Get("PhotonLibraryData")))->GetEntries();
  return out;
}
