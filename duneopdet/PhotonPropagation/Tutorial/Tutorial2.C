

void Tutorial2(int iPD = 0)
{

    // Load the hist file we produced
    TFile *fin = new TFile("dune1x2x6_optical_example_yscan_hist.root");

    // Load the trees
    TTree *pdtree = (TTree*)fin->Get("pmtresponse/OpDets");
    TTree *mctree = (TTree*)fin->Get("AnalysisExample/AnalysisExampleSimulation");

    // Attach variables to the branches we need
    int    EventID;        pdtree->SetBranchAddress("EventID",       &EventID);
    int    OpChannel;      pdtree->SetBranchAddress("OpChannel",     &OpChannel);
    int    CountDetected;  pdtree->SetBranchAddress("CountDetected", &CountDetected);
    double StartXYZT[4];   mctree->SetBranchAddress("StartXYZT",     StartXYZT);

    // Number of events and PDs.
    const int nevents = mctree->GetEntries();
    int numPDs_in_file = pdtree->GetEntries() / nevents;

    // Arrays to hold the y-values and the photon signals
    double axisarray[99];

    // Loop through the events to find the y positions of the muons
    for (int ievent = 0; ievent < nevents; ievent++) {
        mctree->GetEntry(ievent);
        axisarray[ievent] = StartXYZT[1]; // Store the y (1) start position of the muon
    }

    // Add the upper edge of the last bin to the axis
    axisarray[nevents] = 2*axisarray[nevents-1] - axisarray[nevents-2];
    
    // Creat the histogram
    TH1D *hYScan0 = new TH1D("hYScan0", TString::Format("Detector %i;y position of #mu (cm);A.U.",iPD), nevents, axisarray);

    int ientry = 0;

    // Loop through the events and fint eh PD signal which matches up
    for (int ievent = 0; ievent < nevents; ievent++) {
        mctree->GetEntry(ievent);

        // Find the right entry for event/pd
        while (ientry < pdtree->GetEntries()) {
          pdtree->GetEntry(ientry);
          
          
          if (EventID < ievent+1){
            ientry += 1;
            continue;
          }
          if (EventID > ievent+1){
            ientry += 1;
            break;
          }
          if (iPD != OpChannel){
            ientry += 1;
            continue;
          }
          //cout << " ientry " << ientry << " ievent " << ievent+1 << " EventID " << EventID << " iPD " << iPD << " OpChannel " << OpChannel << " " ;
          //cout << "Plot!" << endl;
          
          // Store the y (1) start position of the muon
          double ypos = StartXYZT[1];
          
          // Fill the histogram at this position with the number of counts detected
          hYScan0->Fill(ypos, CountDetected);

          ientry += 1;
          break;
        }
        
    }
    
    hYScan0->Draw("hist");
}
