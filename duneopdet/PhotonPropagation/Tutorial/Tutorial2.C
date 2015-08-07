

void Tutorial2()
{
    // Load the hist file we produced
    TFile *fin = new TFile("dune35t_optical_example_yscan_hist.root");

    // Load the trees
    TTree *pdtree = fin->Get("pmtresponse/OpDets");
    TTree *mctree = fin->Get("AnalysisExample/AnalysisExampleSimulation");

    // Attach variables to the branches we need
    int    EventID;        pdtree->SetBranchAddress("EventID",       &EventID);
    int    CountDetected;  pdtree->SetBranchAddress("CountDetected", &CountDetected);
    double StartXYZT[4];   mctree->SetBranchAddress("StartXYZT",     StartXYZT);

    // Number of events and PDs.
    const int nevents = mctree->GetEntries();
    int numPDs_in_file = pdtree->GetEntries() / nevents;

    // Let's look at the Photon Detector 0
    int iPD = 0;

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
    TH1D *hYScan0 = new TH1D("hYScan0", "Detector 0;y position of #mu (cm);A.U.", nevents, axisarray);


    // Loop through the events and the PD signals
    for (int ievent = 0; ievent < nevents; ievent++) {
        mctree->GetEntry(ievent);
        pdtree->GetEntry(iPD + ievent*numPDs_in_file);

        if (EventID != ievent+1){
            cout << "EventID mismatch: pdtree.EventID =" << EventID << "ievent+1 =" << ievent+1 << endl;
        }
        
        // Store the y (1) start position of the muon
        double ypos = StartXYZT[1];
        
        // Fill the histogram at this position with the number of counts detected
        hYScan0->Fill(ypos, CountDetected);
    }
    
    hYScan0->Draw("hist");
}
