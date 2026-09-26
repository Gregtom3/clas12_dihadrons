#include "../src/MultiplicityBinning.C"

int multiplicity(const char * input_file = "../projects/ana_v2/volatile/data/piplus_piminus/nSidis_RGA_5032.root") {
    
  
    // Define a rootRegion BinRegion with Q2 boundaries
    auto rootRegion = std::make_shared<BinRegion>();
    // Q2_max = x*ymax*s = x*16.623531870004
    rootRegion->addBoundary("[Q2]<[X]*16.623531870004");
    rootRegion->addBoundary("[Q2]> 0.9253 + 3.8613 * [X] + -9.4351 * [X]*[X]+ 16.9292 * [X]*[X]*[X] + 17.2332 * [X]*[X]*[X]*[X]");

    // Create a vector of z bin edges and use createCustomBins to create zBins
    std::vector<float> zBinEdges = {0.0, 0.5, 1.0};
    std::vector<std::shared_ptr<BinRegion>> zBins = createCustomBins("[Z]", zBinEdges);
    // Insert each zBin as a child of rootRegion
    for (const auto& zBin : zBins) {
        rootRegion->insertChild(zBin);
    }

    // Create a BinMap with rootRegion as the root BinRegion
    BinMap binMap;
    binMap.setRoot(rootRegion);

    // Open the input file and get the dihadron TTree
    TFile* tfile = new TFile(input_file);
    TTree* tree = (TTree*) tfile->Get("dihadron");

    // Set the branch addresses for the specified doubles
    double x, y, xF1, xF2, Q2, pT, Mx, z, Mh;
    tree->SetBranchAddress("x", &x);
    tree->SetBranchAddress("y", &y);
    tree->SetBranchAddress("xF1", &xF1);
    tree->SetBranchAddress("xF2", &xF2);
    tree->SetBranchAddress("Q2", &Q2);
    tree->SetBranchAddress("pTtot", &pT);
    tree->SetBranchAddress("Mx", &Mx);
    tree->SetBranchAddress("z", &z);
    tree->SetBranchAddress("Mh", &Mh);

    // Get the number of entries in the TTree
    int nentries = tree->GetEntries();

    // Loop over the entries in the TTree
    cout << "Number of Entries = " << nentries << endl;
    for (int i = 0; i < nentries; ++i) {
        tree->GetEntry(i);
        std::map<std::string, double> point = {{"X", x}, {"Q2", Q2}, {"Z", z}, {"pT", pT}};
        int childId = binMap.findLowestLevelChild(point);
        
        if(i==100){break;}
    }
    
    tfile->Close();
    return 0;
}