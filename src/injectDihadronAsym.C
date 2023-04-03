// Get Psi_i's for the A_LU
std::vector<TF2> get_mods(int L){
    
    int mod = 0;
    std::vector<TF2> mods;
    std::vector<string> str_vec;
    
    // Get strings of all the functions of the azimuthals
    // Useful for removing duplicates later
    for (int l = 0; l <= L; l++)
    {
        for (int m = 1; m <= l; m++)
        {
            string func = Form("sin(%d*x - %d*y)",m,m);
            str_vec.push_back(func);
            
        }
        for (int m = -l; m <= l; m++)
        {
            string func = Form("sin(%d*x + %d*y)",1-m,m);
            str_vec.push_back(func);
        }
    }
    
    // Remove duplicate entries
    sort(str_vec.begin(), str_vec.end());
    str_vec.erase(unique(str_vec.begin(), str_vec.end()), str_vec.end());
    
    // Fill the TF2 vector
    for(unsigned int i = 0 ; i < str_vec.size(); i++){
        mods.push_back(TF2(Form("mod%d",mod),str_vec.at(i).c_str(),-3.15,3.15,-3.15,3.15));
        mod++;
    }
    return mods;
}

// Get the A_LUs
std::vector<TF2> get_ALUs_sig(YAMLbinstruct bs){
    int mod = 0;
    std::vector<TF2> ALUs;
    for(std::string ALU_func: bs.injectsigFuncs){
        ALUs.push_back(TF2(Form("ALUsig%d",mod),ALU_func.c_str(),-100,100,-100,100));
        mod++;
    }
    return ALUs;
}

std::vector<TF2> get_ALUs_bg(YAMLbinstruct bs){
    int mod = 0;
    std::vector<TF2> ALUs;
    for(std::string ALU_func: bs.injectbgFuncs){
        ALUs.push_back(TF2(Form("ALUbg%d",mod),ALU_func.c_str(),-100,100,-100,100));
        mod++;
    }
    return ALUs;
}




// Inject the TTree with a helicity determined by probabilities
void injectDihadronAsym(const char * input_file, YAMLbinstruct binStruct, int pid_h1, int pid_h2, int L){
    
    TFile * f = new TFile(input_file,"UPDATE");
    std::vector<TF2> mods = get_mods(L);
    std::vector<TF2> ALUs_sig = get_ALUs_sig(binStruct);
    std::vector<TF2> ALUs_bg = get_ALUs_bg(binStruct);
    // Get the binNames from the binStruct
    std::vector<string> binNames = binStruct.dimensionNames;
    // Open the TTree
    TTree *tIn = (TTree*)f->Get("dihadron");
    // Dimension of the binning
    int dim = binStruct.numDimensions;
    // Set relevant branch addresses
    int MCmatch;
    int truepid_1;
    int truepid_2;
    double x=0;
    double y=0;
    double phih;
    double phiR;
    
    tIn->SetBranchAddress("truepid_1",&truepid_1);
    tIn->SetBranchAddress("truepid_2",&truepid_2);
    tIn->SetBranchAddress("MCmatch",&MCmatch);
    tIn->SetBranchAddress("truephi_h",&phih);
    tIn->SetBranchAddress("truephi_R0",&phiR);
    tIn->SetBranchAddress(binNames.at(0).c_str(),&x);
    if(dim>1){
        tIn->SetBranchAddress(binNames.at(1).c_str(),&y);
    }
    
    auto injectName = binStruct.injectName;

    TString hel_branchname = TString(injectName) + TString(".hel");
    // Disable the existing branch if it exists
    TBranch *b = tIn->GetBranch(hel_branchname);
    tIn->GetListOfBranches()->Remove(b);
    TLeaf* l= tIn->GetLeaf(hel_branchname);
    tIn->GetListOfLeaves()->Remove(l);
    
    // Create helicity injection
    int hel;
    TBranch * hel_branch = tIn->Branch(hel_branchname,&hel,hel_branchname+TString("/I"));
    // Loop through all events
    float sum=0.0;
    float prob=0.0;
    for ( int i = 0 ; i < tIn->GetEntries() ; i++){

        tIn->GetEntry(i);
        if(i%10000==0)
            cout << i << "/" << tIn->GetEntries()-1 << "    x = " << x << " , y = " << y << " , prob = " << prob << endl;
        
        // If not all Reco Particles had a reco match, skip event
        if(MCmatch!=1){
            hel=0;
            hel_branch->Fill();
            continue;
        }
        
        // Inject either signal or background
        sum=0.0;
        if(truepid_1==pid_h1 && truepid_2==pid_h2){
            for (unsigned int j = 0 ; j < mods.size(); j++){
                sum+=ALUs_sig.at(j).Eval(x,y)*mods.at(j).Eval(phih,phiR);
            }
        }
        else
        {
            for (unsigned int j = 0 ; j < mods.size(); j++){
                sum+=ALUs_bg.at(j).Eval(x,y)*mods.at(j).Eval(phih,phiR);
            }
        }
        
        // Calculate the probability of positive helicity
        prob = (1.0+sum)/2.0; // 0 < prob < 1 for small asymmetries
        
        // Set helicity
        if(gRandom->Uniform()<prob)
            hel=1;
        else
            hel=-1;
        hel_branch->Fill();
    }
    tIn->Print();
    tIn->Write(0,TObject::kOverwrite);
    f->Close();
}