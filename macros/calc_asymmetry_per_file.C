#include "../src/ParseBinYAML.C"
#include "../src/fitTools.C"
#include "../src/Constants.h"
#include "../src/ParseText.C"

void create_sideband(const char * infile, std::string output_subdir, int runNumber, int pid_h1, int pid_h2);

void create_asym(const char *infile, std::string output_subdir, std::string version, int pid_h1, int pid_h2, std::string fit_region);

int calc_asymmetry_per_file(const char *infile = "/work/clas12/users/gmat/clas12/clas12_dihadrons/projects/pipi0_paper_small_w_cuts/volatile/data/piplus_pi0/nSidis_RGA_5032.root", 
                            std::string output_dir = "./tmp"
                           )
{
    // Get the pid's from the input file
    int pid_h1=0;
    int pid_h2=0;
    std::string hadron_pair="";
    getPIDs(std::string(infile),pid_h1,pid_h2,hadron_pair);
    
    // Get the run number from the input file
    int runNumber=extractRunNumberFromDataFile(infile);
    // Get the version (ex: Fall2018_RGA_inbending) from the run number
    std::string version = runPeriod(runNumber);
    
    // Generate the output directory and file
    std::string output_subdir = Form("%s/%s",output_dir.c_str(),hadron_pair.c_str()); 
    gSystem->mkdir(output_subdir.c_str());
    output_subdir += Form("/run_%d",runNumber);
    gSystem->mkdir(output_subdir.c_str());
    // Do sideband subtraction first if necessary
    if(pid_h1==111 || pid_h2==111){
        create_sideband(infile, output_subdir.c_str(), runNumber, pid_h1,pid_h2);
        create_asym(infile, output_subdir, version, pid_h1, pid_h2, "signal");
        create_asym(infile, output_subdir, version, pid_h1, pid_h2, "bkg");
    }else{
        create_asym(infile, output_subdir, version, pid_h1, pid_h2, "");
    }
    
    return 0;
}

void create_asym(const char *infile, std::string output_subdir, std::string version, int pid_h1, int pid_h2, std::string fit_region){
    // Now perform the maximum likelihood asymmetry fits
    FitManager FM;
    FM.SetUp().SetOutDir(Form("%s/asym/",output_subdir.c_str()));
    process_azi_FM(FM,version,"hel");
    if(fit_region=="signal"){
        FM.SetUp().LoadVariable("M2[0.106,0.166]");
    } else if(fit_region=="bkg"){
        FM.SetUp().LoadVariable("M2[0.2,0.4]");
    }
    
    FM.LoadData("dihadron_cuts",infile);
    Here::Go(&FM);
    
    // Open the output TFile
    // Rename the output TFile first
    gSystem->Rename(Form("%s/asym/ResultsHSMinuit2.root",output_subdir.c_str()),Form("%s/asym/%sResultsHSMinuit2.root",output_subdir.c_str(),fit_region.c_str()));
    
    TFile * resultFile = new TFile(Form("%s/asym/%sResultsHSMinuit2.root",output_subdir.c_str(),fit_region.c_str()),"READ");
    TTree * ResultTree = (TTree*)resultFile->Get("ResultTree");
    
    // Save the branch information to a YAML file
    std::ofstream outFile;
    outFile.open(Form("%s/asym/%sResultsHSMinuit2.yaml",output_subdir.c_str(),fit_region.c_str()));

    TObjArray *branches = ResultTree->GetListOfBranches();
    int nBranches = branches->GetEntries();

    for (int i = 0; i < nBranches; ++i) {
        TBranch *branch = (TBranch*) branches->At(i);
        const char* branchName = branch->GetName();

        // Create a double variable to get the branch value
        double branchVal;
        ResultTree->SetBranchAddress(branchName, &branchVal);

        // Ensure we are at the first entry of the TTree
        ResultTree->GetEntry(0);

        outFile << branchName << ": " << std::setprecision(4) << branchVal << "\n";
    }

    outFile.close();    
}

void create_sideband(const char * infile, std::string output_subdir, int runNumber, int pid_h1, int pid_h2){
    
    // Open the TFile
    TFile *tfile = new TFile(infile,"READ");
    // Open the TTree
    TTree *tree = (TTree*)tfile->Get("dihadron_cuts");
    
    // Make sideband info directory
    output_subdir+="/sideband";
    gSystem->mkdir(output_subdir.c_str());
    
    // From PID's determine key features of the sideband fitting
    std::string cvar = "";
    double x1,x2,x3,x4;
    float s1,s2;
    std::vector<std::string> thresh_cuts;
    if((pid_h1==111&&pid_h2!=111) || (pid_h1!=111&&pid_h2==111)){
        x1=0.07;
        x2=0.22;
        x3=0.22;
        x4=0.4;
        if(pid_h1==111){
            cvar="M1";
        }
        else{
            cvar="M2";  
        } 
        s1=0.007;
        s2=0.04;
        
    }else if(pid_h1==111&&pid_h2==111){
        x1=0.14;
        x2=0.44;
        x3=0.44;
        x4=0.8;
        cvar="M12";
        s1=0.014;
        s2=0.08;
    }
    
    
    TH1F *h = new TH1F("Mdiphoton",";Mdiphoton [GeV];",100,x1,x4);
    
    // Draw into the histogram
    tree->Draw(Form("%s>>Mdiphoton",cvar.c_str()),"","goff");

    // Scale the histogram
    TH1F *hnorm = (TH1F*)h->Clone();
    hnorm->SetName("Mdiphoton_normed");
    hnorm->Scale(1/hnorm->Integral());

    ///////////////////////////// Find the minimum and maximum filled bin for bin range
    double xmin = 0;
    double xmax = x4;
    for(int i = 1 ; i <= h->GetNbinsX(); i++){
        float content = h->GetBinContent(i);
        if(content>0 && xmin==0){
            xmin=h->GetBinLowEdge(i);
        }
        if(content==0 && h->GetBinCenter(i)>x2){
            xmax = h->GetBinLowEdge(i+1);
            break;
        }
    }

    ///////////////////////////// Create fit object 
    TF1 * f_sdbnd = new TF1("f_sdbnd","gaus(0)+pol4(3)",xmin,xmax);
    f_sdbnd->SetParameter(1,(x1+x2)/2);
    f_sdbnd->SetParLimits(0,0.001,100);
    f_sdbnd->SetParLimits(1,x1,x2);
    f_sdbnd->SetParLimits(2,s1,s2);

    ///////////////////////////// Perform fit
    hnorm->Fit(f_sdbnd,"R");

    // Calculate mean and +- 2sigma
    double Mggmin = 0.106;
    double Mggmax = 0.166;
    if(pid_h1==pid_h2&&pid_h1==111){
        Mggmin*=2;
        Mggmax*=2;
    }


    ///////////////////////////// Calculate purity
    TF1 * sig = new TF1("sig","gaus(0)",Mggmin,Mggmax);
    TF1 * bg = new TF1("bg","pol4(0)",Mggmin,Mggmax);
    for(int j = 0 ; j < 3 ; j++){
        sig->SetParameter(j,f_sdbnd->GetParameter(j));
    }
    for(int k = 0 ; k < 5 ; k++){
        bg->SetParameter(k,f_sdbnd->GetParameter(k+3));
    }

    // Calculate integrals
    double signal_integral = sig->Integral(Mggmin, Mggmax);
    double background_integral = bg->Integral(Mggmin, Mggmax);
    double normalization_integral = hnorm->Integral(hnorm->FindBin(Mggmin), hnorm->FindBin(Mggmax));
    double f_sdbnd_integral = f_sdbnd->Integral(Mggmin, Mggmax);

    // Calculate ratios
    double u1 = signal_integral / normalization_integral;
    double u2 = signal_integral / f_sdbnd_integral;
    double u3 = 1 - (background_integral / normalization_integral);
    double u4 = 1 - (background_integral / f_sdbnd_integral);

    // Create a YAML file for the output
    std::ofstream outFile;
    outFile.open(output_subdir + "/sideband.yaml");

    outFile << "purity_1: " << std::setprecision(4) << u1 << "\n";
    outFile << "purity_2: " << std::setprecision(4) << u2 << "\n";
    outFile << "purity_3: " << std::setprecision(4) << u3 << "\n";
    outFile << "purity_4: " << std::setprecision(4) << u4 << "\n";
    outFile << "signal_region: [" << Mggmin << ", " << Mggmax << "]\n";

    // Assuming that f_sdbnd is not null and properly initialized.
    std::vector<double> parameters(8);
    std::vector<double> errors(8);
    for (int i = 0; i < 8; ++i) {
        parameters[i] = f_sdbnd->GetParameter(i);
        errors[i] = f_sdbnd->GetParError(i);
    }

    outFile << "parameters: [";
    for (const auto& p : parameters) {
        outFile << p << ", ";
    }
    outFile.seekp(-2, std::ios_base::end);  // To remove the last comma and space.
    outFile << "]\n";

    outFile << "errors: [";
    for (const auto& e : errors) {
        outFile << e << ", ";
    }
    outFile.seekp(-2, std::ios_base::end);  // To remove the last comma and space.
    outFile << "]\n";

    outFile.close();
}


