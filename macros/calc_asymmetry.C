#include "../src/ParseBinYAML.C"
#include "../src/ParseCutYAML.C"
#include "../src/fitTools.C"
#include "../src/Constants.h"
#include "../src/ParseText.C"
#include "../src/injectDihadronAsym.C"

enum ASYM_TYPE {
    AZI,
    TWOH
};

enum ANA_TYPE {
    SPLOT,
    SIDEBAND,
    STANDARD
};


void create_sweights(const char * infile, const char * brudir, YAMLbinstruct binStruct, Block cut, int pid_h1, int pid_h2, bool use_ML);
void create_sideband(const char * infile, const char * brudir, YAMLbinstruct binStruct, Block cut, int pid_h1, int pid_h2, bool use_ML);
void asym(const char * infile, const char * brudir, std::string version, YAMLbinstruct binStruct, Block cut, int pid_h1, int pid_h2, bool use_ML, ANA_TYPE ana_type,  ASYM_TYPE asym_type);

int calc_asymmetry(//const char * infile = "/volatile/clas12/users/gmat/clas12analysis.sidis.data/clas12_dihadrons/projects/ana_v0/data/piplus_piminus/MC_RGA_inbending_merged.root",
                   //const char * infile = "./MC_RGA_inbending_merged.root",
		   const char * infile = "/volatile/clas12/users/gmat/clas12analysis.sidis.data/clas12_dihadrons/projects/ana_v2/data/piplus_pi0/nSidis_RGA_5036.root",
		   const char * binfile = "/work/clas12/users/gmat/clas12/clas12_dihadrons/utils/binning_files/Binning.yaml",
		   //                   std::string brudir  = "/work/clas12/users/gmat/clas12/clas12_dihadrons/projects/ana_v0/asym",
		   std::string brudir  = "./bru",
                   const int binscheme = 0,
                   const std::string cut_title = "v1",
                   bool use_ML = false,
                   bool do_inject = false,
                   bool create_splot = false,
                   bool do_sweighted = false,
                   bool do_sideband = false,
                   bool remove_inject = false){
    
    
    // Pull the /<h1>_<h2>/ from the "infile"
    // and determine the PID's of the hadron from the file path
    int pid_h1=0;
    int pid_h2=0;
    std::string hadron_pair="";
    getPIDs(std::string(infile),pid_h1,pid_h2,hadron_pair);
    // Pull the <version>_merged.root from "infile"
    std::string version=getVersion(infile);
    
    // Make the first subdirectory 
    brudir+="/"+version;
    gSystem->mkdir(TString(brudir)); // <project>/asym/<version>
    
    // Make the second subdirectory 
    brudir+="/cut_"+cut_title;
    gSystem->mkdir(TString(brudir)); // <project>/asym/<version>/<cut>
    
    // Create a directory for the hadron pair
    brudir+="/"+hadron_pair; 
    gSystem->mkdir(TString(brudir)); // <project>/asym/<version>/<cut>/<hadron_pair>
    // Pull the bin structures from the yaml file
    auto binStructs = get_structs(binfile);
    auto binStruct = binStructs[binscheme];
    // If this is a Monte Carlo file, then we must inject the asymmetries
    std::string INFILE="";
    if (strstr(infile, "MC_") != NULL){
      INFILE=std::string(infile) + ".inject." + binStruct.name + std::string(use_ML == true ? ".ML." : ".noML.")  +  std::string(cut_title) + ".splot.root";
        // Delete the injected Monte carlo file after use
        if(remove_inject){   
            if (strstr(infile, "MC_") != NULL){
                gSystem->Unlink(TString(INFILE));
            }
            return 0;
        }
        else if(do_inject){
            injectDihadronAsym(infile,INFILE.c_str(),binStruct,pid_h1,pid_h2,2);
        }
        
    }else{
        INFILE=std::string(infile);
    }


    
    // Create a directory for the binning scheme selected by binscheme
    brudir+="/";
    for (int i = 0; i < binStruct.numDimensions; i++) {
        std::string binName = binStruct.dimensionNames[i];
        brudir+=binName;
        if(i!=binStruct.numDimensions-1){
            brudir+="_";
        }
    }
    gSystem->mkdir(TString(brudir)); // <project>/asym/<version>/<hadron_pair>/<binning>
    
    // Create a directory for ML usage
    if(use_ML==true){
        brudir+="/ML";
    }else{
        brudir+="/noML";
    }
    gSystem->mkdir(TString(brudir)); // <project>/asym/<version>/<cut>/<hadron_pair>/<binning>/<ML>
    // Save the YAML bin into the brudir
    serializeYAMLbinstruct(binStruct,brudir+"/Binning.yaml");
    
    // Determine the auxillary cuts to place
    const std::string cut_file = "/work/clas12/users/gmat/clas12/clas12_dihadrons/utils/cut_library.yaml";
    auto cuts = readYaml(cut_file);
    Block cut;
    for (auto c: cuts){
        if(c.title==cut_title){
            cut=c;
            break;
        }
    }
    
    if(cut.title!=cut_title || cut_title==""){
        cout << "Error...unsure what cuts to use if any...Aborting..." << endl;
        return -1;
    }
    
    // Save auxillary cut file to folder
    writeBlockToFile(brudir+"/Cut.yaml",cut);
    
    
    
    // ----------------------------------------------------------------------------------
    // Fitting Code
    // ----------------------------------------------------------------------------------
    
    
    
    // For loop over dihadron asym approaches
    //std::vector<ASYM_TYPE> asym_types = {AZI,TWOH};
    std::vector<ASYM_TYPE> asym_types = {AZI};
    for(ASYM_TYPE asym_type: asym_types){
        
        // Create subdirectories for each of the asymmetry types
        std::string BRUDIR=brudir+"/";
        if(asym_type==AZI) BRUDIR+="AZI";
        else if(asym_type==TWOH) BRUDIR+="TWOH";
        gSystem->mkdir(TString(BRUDIR)); // <project>/asym/<version>/<cut>/<hadron_pair>/<binning>/<ML>/<asym_type>
        
        
        
        // Based on the pid combination, determine if we need background subtraction
        // If we are using pi0's
        if(pid_h1==111 || pid_h2==111){
            if(do_sideband){
                create_sideband(INFILE.c_str(),BRUDIR.c_str(),binStruct, cut,pid_h1,pid_h2,use_ML);
                asym(INFILE.c_str(), BRUDIR.c_str(), version, binStruct, cut, pid_h1, pid_h2, use_ML, SIDEBAND,  asym_type);
            }
            if(create_splot){
                create_sweights(INFILE.c_str(),BRUDIR.c_str(),binStruct, cut,pid_h1,pid_h2,use_ML);
            }
            if(do_sweighted){
                asym(INFILE.c_str(), BRUDIR.c_str(), version, binStruct, cut, pid_h1, pid_h2, use_ML, SPLOT,  asym_type);
            }
        }
        else{
            asym(INFILE.c_str(), BRUDIR.c_str(), version, binStruct, cut, pid_h1, pid_h2, use_ML, STANDARD, asym_type);
        }
        
        break;
    }
    

    
    return 0;
}







void create_sweights(const char * infile, const char * brudir, YAMLbinstruct binStruct, Block cut, int pid_h1, int pid_h2, bool use_ML){
    
    const double p_thresh = 0.9;
    
    // Create sPlot
    sPlot RF;
    RF.SetUp().SetOutDir(Form("%s/outsPlotBins/",brudir));
    
    //////////////////////////////////// Load variables
    RF.SetUp().LoadVariable("phi_h[-100,100]");
    RF.SetUp().LoadVariable("phi_R0[-100,100]");
    RF.SetUp().LoadVariable("th[-100,100]");
    RF.SetUp().LoadVariable("phi_h1[-100,100]");
    RF.SetUp().LoadVariable("phi_h2[-100,100]");
    RF.SetUp().SetIDBranchName("fggID");
    
    // From PID's determine several properties of the sWeight fitting
    std::string cvar = "";
    float u1,u2,u3; // mean vars
    float s1,s2,s3; // sigma vars
    if(pid_h1==111 && pid_h2!=111){
        cvar = "M1";
        RF.SetUp().LoadVariable("M1[0.07,0.4]");
        u1=0.131; u2=0.129; u3=0.15;
        s1=0.01;  s2=0.0001; s3=0.1;
        if(use_ML){
            RF.SetUp().LoadVariable(Form("p_11[%f,1]",p_thresh));   
            RF.SetUp().LoadVariable(Form("p_12[%f,1]",p_thresh)); 
            RF.SetUp().AddCut(Form("p_11>%f",p_thresh));
            RF.SetUp().AddCut(Form("p_12>%f",p_thresh));
        }else{
            RF.SetUp().LoadVariable("isGoodEventWithoutML[0.9,1.1]");
            RF.SetUp().AddCut("isGoodEventWithoutML==1");
        }
    }
    else if(pid_h1!=111 && pid_h2==111){
        cvar = "M2";
        RF.SetUp().LoadVariable("M2[0.07,0.4]");
        u1=0.131; u2=0.129; u3=0.15;
        s1=0.01;  s2=0.0001; s3=0.1;
        if(use_ML){
            RF.SetUp().LoadVariable(Form("p_21[%f,1]",p_thresh)); 
            RF.SetUp().LoadVariable(Form("p_22[%f,1]",p_thresh)); 
            RF.SetUp().AddCut(Form("p_21>%f",p_thresh));
            RF.SetUp().AddCut(Form("p_22>%f",p_thresh));
        }else{
            RF.SetUp().LoadVariable("isGoodEventWithoutML[0.9,1.1]");
            RF.SetUp().AddCut("isGoodEventWithoutML==1");
        }
    }
    else if(pid_h1==111 && pid_h2==111){
        cvar = "M12";
        RF.SetUp().LoadVariable("M12[0.14,0.8]");
        u1=0.262; u2=0.258; u3=0.3;
        s1=0.02;  s2=0.0002; s3=0.2;
        if(use_ML){
            RF.SetUp().LoadVariable(Form("p_11[%f,1]",p_thresh));   
            RF.SetUp().LoadVariable(Form("p_12[%f,1]",p_thresh)); 
            RF.SetUp().LoadVariable(Form("p_21[%f,1]",p_thresh)); 
            RF.SetUp().LoadVariable(Form("p_22[%f,1]",p_thresh)); 
            RF.SetUp().AddCut(Form("p_11>%f",p_thresh));
            RF.SetUp().AddCut(Form("p_12>%f",p_thresh));
            RF.SetUp().AddCut(Form("p_21>%f",p_thresh));
            RF.SetUp().AddCut(Form("p_22>%f",p_thresh));
        }else{
            RF.SetUp().LoadVariable("isGoodEventWithoutML[0.9,1.1]");
            RF.SetUp().AddCut("isGoodEventWithoutML==1");
        }
    }
    
    // For loop over the cuts desired to apply them
    for(int j = 0 ; j < cut.vmin.size() ; j++){                    
        RF.SetUp().LoadVariable(Form("%s[%f,%f]",cut.var.at(j).c_str(),cut.vmin.at(j),cut.vmax.at(j)));
        RF.SetUp().AddCut(Form("%s>%f",cut.var.at(j).c_str(),cut.vmin.at(j)));
        RF.SetUp().AddCut(Form("%s<%f",cut.var.at(j).c_str(),cut.vmax.at(j)));
    }
    
    //////////////////////////////////// Gaussian
    RF.SetUp().FactoryPDF(Form("Gaussian::Signal( %s, mean[%f,%f,%f], sigma[%f,%f,%f] )",cvar.c_str(),u1,u2,u3,s1,s2,s3));
    RF.SetUp().LoadSpeciesPDF("Signal",1);
  //////////////////////////////////// Background
    RF.SetUp().FactoryPDF(Form("Chebychev::BG( %s, {a0[-0.1,-1,1], a1[0.1,-1,1], a2[-0.1,-1,1], a3[-0.1,-1,1]} )",cvar.c_str()));
    RF.SetUp().LoadSpeciesPDF("BG",1);
    
    // Load data
    RF.LoadData("dihadron",infile);
    
    //////////////////////////////////// Bins
    // Not sure if this is still needed?
    // Load bin variables
    for (int i = 0; i < binStruct.numDimensions; i++) {
      std::string binName = binStruct.dimensionNames[i];
      std::vector<double> binEdges = binStruct.binEdges[i];
      int numBins = binEdges.size();
      Double_t binEdgesArr[numBins];
      for (int j = 0; j < numBins; j++) {
        binEdgesArr[j] = binEdges[j];
      }
      RF.Bins().LoadBinVar(binName, numBins-1, binEdgesArr);
    }
  
    // Run sWeighting procedure
    Here::Go(&RF);  
    
    // Create TCanvas and Plot Signal distribution
    new TCanvas();
    RF.DrawWeighted(TString(cvar),"Signal");
    RF.DeleteWeightedTree();
}


void create_sideband(const char * infile, const char * brudir, YAMLbinstruct binStruct, Block cut, int pid_h1, int pid_h2, bool use_ML){
    
    const double p_thresh = 0.9; // ML cut
    std::string hel_str="hel";  // changes if we inject the Monte Carlo
    if (strstr(infile, "MC_") != NULL){
        hel_str=binStruct.injectName+"."+hel_str;
    }
    
    // Open the TFile
    TFile *tfile = new TFile(infile,"READ");
    // Open the TTree
    TTree *tree = (TTree*)tfile->Get("dihadron");
    
    // Make sideband info directory
    gSystem->mkdir(Form("%s/outSdbndBins/",brudir));
    
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
            thresh_cuts.push_back(Form("p_11>%f",p_thresh));
            thresh_cuts.push_back(Form("p_12>%f",p_thresh));
        }
        else{
            cvar="M2";  
            thresh_cuts.push_back(Form("p_21>%f",p_thresh));
            thresh_cuts.push_back(Form("p_22>%f",p_thresh));
        } 
        s1=0.007;
        s2=0.04;
        
    }else if(pid_h1==111&&pid_h2==111){
        x1=0.14;
        x2=0.44;
        x3=0.44;
        x4=0.8;
        cvar="M12";
        thresh_cuts.push_back(Form("p_11>%f",p_thresh));
        thresh_cuts.push_back(Form("p_12>%f",p_thresh));
        thresh_cuts.push_back(Form("p_21>%f",p_thresh));
        thresh_cuts.push_back(Form("p_22>%f",p_thresh));
        s1=0.014;
        s2=0.08;
    }
    
    
    TH1F *h = new TH1F("Mdiphoton",";Mdiphoton [GeV];",100,x1,x4);
    
    
    // perform sideband for each of the binnings
    const int N=binStruct.nBins;
    for(int n = 0; n < N; n++){
        
        h->Reset(); // Clear binContents but keep bin sizes and range'
        
        std::string subdir = get_dir_from_binstruct_idx(binStruct,n); // x_<min>_<max>_y_<min>_<max>_z_<min>_<max>
        std::string cutstr = get_cut_from_binstruct_idx(binStruct,n); // x>min&&x<max&&...
        cout << "Parsing " << n+1 << "/" << N << "    " << subdir << endl;
        
        // Make directory to store sideband fit info
        gSystem->mkdir(Form("%s/outSdbndBins/%s",brudir,subdir.c_str()));
        
        // Add helicity requirements
        cutstr += "&&("+hel_str+"==1||"+hel_str+"==-1)";
        
        // Determine additional cuts based on ML usage
        if(use_ML){
            for(auto s: thresh_cuts){
                cutstr+="&&"+s;
            }
        }else{
            cutstr+="&&isGoodEventWithoutML==1";
        }
        
        // For loop over the yaml cuts desired to apply them
        for(int j = 0 ; j < cut.vmin.size() ; j++){  
            cutstr+="&&"+cut.var.at(j)+">"+std::to_string(cut.vmin.at(j))+"&&"+cut.var.at(j)+"<"+std::to_string(cut.vmax.at(j));
        }
        
        // Draw into the histogram
        tree->Draw(Form("%s>>Mdiphoton",cvar.c_str()),TString(cutstr),"goff");
        
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
        
        ///////////////////////////// Create TFile for saving output
        TFile *sdbnd_out = new TFile(Form("%s/outSdbndBins/%s/sideband.root",brudir,subdir.c_str()),"RECREATE");
        TVectorD purity_1(1);
        TVectorD purity_2(1);
        TVectorD purity_3(1);
        TVectorD purity_4(1);
        purity_1[0]=u1;
        purity_2[0]=u2;
        purity_3[0]=u3;
        purity_4[0]=u4;
        TVectorD sigregion(2);
        sigregion[0]=Mggmin;
        sigregion[1]=Mggmax;
        sdbnd_out->WriteObject(&purity_1,"purity_1");
        sdbnd_out->WriteObject(&purity_2,"purity_2");
        sdbnd_out->WriteObject(&purity_3,"purity_3");
        sdbnd_out->WriteObject(&purity_4,"purity_4");
        sdbnd_out->WriteObject(&sigregion,"signal_region");
        f_sdbnd->Write();
        hnorm->Write();
        h->Write();
        sdbnd_out->Close();
        delete hnorm;
        delete f_sdbnd;
        delete sig;
        delete bg;
        hnorm=NULL;
        f_sdbnd=NULL;
        sig=NULL;
        bg=NULL;
    }
}


void asym(const char * infile, const char * brudir, std::string version, YAMLbinstruct binStruct, Block cut, int pid_h1, int pid_h2, bool use_ML, ANA_TYPE ana_type,  ASYM_TYPE asym_type){

    const double p_thresh = 0.9; // ML cut
    std::string hel_str="hel";  // changes if we inject the Monte Carlo
    if (strstr(infile, "MC_") != NULL){
        hel_str=binStruct.injectName+"."+hel_str;
    }
    
    if(ana_type==SPLOT){
        
        for(int reg = 0 ; reg < 2; reg++){
            FitManager FM;
            if(reg==0)
                FM.SetUp().SetOutDir(Form("%s/outObsBins_splot_sig/",brudir));
            else
                FM.SetUp().SetOutDir(Form("%s/outObsBins_splot_bg/",brudir));
            if(asym_type==AZI)
                process_azi_FM(FM,version,hel_str);
            else if(asym_type==TWOH)
                process_2h_FM(FM,version,hel_str);
            
            //////////////////////////////////// Bins
            //Load bin variables
            for (int i = 0; i < binStruct.numDimensions; i++) {
                std::string binName = binStruct.dimensionNames[i];
                std::vector<double> binEdges = binStruct.binEdges[i];
                int numBins = binEdges.size();
                Double_t binEdgesArr[numBins];
                for (int j = 0; j < numBins; j++) {
                    binEdgesArr[j] = binEdges[j];
                }
                FM.Bins().LoadBinVar(binName, numBins-1, binEdgesArr);
            }
            
            FM.LoadData("dihadron",infile);
            if(reg==0)
                FM.Data().LoadWeights("Signal",Form("%s/outsPlotBins/Tweights.root",brudir));
            else
                FM.Data().LoadWeights("BG",Form("%s/outsPlotBins/Tweights.root",brudir));
            
            
            Here::Go(&FM);
        }
    }
    
    else if(ana_type==SIDEBAND){
        // Do signal+bkg and bkg regions separately
        for(int reg=0;reg<2;reg++){
            FitManager FM;

            
            // We need to reapply the cuts made to generate the sideband plots
            if(reg==0){
                FM.SetUp().SetOutDir(Form("%s/outObsBins_sdbnd_sigbg/",brudir));
            }else{
                FM.SetUp().SetOutDir(Form("%s/outObsBins_sdbnd_bg/",brudir));
            }
            
            if(asym_type==AZI)
                process_azi_FM(FM,version,hel_str);
            else if(asym_type==TWOH)
                process_2h_FM(FM,version,hel_str);
            
            if(reg==0 && pid_h1==111 && pid_h2!=111){
                FM.SetUp().LoadVariable("M1[0.106,0.166]");
            }else if(reg==1 && pid_h1==111 && pid_h2!=111){
                FM.SetUp().LoadVariable("M1[0.2,0.4]");
            }else if(reg==0 && pid_h1!=111 && pid_h2==111){
                FM.SetUp().LoadVariable("M2[0.106,0.166]");
            }else if(reg==1 && pid_h1!=111 && pid_h2==111){
                FM.SetUp().LoadVariable("M2[0.2,0.4]");
            }else if(reg==0 && pid_h1==111 && pid_h2==111){
                FM.SetUp().LoadVariable("M12[0.212,0.332]");
            }else if(reg==1 && pid_h1==111 && pid_h2==111){
                FM.SetUp().LoadVariable("M12[0.4,0.8]");
            }
            // From PID's determine ML cuts
            if(pid_h1==111 && pid_h2!=111){
                if(use_ML){
                    FM.SetUp().LoadVariable(Form("p_11[%f,1]",p_thresh));   
                    FM.SetUp().LoadVariable(Form("p_12[%f,1]",p_thresh)); 
                    FM.SetUp().AddCut(Form("p_11>%f",p_thresh));
                    FM.SetUp().AddCut(Form("p_12>%f",p_thresh));
                }else{
                    FM.SetUp().LoadVariable("isGoodEventWithoutML[0.9,1.1]");
                    FM.SetUp().AddCut("isGoodEventWithoutML==1");
                }
            }
            else if(pid_h1!=111 && pid_h2==111){
                if(use_ML){
                    FM.SetUp().LoadVariable(Form("p_21[%f,1]",p_thresh)); 
                    FM.SetUp().LoadVariable(Form("p_22[%f,1]",p_thresh)); 
                    FM.SetUp().AddCut(Form("p_21>%f",p_thresh));
                    FM.SetUp().AddCut(Form("p_22>%f",p_thresh));
                }else{
                    FM.SetUp().LoadVariable("isGoodEventWithoutML[0.9,1.1]");
                    FM.SetUp().AddCut("isGoodEventWithoutML==1");
                }
            }
            else if(pid_h1==111 && pid_h2==111){
                if(use_ML){
                    FM.SetUp().LoadVariable(Form("p_11[%f,1]",p_thresh));   
                    FM.SetUp().LoadVariable(Form("p_12[%f,1]",p_thresh)); 
                    FM.SetUp().LoadVariable(Form("p_21[%f,1]",p_thresh)); 
                    FM.SetUp().LoadVariable(Form("p_22[%f,1]",p_thresh)); 
                    FM.SetUp().AddCut(Form("p_11>%f",p_thresh));
                    FM.SetUp().AddCut(Form("p_12>%f",p_thresh));
                    FM.SetUp().AddCut(Form("p_21>%f",p_thresh));
                    FM.SetUp().AddCut(Form("p_22>%f",p_thresh));
                }else{
                    FM.SetUp().LoadVariable("isGoodEventWithoutML[0.9,1.1]");
                    FM.SetUp().AddCut("isGoodEventWithoutML==1");
                }
            }


            // For loop over the cuts desired to apply them
            for(int j = 0 ; j < cut.vmin.size() ; j++){                    
               FM.SetUp().LoadVariable(Form("%s[%f,%f]",cut.var.at(j).c_str(),cut.vmin.at(j),cut.vmax.at(j)));    
               FM.SetUp().AddCut(Form("%s>%f",cut.var.at(j).c_str(),cut.vmin.at(j)));
               FM.SetUp().AddCut(Form("%s<%f",cut.var.at(j).c_str(),cut.vmax.at(j)));  
            }
            
            //////////////////////////////////// Bins
            //Load bin variables
            for (int i = 0; i < binStruct.numDimensions; i++) {
                std::string binName = binStruct.dimensionNames[i];
                std::vector<double> binEdges = binStruct.binEdges[i];
                int numBins = binEdges.size();
                Double_t binEdgesArr[numBins];
                for (int j = 0; j < numBins; j++) {
                    binEdgesArr[j] = binEdges[j];
                }
                FM.Bins().LoadBinVar(binName, numBins-1, binEdgesArr);
            }
        
            FM.LoadData("dihadron",infile);
            Here::Go(&FM);
        }
    }
    
    
    // No background subtraction needed
    else if(ana_type==STANDARD){
        FitManager FM;
        FM.SetUp().SetOutDir(Form("%s/outObsBins/",brudir));
        if(asym_type==AZI)
            process_azi_FM(FM,version,hel_str);
        else if(asym_type==TWOH)
            process_2h_FM(FM,version,hel_str);
        
        // For loop over the cuts desired to apply them
        for(int j = 0 ; j < cut.vmin.size() ; j++){
            FM.SetUp().LoadVariable(Form("%s[%f,%f]",cut.var.at(j).c_str(),cut.vmin.at(j),cut.vmax.at(j)));
            FM.SetUp().AddCut(Form("%s>%f",cut.var.at(j).c_str(),cut.vmin.at(j)));
            FM.SetUp().AddCut(Form("%s<%f",cut.var.at(j).c_str(),cut.vmax.at(j)));
        }
        
        //////////////////////////////////// Bins
        //Load bin variables
        for (int i = 0; i < binStruct.numDimensions; i++) {
            std::string binName = binStruct.dimensionNames[i];
            std::vector<double> binEdges = binStruct.binEdges[i];
            int numBins = binEdges.size();
            Double_t binEdgesArr[numBins];
            for (int j = 0; j < numBins; j++) {
                binEdgesArr[j] = binEdges[j];
            }
            FM.Bins().LoadBinVar(binName, numBins-1, binEdgesArr);
        }
        
        FM.LoadData("dihadron",infile);
        Here::Go(&FM);
    }
    
}
