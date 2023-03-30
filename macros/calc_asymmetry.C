#include "../src/ParseBinYAML.C"
#include "../src/ParseCutYAML.C"
#include "../src/fitTools.C"
#include "../src/Constants.h"

void create_sweights(const char * infile, const char * brudir, YAMLbinstruct bs, Block cut, int pid_h1, int pid_h2, bool use_ML);
void create_sideband();


int calc_asymmetry(const char * infile = "/volatile/clas12/users/gmat/clas12analysis.sidis.data/clas12_dihadrons/projects/ana_v0/data/pi0_pi0/Fall2018_RGA_inbending_merged.root",
                   const char * binfile = "/work/clas12/users/gmat/scipio/utils/Binning.yaml",
                   std::string brudir  = "/work/clas12/users/gmat/clas12/clas12_dihadrons/projects/ana_v0/asym"
                   const int binnum = 0,
                   const std::string cut_title = "v1",
                   bool use_ML = true){
    
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
    
    // Create a directory for the hadron pair
    brudir+="/"+hadron_pair; 
    gSystem->mkdir(TString(brudir)); // <project>/asym/<version>/<hadron_pair>
    
    // Pull the bin structures from the yaml file
    auto binStructs = get_structs(binfile);
    auto binStruct = binStructs[binnum];
    
    // Create a directory for the binning scheme selected by binnum
    brudir+="/";
    for (int i = 0; i < binStruct.numDimensions; i++) {
        std::string binName = binStruct.dimensionNames[i];
        brudir+=binName;
        if(i!=binStruct.numDimensions){
            brudir+="_";
        }
    }
    gSystem->mkdir(TString(brudir)); // <project>/asym/<version>/<hadron_pair>/<binning>
    
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
    
    // Based on the pid combination, determine if we need background subtraction
    // If we are using pi0's
    if(pid_h1==111 || pid_h2==111){
        create_sweights(infile,brudir,binStruct,cut,pid_h1,pid_h2,use_ML);
    }
    
    return 0;
}







void create_sweights(const char * infile, const char * brudir, YAMLbinstruct binStruct, Block cut, int pid_h1, int pid_h2, bool use_ML){
    
    // Create sPlot
    sPlot RF;
    RF.SetUp().SetOutDir(Form("%s/outsPlotBins/",brudir));
    
    //////////////////////////////////// Load variables
    RF.SetUp().LoadVariable("M1[0,2]");
    RF.SetUp().LoadVariable("M2[0,2]");
    RF.SetUp().LoadVariable("Mh[-100,100]");
    RF.SetUp().LoadVariable("z[-100,100]");
    RF.SetUp().LoadVariable("z1[-100,100]");
    RF.SetUp().LoadVariable("z2[-100,100]");
    RF.SetUp().LoadVariable("phi_h[-100,100]");
    RF.SetUp().LoadVariable("phi_R0[-100,100]");
    RF.SetUp().LoadVariable("th[-100,100]");
    RF.SetUp().LoadVariable("x[-100,100]");
    RF.SetUp().LoadVariable("Mx[-100,100]");
    RF.SetUp().LoadVariable("xF1[-100,100]");
    RF.SetUp().LoadVariable("xF2[-100,100]");
    RF.SetUp().LoadVariable("phi_h1[-100,100]");
    RF.SetUp().LoadVariable("phi_h2[-100,100]");
    RF.SetUp().SetIDBranchName("fgID");
    
    // From PID's determine several properties of the sWeight fitting
    std::string cvar = "";
    float u1,u2,u3; // mean vars
    float s1,s2,s3; // sigma vars
    if(pid_h1==111 && pid_h2!=111){
        cvar = "M1";
        u1=0.131; u2=0.129; u3=0.15;
        s1=0.01;  s2=0.0001; s3=0.1;
        if(use_ML){
            RF.SetUp().LoadVariable("p_11[0,1]");
            RF.SetUp().LoadVariable("p_12[0,1]");
            RF.SetUp().AddCut("p_11>0.9");
            RF.SetUp().AddCut("p_12>0.9");
        }else{
            RF.SetUp().LoadVariable("isGoodEventWithoutML[0,1]");
            RF.SetUp().AddCut("isGoodEventWithoutML==1");
        }
    }
    else if(pid_h1!=111 && pid_h2==111){
        cvar = "M2";
        u1=0.131; u2=0.129; u3=0.15;
        s1=0.01;  s2=0.0001; s3=0.1;
        if(use_ML){
            RF.SetUp().LoadVariable("p_21[0,1]");
            RF.SetUp().LoadVariable("p_22[0,1]");
            RF.SetUp().AddCut("p_21>0.9");
            RF.SetUp().AddCut("p_22>0.9");
        }else{
            RF.SetUp().LoadVariable("isGoodEventWithoutML[0,1]");
            RF.SetUp().AddCut("isGoodEventWithoutML==1");
        }
    }
    else if(pid_h1==111 && pid_h2==111){
        cvar = "M12";
        u1=0.262; u2=0.258; u3=0.3;
        s1=0.02;  s2=0.0002; s3=0.2;
        if(use_ML){
            RF.SetUp().LoadVariable("p_11[0,1]");
            RF.SetUp().LoadVariable("p_12[0,1]");
            RF.SetUp().AddCut("p_11>0.9");
            RF.SetUp().AddCut("p_12>0.9");
            RF.SetUp().LoadVariable("p_21[0,1]");
            RF.SetUp().LoadVariable("p_22[0,1]");
            RF.SetUp().AddCut("p_21>0.9");
            RF.SetUp().AddCut("p_22>0.9");
        }else{
            RF.SetUp().LoadVariable("isGoodEventWithoutML[0,1]");
            RF.SetUp().AddCut("isGoodEventWithoutML==1");
        }
    }
    
    // For loop over the cuts desired to apply them
    for(int j = 0 ; j < cut.vmin.size() ; j++){                     RF.SetUp().AddCut(Form("%s>%f",cut.var.at(j).c_str(),cut.vmin.at(j)).Data());
   RF.SetUp().AddCut(Form("%s<%f",cut.var.at(j).c_str(),cut.vmax.at(j)).Data());     
    }
    
    //////////////////////////////////// Gaussian
    RF.SetUp().FactoryPDF(Form("Gaussian::Signal( %s, mean[%f,%f,%f], sigma[%f,%f,%f] )",cvar.c_str(),u1,u2,u3,s1,s2,s3).Data());
    RF.SetUp().LoadSpeciesPDF("Signal",1);
  //////////////////////////////////// Background
    RF.SetUp().FactoryPDF(Form("Chebychev::BG( %s, {a0[-0.1,-1,1], a1[0.1,-1,1], a2[-0.1,-1,1], a3[-0.1,-1,1]} )",cvar.c_str()).Data());
    RF.SetUp().LoadSpeciesPDF("BG",1);
    
    // Load data
    RF.LoadData("dihadron",infile);
    
    //////////////////////////////////// Bins
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
    RF.DrawWeighted("Mgg","Signal");
    RF.DeleteWeightedTree();
}