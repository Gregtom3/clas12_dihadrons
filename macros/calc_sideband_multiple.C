#include "../src/ParseBinYAML.C"
#include "../src/fitTools.C"
#include "../src/Constants.h"
#include "../src/ParseText.C"


void asym(const char *infile, std::string new_outdir, std::string version, double Mgg_min, double Mgg_max);

//
//
// Program: calc_sideband_multiple(const char *, std::string)
// Author:  Gregory Matousek
// Date:    7/17/2023
//
//


// Program to perform multiple sideband asymmetry extractions in different regions
// infile --> Input "merged_cuts" root file
// brudir --> Main output directory path

int calc_sideband_multiple(const char *infile = "/volatile/clas12/users/gmat/clas12analysis.sidis.data/clas12_dihadrons/projects/pipi0_paper_RGA_only/data/piminus_pi0/Fall2018_RGA_outbending_merged_cuts.root", 
                           std::string outdir  = "./bru"
                           )
{
    // Define a set of Mgg Ranges
    const std::vector<std::vector<double>> Mgg_pairs = {{0.00,0.05},
                                                        {0.05,0.08},
                                                        {0.08,0.11},
                                                        {0.11,0.14},
                                                        {0.14,0.17},
                                                        {0.17,0.20},
                                                        {0.2,0.25},
                                                        {0.25,0.3},
                                                        {0.3,0.35},
                                                        {0.35,0.4},
                                                        {0.4,0.475},
                                                        {0.475,0.55},
                                                        {0.55,0.65}};
    
    // Get the pid's from the input file
    int pid_h1=0;
    int pid_h2=0;
    std::string hadron_pair="";
    getPIDs(std::string(infile),pid_h1,pid_h2,hadron_pair);   
    // This code should only work for PiPlusPi0 and PiMinusPi0
    if(pid_h1==111 || pid_h2!=111){
        cout << "ERROR: Code only intended for PiPlusPi0 and PiMinusPi0...Aborting..." << endl;
        return -1;
    }
    
    // Pull the ..._<version>_merged_cuts.root from "infile"
    std::string version=getVersion2(infile);
    
    // Make the first subdirectory 
    outdir+="/"+version;
    gSystem->mkdir(TString(outdir)); // <project>/asym/<version>
    
    // Make the second subdirectory 
    outdir+="/precut";
    gSystem->mkdir(TString(outdir)); // <project>/asym/<version>/<cut>
    
    // Create a directory for the hadron pair
    outdir+="/"+hadron_pair; 
    gSystem->mkdir(TString(outdir)); // <project>/asym/<version>/<cut>/<hadron_pair>
    
    
    // Loop over the Mgg_pairs
    double Mgg_min = 0;
    double Mgg_max = 0;
    for(const auto Mgg_pair: Mgg_pairs){
        // Get the minimum and maximum Mgg from the pair
        Mgg_min = Mgg_pair.at(0);
        Mgg_max = Mgg_pair.at(1);
        // Create new subdirectory for this specific pair
        std::string new_outdir = outdir+"/"+Form("%f_%f",Mgg_min,Mgg_max);
        gSystem->mkdir(TString(new_outdir));
        // Create the sideband
        asym(infile, new_outdir, version, Mgg_min, Mgg_max);
        // Perform asymmetry calculation
        
    
    }

    return 0;
}




void asym(const char *infile, std::string new_outdir, std::string version, double Mgg_min, double Mgg_max){

    std::string hel_str="hel";  // changes if we inject the Monte Carlo

    FitManager FM;
    FM.SetUp().SetOutDir(Form("%s/outObsBins_sdbnd/",new_outdir.c_str()));

    //process_azi_FM(FM,version,hel_str);
    process_2h_FM(FM,version,hel_str);

    FM.SetUp().LoadVariable(Form("M2[%f,%f]",Mgg_min,Mgg_max));
    FM.SetUp().LoadVariable("isGoodEventWithoutML[0.5,1.5]");

    //////////////////////////////////// Bins
    //Load bin variables
//     for (int i = 0; i < binStruct.numDimensions; i++) {
//         std::string binName = binStruct.dimensionNames[i];
//         std::vector<double> binEdges = binStruct.binEdges[i];
//         int numBins = binEdges.size();
//         Double_t binEdgesArr[numBins];
//         for (int j = 0; j < numBins; j++) {
//             binEdgesArr[j] = binEdges[j];
//         }
//         FM.Bins().LoadBinVar(binName, numBins-1, binEdgesArr);
//     }

    FM.LoadData("dihadron_cuts",infile);
    Here::Go(&FM);
    
}
