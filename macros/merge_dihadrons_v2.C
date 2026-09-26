#include <dirent.h>
#include <sys/types.h>
#include "../src/Constants.h"

// Function prototypes
int getRunNumber(const std::string &filename, const std::string &version);
bool shouldKeepFile(const std::string &filename, int runNumber, const std::string &version);

std::vector<TString> getFilesInDir(const char* path, const std::string &version)
{
    std::vector<TString> files;
    DIR *dir;
    struct dirent *ent;

    if ((dir = opendir (path)) != NULL) {
        while ((ent = readdir (dir)) != NULL) {

            std::string filename = ent->d_name;
	    cout << filename << endl;
            // Skip short filenames and non-root files
            if (filename.length() <= 5 || filename.substr(filename.length() - 5) != ".root") continue;
	    if (filename.find("merged") != std::string::npos) continue;

            // Get the run number based on the filename and version
	    //            int runNumber = getRunNumber(filename, version);
	    try{
	      TFile * fIn = TFile::Open(TString(std::string(path) + std::string("/") + filename));
	      if(!fIn || fIn->IsZombie()){
		delete fIn;
		continue;
	      }
	    //	    TFile * fIn = new TFile(TString(std::string(path) + std::string("/") + filename));
	      if (!fIn->GetListOfKeys()->Contains("dihadron")) {
		std::cout << "Bad: TTree 'dihadron' not found in file " << filename << "...skipping file..." << std::endl;
		continue;
	      }
	      TTree * tIn = (TTree*)fIn->Get("dihadron");
	      int runNumber;
	      tIn->SetBranchAddress("run", &runNumber);
	      // Get the first event's "run" value
	      tIn->GetEntry(0);
	      // Determine if the file should be kept based on the run number and other conditions
	      if (shouldKeepFile(filename, runNumber, version)) {
                files.push_back(std::string(path) + std::string("/") + filename);
                cout << filename << endl;
	      }
	      fIn->Close();
	      delete fIn;
	    }
	    catch(const std::exception&e){
	      continue;
	    }

	}
	closedir(dir);
    }

    return files;
}

int getRunNumber(const std::string &filename, const std::string &version)
{
    int i1_values[] = {0, 11, 7, 14, 14};
    int i2_values[] = {4, 4, 4, 5, 4};

    for (int i = 0; i < 4; i++) {
        int i1 = i1_values[i];
        int i2 = i2_values[i];
        std::string numberString = filename.substr(i1, i2);
        try {
            return std::stoi(numberString);
        } catch (...) {
            // std::stoi failed, try the next i1 i2 combination
        }
    }
    
    return -1;
}

bool shouldKeepFile(const std::string &filename, int runNumber, const std::string &version)
{
    if (filename.find("merged") != std::string::npos) return false;

    if (filename.find("MC_RGA") != std::string::npos && version.find("MC_RGA") == std::string::npos) return false;

    if (filename.find("MC_RGB") != std::string::npos && version.find("MC_RGB") == std::string::npos) return false;

    if (filename.find("nSidis_RGA") != std::string::npos && (version.find("2018_RGA") == std::string::npos && version.find("2019_RGA") == std::string::npos)) return false;

    if (filename.find("sidisdvcs_RGC") != std::string::npos && (version.find("Data_RGC") == std::string::npos)) return false;

    if (filename.find("sidisdvcs_RGB") != std::string::npos && (version.find("2019_RGB") == std::string::npos) && (version.find("2020_RGB") == std::string::npos)) return false;

    if (filename.find("MC_RGC") != std::string::npos && (version.find("MC_RGC") == std::string::npos)) return false;
    

    string runPeriodFromConstants = runPeriod(runNumber);
    if (version == runPeriodFromConstants) return true;
    else if (version == "MC_RGB_inbending" && (runPeriodFromConstants=="MC_RGA_inbending")) return true;
    else if (version == "MC_RGB_outbending" && (runPeriodFromConstants=="MC_RGA_outbending")) return true;
    else if (version == "Fall2018Spring2019_RGA_inbending" && (runPeriodFromConstants=="Fall2018_RGA_inbending" || runPeriodFromConstants=="Spring2019_RGA_inbending")) return true;
    else if (version == "Fall2018Spring2019_RGA_inbendingoutbending" && (runPeriodFromConstants=="Fall2018_RGA_inbending"||runPeriodFromConstants=="Fall2018_RGA_outbending"||runPeriodFromConstants=="Spring2019_RGA_inbending")) return true;
    
    return false;
}



int merge_dihadrons_v2(
	   const char * rootdir = "/volatile/clas12/users/gmat/clas12analysis.sidis.data/clas12_dihadrons/projects/ana_v0/data/piplus_pi0",
	   const char * version       = "Fall2018Spring2019_RGA_inbending"
	   )

{
  string outfile = string(rootdir)+"/"+string(version)+"_merged.root";
  // Get files to merge
  auto files = getFilesInDir(rootdir,string(version));
    
    
  //------------------------------------------
  // MERGE THE DIHADRON TTREE
  // -----------------------------------------

  // Create large TChain
  TChain *chain = new TChain("dihadron");
  for(int i = 0 ; i < files.size() ; i++){
    cout << i+1 << " of " << files.size() << endl;
    chain->Add(files.at(i));
  }

  // Merge the TTrees into the outfile
  chain->Merge(outfile.c_str());
    
    
  //------------------------------------------
  // MERGE THE DIHADRON_CUTS TTREE
  // This TTree contains less events because the default cuts are placed (including ML)
  // -----------------------------------------
  string outfile_cuts = string(rootdir)+"/"+string(version)+"_merged_cuts.root";
  // Create large TChain
  TChain *chain_cuts = new TChain("dihadron_cuts");
  for(int i = 0 ; i < files.size() ; i++){
    cout << i+1 << " of " << files.size() << endl;
    chain_cuts->Add(files.at(i));
  }

  // Merge the TTrees into the outfile
  chain_cuts->Merge(outfile_cuts.c_str());

  //------------------------------------------
  // MERGE THE DIHADRON_EXCLUSIVE_CUTS TTREE
  // This TTree contains less events because the default cuts are placed (including ML)
  // -----------------------------------------
  string outfile_exclusive_cuts = string(rootdir)+"/"+string(version)+"_merged_exclusive_cuts.root";
  // Create large TChain
  TChain *chain_exclusive_cuts = new TChain("dihadron_exclusive_cuts");
  for(int i = 0 ; i < files.size() ; i++){
    cout << i+1 << " of " << files.size() << endl;
    chain_exclusive_cuts->Add(files.at(i));
  }

  // Merge the TTrees into the outfile
  chain_exclusive_cuts->Merge(outfile_exclusive_cuts.c_str());

        
  return 0;
} 
