#include <dirent.h>
#include <sys/types.h>
#include "../src/Constants.h"
// std::vector<TString> getFilesInDir(const char* path, string version)
// {

//   std::vector<TString> files;
//   DIR *dir;
//   struct dirent *ent;
//   int nmin = 0;
//   int nmax = 0;
//   int i1 = 0;
//   int i2 = 4;
//   if(version=="Fall2018_RGA_inbending"){
//     nmin = 5032;
//     nmax = 5332;
//     i1 = 11;
//   }
//   else if(version=="Fall2018_RGA_outbending"){
//     nmin = 5333;
//     nmax = 5666;
//     i1 = 11;
//   }
//   else if(version=="Spring2019_RGA_inbending"){
//     nmin = 6616;
//     nmax = 6783;
//     i1 = 11;
//   }
//   else if(version=="MC_RGA_inbending"){
//     nmin = 3051;
//     nmax = 3304;
//     i1 = 7;
//   }
//   else if(version=="MC_RGA_outbending"){
//     nmin = 3313;
//     nmax = 3327;
//     i1 = 7;
//   }
//   else if(version=="MC_RGC"){
//     nmin = 0;
//     nmax = 1000;
//     i1 = 7;
//   }
//   else if(version=="Data_RGC"){
//     nmin = 0;
//     nmax = 20000;
//     i1 = 14;
//     i2 = 5;
//   }

//   if ((dir = opendir (path)) != NULL) {
//     while ((ent = readdir (dir)) != NULL) {
//       std::string filename = ent->d_name;
//       cout << filename << endl;
//       if (filename.length()<=5) continue;
//       if (filename.substr(filename.length() - 5) != ".root") continue;
//       if (filename.find("merged") != std::string::npos) continue;
//       if (filename.find("MC_RGA") != std::string::npos && version.find("MC_RGA") == std::string::npos)
//             continue;
//       if (filename.find("nSidis_RGA") != std::string::npos && (version.find("2018_RGA")==std::string::npos && version.find("2019_RGA")==std::string::npos))
//           continue;
//       if (filename.find("sidisdvcs_RGC") != std::string::npos && (version.find("Data_RGC")==std::string::npos))
//           continue;
//       if (filename.find("MC_RGC") != std::string::npos && (version.find("MC_RGC")==std::string::npos))
//           continue;
//       std::string numberString = filename.substr(i1, i2);
//       int number = std::stoi(numberString);
//       if (number <= nmax && number >= nmin) {
//           files.push_back(string(path)+string("/") + filename);
//       }
//     }
//     closedir (dir);
//   } 
  
//   return files;
// }

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

            // Skip short filenames and non-root files
            if (filename.length() <= 5 || filename.substr(filename.length() - 5) != ".root") continue;

            // Get the run number based on the filename and version
            int runNumber = getRunNumber(filename, version);

            // Determine if the file should be kept based on the run number and other conditions
            if (shouldKeepFile(filename, runNumber, version)) {
                files.push_back(std::string(path) + std::string("/") + filename);
                cout << filename << endl;
            }
        }
        closedir(dir);
    }

    return files;
}

int getRunNumber(const std::string &filename, const std::string &version)
{
    int i1_values[] = {0, 11, 7, 14};
    int i2_values[] = {4, 4, 4, 5};

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

    if (filename.find("nSidis_RGA") != std::string::npos && (version.find("2018_RGA") == std::string::npos && version.find("2019_RGA") == std::string::npos)) return false;

    if (filename.find("sidisdvcs_RGC") != std::string::npos && (version.find("Data_RGC") == std::string::npos)) return false;

    if (filename.find("MC_RGC") != std::string::npos && (version.find("MC_RGC") == std::string::npos)) return false;

    string runPeriodFromConstants = runPeriod(runNumber);
    
    if (version == runPeriodFromConstants) return true;
    else if (version == "Fall2018Spring2019_RGA_inbending" && (runPeriodFromConstants=="Fall2018_RGA_inbending" || runPeriodFromConstants=="Spring2019_RGA_inbending")) return true;
    else if (version == "Fall2018Spring2019_RGA_inbendingoutbending" && (runPeriodFromConstants=="Fall2018_RGA_inbending"||runPeriodFromConstants=="Fall2018_RGA_outbending"||runPeriodFromConstants=="Spring2019_RGA_inbending")) return true;
    
    return false;
}



int merge_dihadrons(
	   const char * rootdir = "/volatile/clas12/users/gmat/clas12analysis.sidis.data/clas12_dihadrons/projects/ana_v0/data/pi0_pi0",
	   const char * version       = "Fall2018Spring2019_RGA_inbending"
	   )

{
  string outfile = string(rootdir)+"/"+string(version)+"_merged.root";
  // Get files to merge
  auto files = getFilesInDir(rootdir,string(version));
  // Create large TChain
  TChain *chain = new TChain("dihadron");
  for(int i = 0 ; i < files.size() ; i++){
    cout << i+1 << " of " << files.size() << endl;
    chain->Add(files.at(i));
  }

  // Merge the TTrees into the outfile
  chain->Merge(outfile.c_str());

  // Reopen the outfile
  TFile *F = new TFile(outfile.c_str(),"UPDATE");

  // Load in the merged TTree
  TTree *t = (TTree*)F->Get("dihadron");

  // Create ID branch for brufit
  double fgID=0;
  TBranch *bfgId = t->Branch("fggID",&fgID,"fggID/D"); // create new branch
  const int N = t->GetEntries();
  for(int i = 0; i < N; ++i){
    t->GetEntry(i); // load in all TBranches
    bfgId->Fill(); // Fill only the new TBranch
    fgID+=1;
  }  
  
  // Print final TTree
  t->Print();

  // Overwrite merged TTree without the fgID branch
  t->Write(0,TObject::kOverwrite);

  // Close TFile
  F->Close();
        
  return 0;
} 
