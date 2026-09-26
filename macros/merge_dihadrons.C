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
    else if (version == "Data_RGC" && (runPeriodFromConstants=="Fall2022_RGC"||runPeriodFromConstants=="Summer2022_RGC")) return true;
    return false;
}



int merge_dihadrons(
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

  // Reopen the outfile
  TFile *F_cuts = new TFile(outfile_cuts.c_str(),"UPDATE");

  // Load in the merged TTree
  TTree *t_cuts = (TTree*)F_cuts->Get("dihadron_cuts");

  // Create ID branch for brufit
  double fgID_cuts=0;
  TBranch *bfgId_cuts = t_cuts->Branch("fggID",&fgID_cuts,"fggID/D"); // create new branch
  const int N_cuts = t_cuts->GetEntries();
  for(int i = 0; i < N_cuts; ++i){
    t_cuts->GetEntry(i); // load in all TBranches
    bfgId_cuts->Fill(); // Fill only the new TBranch
    fgID_cuts+=1;
  }  
  
  // Print final TTree
  t_cuts->Print();

  // Overwrite merged TTree without the fgID branch
  t_cuts->Write(0,TObject::kOverwrite);

  // Close TFile
  F_cuts->Close();
    
    
    //------------------------------------------
    // MERGE THE DIHADRON_CUTS_NO_PMIN TTREE
    //------------------------------------------
    string outfile_cuts_noPmin = string(rootdir) + "/" + string(version) + "_merged_cuts_noPmin.root";
    TChain *chain_cuts_noPmin = new TChain("dihadron_cuts_noPmin");
    for (int i = 0; i < files.size(); ++i) {
        cout << i + 1 << " of " << files.size() << endl;
        chain_cuts_noPmin->Add(files.at(i));
    }
    chain_cuts_noPmin->Merge(outfile_cuts_noPmin.c_str());
    
    TFile *F_cuts_noPmin = new TFile(outfile_cuts_noPmin.c_str(), "UPDATE");
    TTree *t_cuts_noPmin = (TTree *)F_cuts_noPmin->Get("dihadron_cuts_noPmin");
    
    double fgID_cuts_noPmin = 0;
    TBranch *bfgId_cuts_noPmin = t_cuts_noPmin->Branch("fggID", &fgID_cuts_noPmin, "fggID/D");
    const int N_cuts_noPmin = t_cuts_noPmin->GetEntries();
    for (int i = 0; i < N_cuts_noPmin; ++i) {
        t_cuts_noPmin->GetEntry(i);
        bfgId_cuts_noPmin->Fill();
        fgID_cuts_noPmin += 1;
    }
    t_cuts_noPmin->Write(0, TObject::kOverwrite);
    F_cuts_noPmin->Close();  
    
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

  // Reopen the outfile
  TFile *F_exclusive_cuts = new TFile(outfile_exclusive_cuts.c_str(),"UPDATE");

  // Load in the merged TTree
  TTree *t_exclusive_cuts = (TTree*)F_exclusive_cuts->Get("dihadron_exclusive_cuts");

  // Create ID branch for brufit
  double fgID_exclusive_cuts=0;
  TBranch *bfgId_exclusive_cuts = t_exclusive_cuts->Branch("fggID",&fgID_exclusive_cuts,"fggID/D"); // create new branch
  const int N_exclusive_cuts = t_exclusive_cuts->GetEntries();
  for(int i = 0; i < N_exclusive_cuts; ++i){
    t_exclusive_cuts->GetEntry(i); // load in all TBranches
    bfgId_exclusive_cuts->Fill(); // Fill only the new TBranch
    fgID_exclusive_cuts+=1;
  }  
  
  // Print final TTree
  t_exclusive_cuts->Print();

  // Overwrite merged TTree without the fgID branch
  t_exclusive_cuts->Write(0,TObject::kOverwrite);

  // Close TFile
  F_exclusive_cuts->Close();


    // ————————————————
    // merge dihadron_legacy_cuts
    // ————————————————
    string outfile_legacy = string(rootdir) + "/" + string(version) + "_merged_legacy_cuts.root";
    TChain *chain_legacy = new TChain("dihadron_legacy_cuts");
    for (auto &fn : files) chain_legacy->Add(fn);
    chain_legacy->Merge(outfile_legacy.c_str());

    TFile *F_legacy = new TFile(outfile_legacy.c_str(),"UPDATE");
    TTree *t_legacy = (TTree*)F_legacy->Get("dihadron_legacy_cuts");
    double fgID_legacy = 0;
    TBranch *bfgId_legacy = t_legacy->Branch("fggID", &fgID_legacy, "fggID/D");
    const int N_legacy = t_legacy->GetEntries();
    for (int i = 0; i < N_legacy; ++i) {
        t_legacy->GetEntry(i);
        bfgId_legacy->Fill();
        fgID_legacy += 1;
    }
    t_legacy->Write(0, TObject::kOverwrite);
    F_legacy->Close();

  return 0;
} 
