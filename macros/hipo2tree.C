#include "../src/CutManager.C"
#include "../src/CLAS12Ana.C"
#include "../src/TreeManager.C"
#include "../src/HipoBankInterface.C"
#include "../src/Constants.h"
#include "../src/Structs.h"
#include "../src/Kinematics.C"
#include "../src/ParseBinYAML.C"
#include "../src/ParseText.C"
#include "QADB.h"
using namespace QA;


int hipo2tree(
	      const char * hipoFile = "/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass2/main/train/nSidis/nSidis_005038.hipo",
	      //const char * hipoFile = "/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v1/dst/train/nSidis/nSidis_00503*.hipo",
	      //const char * hipoFile = "/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass1/v1/dst/train/nSidis/nSidis_005455.hipo",
	      //const char * hipoFile = "/cache/clas12/rg-a/production/montecarlo/clasdis/fall2018/torus+1/v1/bkg50nA_10604MeV/50nA_OB_job_3313_0.hipo",
	      //const char * hipoFile = "/cache/clas12/rg-b/production/recon/spring2020/torus-1/pass1/v1/dst/train/sidisdvcs/sidisdvcs_011494.hipo",
	      //const char * hipoFile = "/cache/hallb/scratch/rg-c/dst/train/sidisdvcs/sidisdvcs*.hipo",
	      //	      const char * hipoFile = "/work/cebaf24gev/sidis/reconstructed/polarized-plus-10.5GeV-proton/hipo/0051.hipo",
              const char * outputFile = "nSidis_005038.root",
              const double _electron_beam_energy = 10.6041,
              const int pid_h1=211,
              const int pid_h2=111,
              const int maxEvents = 5000000000,
              bool hipo_is_mc = false)
{

  // Create a TFile to save the data
  TFile* fOut = new TFile(outputFile, "RECREATE");

  // Create a TTree to store the data
  EventTree * tree = new EventTree("EventTree");
    
  // Configure CLAS12 Reader and HipoChain
  // -------------------------------------
  std::cout << "Configuring CLAS12 Reader and HipoChain." << std::endl;
  clas12root::HipoChain _chain;
  // clas12::clas12reader *_config_c12{nullptr};

  _chain.Add(hipoFile);
  auto _config_c12=_chain.GetC12Reader();

  std::cout << "HipoChain initialized with " << _chain.GetNFiles() << " files." << std::endl;
  // If not monte carlo, enforce QADB
  // -------------------------------------
  bool do_QADB=(hipo_is_mc==false && std::string(hipoFile).find("/rg-c/")==std::string::npos);
  
  std::cout << "Enabling QADB." << std::endl;
  QADB * qa = new QADB("pass2");
  qa->CheckForDefect("TotalOutlier");     // these choices match the criteria of `OkForAsymmetry`
  qa->CheckForDefect("TerminalOutlier");
  qa->CheckForDefect("MarginalOutlier");
  qa->CheckForDefect("SectorLoss");
  qa->CheckForDefect("Misc");
  for(int run : { // list of runs with `Misc` defect that are allowed by `OkForAsymmetry`
    5046, 5047, 5051, 5128, 5129, 5130, 5158, 5159,
    5160, 5163, 5165, 5166, 5167, 5168, 5169, 5180,
    5181, 5182, 5183, 5400, 5448, 5495, 5496, 5505,
    5567, 5610, 5617, 5621, 5623, 6736, 6737, 6738,
    6739, 6740, 6741, 6742, 6743, 6744, 6746, 6747,
    6748, 6749, 6750, 6751, 6753, 6754, 6755, 6756,
    6757})
  qa->AllowMiscBit(run);
  

  // Configure PIDs for final state
  // -------------------------------------
  FS fs = get_FS(pid_h1,pid_h2);
  _config_c12->addAtLeastPid(11,1);     // At least 1 electron
  if(fs.pid_h1!=0)    _config_c12->addAtLeastPid(fs.pid_h1,fs.num_h1);
  if(fs.pid_h2!=0)    _config_c12->addAtLeastPid(fs.pid_h2,fs.num_h2); // Doesn't run if duplicate final state

  // Establish CLAS12 event parser
  // -------------------------------------
  auto &_c12=_chain.C12ref();
 
  // Create RCDB Connection
  // -------------------------------------
  std::cout << "Setting RCDB root connection." << std::endl;
  clas12::clas12databases::SetRCDBRootConnection("/work/clas12/users/gmat/clas12/clas12_dihadrons/utils/rcdb.root"); 
  clas12::clas12databases db;
  
  // Add Analysis Objects
  // -------------------------------------
  CutManager _cm = CutManager();
  CLAS12Ana clas12ana = CLAS12Ana(_c12,_electron_beam_energy);

  clas12ana.set_run_config(_c12);
  
  // Add Analysis Structs
  // -------------------------------------  
  std::vector<part> vec_particles;
  std::vector<part> vec_mcparticles;
  EVENT_INFO event_info;
  EVENT event;

    
  int whileidx=0;
  int _ievent=0;
  int badAsym=0;

  while(_chain.Next()==true && (whileidx < maxEvents || maxEvents < 0)){

    if(whileidx%10000==0 && whileidx!=0){
      std::cout << whileidx << " events read | " << _ievent*100.0/whileidx << "% passed event selection | " << badAsym << " events skipped from QADB" << std::endl;
    }
      
    clas12ana.get_event_info(_c12,event_info);
    event_info.uID = whileidx;
    whileidx++;
    
    // Set run specific information
    // -------------------------------------
    _cm.set_run(event_info.run);
    _cm.set_run_period(std::string(hipoFile));
    
    // Skip events that are not ok for asymmetry analysis based on QADB
    if(do_QADB){
        if(!qa->Pass(event_info.run,event_info.evnum)) {
            badAsym++;
            continue;
          }
        }

    // Skip helicity==0 events
    // -------------------------------------
    if(!hipo_is_mc && event_info.hel==0)
        continue;
    
    // *******************************************************************
    //     Reconstructed Particles
    //

    vec_particles = clas12ana.load_reco_particles(_c12);
    int idx_scattered_ele = clas12ana.find_reco_scattered_electron(vec_particles);
    if(idx_scattered_ele==-1)
      continue; // No scattered electron found
    vec_particles[idx_scattered_ele].is_scattered_electron=1;
    clas12ana.fill_reco_event_variables(event, vec_particles);
    if(event.y > 0.8 || event.Q2 < 1)
      continue; // Maximum y cut
    vec_particles = _cm.filter_particles(vec_particles); // Apply Cuts
    if(clas12ana.reco_event_contains_final_state(vec_particles,fs)==false)
      continue; // Missing final state particles needed for event

    //
    //
    // *******************************************************************


    // *******************************************************************
    //     Monte Carlo Generated Particles
    //
      
    if(hipo_is_mc){
        vec_mcparticles = clas12ana.load_mc_particles(_c12);
        clas12ana.fill_mc_event_variables(event,vec_mcparticles);
        clas12ana.match_mc_to_reco(vec_particles, vec_mcparticles);
    }
      
    // 
    //
    // *******************************************************************
    tree->FillTree(vec_particles,event,event_info);

    _ievent++;
  }
  fOut->cd();
  tree->Write();
  fOut->Close(); 
  return 0;
}
