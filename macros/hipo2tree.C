#include "../src/CutManager.C"
#include "../src/CLAS12Analysis.C"
#include "../src/TreeManager.C"
#include "../src/HipoBankInterface.C"
#include "../src/Constants.h"
#include "../src/Structs.h"
#include "../src/Kinematics.C"
#include "../src/ParseBinYAML.C"
#include "../src/ParseText.C"


int hipo2tree(
	                   const char * hipoFile = "/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v1/dst/train/nSidis/nSidis_005036.hipo",
	      //const char * hipoFile = "/cache/clas12/rg-a/production/montecarlo/clasdis/fall2018/torus+1/v1/bkg50nA_10604MeV/50nA_OB_job_3313_0.hipo",
	      //const char * hipoFile = "/cache/clas12/rg-b/production/recon/spring2020/torus-1/pass1/v1/dst/train/sidisdvcs/sidisdvcs_011494.hipo",
              const char * outputFile = "hipo2tree.root",
              const double _electron_beam_energy = 10.6041,
              const int pid_h1=211,
              const int pid_h2=-211,
              const int maxEvents = 500000000,
              bool hipo_is_mc = false)
{



  // Create a TFile to save the data
  TFile* fOut = new TFile(outputFile, "RECREATE");

  // Create a TTree to store the data
  EventTree * tree = new EventTree("EventTree");
    
  // Configure CLAS12 Reader and HipoChain
  // -------------------------------------
  clas12root::HipoChain _chain;
  clas12::clas12reader *_config_c12{nullptr};

  _chain.Add(hipoFile);
  _config_c12=_chain.GetC12Reader();

  // If not monte carlo, enforce QADB
  // -------------------------------------
  bool do_QADB=(hipo_is_mc==false && std::string(hipoFile).find("/rg-c/")==std::string::npos);
  if(!do_QADB)
    _config_c12->db()->turnOffQADB();


  // Configure PIDs for final state
  // -------------------------------------
  FS fs = get_FS(pid_h1,pid_h2);
  _config_c12->addAtLeastPid(11,1);     // At least 1 electron
  if(fs.pid_h1!=0)    _config_c12->addAtLeastPid(fs.pid_h1,fs.num_h1);
  if(fs.pid_h2!=0)    _config_c12->addAtLeastPid(fs.pid_h2,fs.num_h2); // Doesn't run if duplicate final state

  // Establish CLAS12 event parser
  // -------------------------------------
  auto &_c12=_chain.C12ref();
  if(do_QADB)
    _c12->db()->qadb_requireOkForAsymmetry(true);  
    
  // Create RCDB Connection
  // -------------------------------------
  clas12::clas12databases::SetRCDBRootConnection("/work/clas12/users/gmat/clas12/clas12_dihadrons/utils/rcdb.root"); 
  clas12::clas12databases db;
  
  // Add Analysis Objects
  // -------------------------------------
  CutManager _cm = CutManager();
  CLAS12Analysis clas12ana = CLAS12Analysis(_c12,_electron_beam_energy);
  clas12ana.set_run_config(_c12);
    
  // Add Analysis Structs
  // -------------------------------------  
  std::vector<part> vec_particles;
  std::vector<part> vec_mcparticles;
  EVENT_INFO event_info;
  EVENT reco_event;
  EVENT   mc_event;

    
  int whileidx=0;
  int _ievent=0;
  int badAsym=0;
    
  while(_chain.Next()==true && (whileidx < maxEvents || maxEvents < 0)){
      
    if(whileidx%10000==0 && whileidx!=0){
      std::cout << whileidx << " events read | " << _ievent*100.0/whileidx << "% passed event selection | " << badAsym << " events skipped from QADB" << std::endl;
    }
      
    event_info = clas12ana.get_event_info(_c12);
    event_info.uID = whileidx;
    whileidx++;
    
    // Set run specific information
    // -------------------------------------
    _cm.set_run(event_info.run);
    _cm.set_run_period(std::string(hipoFile));
    
    
    // Skip events that are not ok for asymmetry analysis based on QADB
    if(do_QADB){
        if(!_c12->db()->qa()->isOkForAsymmetry(event_info.run,event_info.evnum)){
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
    reco_event = clas12ana.calc_reco_event_variables(vec_particles);
    if(reco_event.y > 0.8)
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
        mc_event = clas12ana.calc_mc_event_variables(vec_mcparticles);
        clas12ana.match_mc_to_reco(vec_particles, vec_mcparticles);
    }
      
    // 
    //
    // *******************************************************************
    tree->FillTree(vec_particles,reco_event,mc_event,event_info);

    _ievent++;
  }
  fOut->cd();
  tree->Write();
  fOut->Close(); 
  return 0;
}
