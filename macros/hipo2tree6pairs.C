#include "../src/CutManager.C"
#include "../src/CLAS12Ana.C"
#include "../src/TreeManager.C"
#include "../src/HipoBankInterface.C"
#include "../src/Constants.h"
#include "../src/Structs.h"
#include "../src/Kinematics.C"
#include "../src/ParseBinYAML.C"
#include "../src/ParseText.C"


int hipo2tree6pairs(
	  const char * hipoFile = "/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass2/train/nSidis/nSidis_005555.hipo",
          const char *outputPath = ".",
          const char *fileName = "nSidis_RGA_5555.root",
          const double _electron_beam_energy = 10.6041,
          const int maxEvents = 5000,
          bool hipo_is_mc = false)
{
  std::vector<std::vector<int>> pids = {{211,211},
                                   {211,111},
                                   {211,-211},
                                   {111,111},
                                   {-211,111},
                                   {-211,-211}};
  
  std::vector<std::string> pion_pairs = {"piplus_piplus",
                              "piplus_pi0",
                              "piplus_piminus",
                              "pi0_pi0",
                              "piminus_pi0",
                              "piminus_piminus"};
  std::vector<std::string> outputFiles;
  std::vector<TFile*> fOuts;
  std::vector<EventTree*> trees;
  std::vector<FS> fss;
  for(int i = 0 ; i < 6 ; i++){
      outputFiles.push_back(std::string(outputPath)+std::string("/")+pion_pairs.at(i)+std::string("/")+std::string(fileName));
      fOuts.push_back(new TFile(outputFiles.at(i).c_str(),"RECREATE"));
      trees.push_back(new EventTree("EventTree"));
      fss.push_back(get_FS(pids.at(i).at(0),pids.at(i).at(1)));
  }

  // Configure CLAS12 Reader and HipoChain
  // -------------------------------------
  clas12root::HipoChain _chain;
  clas12::clas12reader *_config_c12{nullptr};

  _chain.Add(hipoFile);
  _config_c12=_chain.GetC12Reader();

  // If not monte carlo, enforce QADB
  // -------------------------------------
  bool do_QADB=(hipo_is_mc==false && std::string(hipoFile).find("/rg-c/")==std::string::npos);
  // Default pass
  std::string pass = "pass1";

  // Check if the path contains "/pass2/"
  if (std::string(hipoFile).find("/pass2/") != std::string::npos) {
    pass = "pass2";
  }
  if(pass=="pass2"){
    do_QADB=false;
  }
  if(!do_QADB)
    _config_c12->db()->turnOffQADB();

  
  _config_c12->addAtLeastPid(11,1);     // At least 1 electron
  if (do_QADB) {
    _config_c12->applyQA(pass);
    _config_c12->db()->qadb_addQARequirement("OkForAsymmetry");
  }

  // Establish CLAS12 event parser
  // -------------------------------------
  auto &_c12=_chain.C12ref();

  //if(do_QADB)
  //  _c12->db()->qadb_requireOkForAsymmetry(true);  

  // Create RCDB Connection
  // -------------------------------------
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
    // if(do_QADB){
    //     if(!_c12->db()->qa()->isOkForAsymmetry(event_info.run,event_info.evnum)){
    //         badAsym++;
    //         continue;
    //     }
    // }

    // Skip helicity==0 events
    // -------------------------------------
    if(!hipo_is_mc && event_info.hel==0){
        continue;
    }
    
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
    for(int j = 0; j < 6; j++){
        if(clas12ana.reco_event_contains_final_state(vec_particles,fss.at(j))==false)
          continue; // Missing final state particles needed for event
        trees.at(j)->FillTree(vec_particles,event,event_info);
    }
      
    
    _ievent++;
  }
  for(int j = 0; j<6; j++){
      fOuts.at(j)->cd();
      trees.at(j)->Write();
      cout << "Writing " << outputFiles.at(j).c_str() << " with " << trees.at(j)->GetEntries() << " entries" << endl;
      fOuts.at(j)->Close();

  }
//   fOut->cd();
//   tree->Write();
//   fOut->Close(); 
  return 0;
}
