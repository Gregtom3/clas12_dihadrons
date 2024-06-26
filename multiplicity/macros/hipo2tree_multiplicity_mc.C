#include "../../src/CutManager.C"
#include "../../src/CLAS12Analysis.C"
#include "../../src/TreeManager.C"
#include "../../src/HipoBankInterface.C"
#include "../../src/Constants.h"
#include "../../src/Structs.h"
#include "../../src/Kinematics.C"
#include "../../src/ParseBinYAML.C"
#include "../../src/ParseText.C"

double calculate_angle(part p1, part p2, bool use_true=false){
    TVector3 v1, v2;
    if(use_true){
        v1.SetXYZ(p1.truepx,p1.truepy,p1.truepz);
        v2.SetXYZ(p2.truepx,p2.truepy,p2.truepz);
    } else {
        v1.SetXYZ(p1.px,p1.py,p1.pz);
        v2.SetXYZ(p2.px,p2.py,p2.pz);
    }
    return v1.Angle(v2)*180/3.14159265;
}

bool checkDihadronCuts(int pid_h1, int pid_h2, double Mx, double z, double xF1, double xF2, int idx_e, std::vector<part> vec_parts, int i, int ii, int j, int jj, bool use_true=false) {
    
    bool passDihadron = true;
    
    // Minimum momentum cuts for hadrons
    if(pid_h1!=111 && (use_true ? vec_parts[i].truep<1.25 : vec_parts[i].p<1.25)) passDihadron=false;
    if(pid_h2!=111 && (use_true ? vec_parts[j].truep<1.25 : vec_parts[j].p<1.25)) passDihadron=false;
    
    // Missing mass cut for pi+pi- and pi+pi0 final states
    if((pid_h1==211 && pid_h2==-211) || (pid_h1==211 && pid_h2==111)){
        if(Mx<1.5) passDihadron=false;
    }
    
    // Maximum z to avoid exclusive events
    if(z>0.95) passDihadron=false;
    
    // CFR hadrons
    if(xF1<0 || xF2<0) passDihadron=false;
    
    // Minimum photon energy
    if(pid_h1==111 && ((use_true ? vec_parts[i].trueE<0.2 : vec_parts[i].E<0.2) || (use_true ? vec_parts[ii].trueE<0.2 : vec_parts[ii].E<0.2))) passDihadron=false;
    if(pid_h2==111 && ((use_true ? vec_parts[j].trueE<0.2 : vec_parts[j].E<0.2) || (use_true ? vec_parts[jj].trueE<0.2 : vec_parts[jj].E<0.2))) passDihadron=false;
    
    // Minimum angle between scattered e- and gamma, and FD cut on photons
    if(pid_h1==111){
        if(calculate_angle(vec_parts[idx_e],vec_parts[i],use_true)<8 || calculate_angle(vec_parts[idx_e],vec_parts[ii],use_true)<8) passDihadron=false;
        if((use_true ? vec_parts[i].truetheta*180/3.141592<5 : vec_parts[i].theta*180/3.141592<5) || (use_true ? vec_parts[i].truetheta*180/3.14159265>35 : vec_parts[i].theta*180/3.14159265>35)) passDihadron=false;
        if((use_true ? vec_parts[ii].truetheta*180/3.141592<5 : vec_parts[ii].theta*180/3.141592<5) || (use_true ? vec_parts[ii].truetheta*180/3.14159265>35 : vec_parts[ii].theta*180/3.14159265>35)) passDihadron=false;
    }
    if(pid_h2==111){
        if(calculate_angle(vec_parts[idx_e],vec_parts[j],use_true)<8 || calculate_angle(vec_parts[idx_e],vec_parts[jj],use_true)<8) passDihadron=false;
        if((use_true ? vec_parts[j].truetheta*180/3.141592<5 : vec_parts[j].theta*180/3.141592<5) || (use_true ? vec_parts[j].truetheta*180/3.14159265>35 : vec_parts[j].theta*180/3.14159265>35)) passDihadron=false;
        if((use_true ? vec_parts[jj].truetheta*180/3.141592<5 : vec_parts[jj].theta*180/3.141592<5) || (use_true ? vec_parts[jj].truetheta*180/3.14159265>35 : vec_parts[jj].theta*180/3.14159265>35)) passDihadron=false;
    }
    
    // FD cut on hadrons
    if(pid_h1!=111 && ((use_true ? vec_parts[i].truetheta*180/3.14159265<5 : vec_parts[i].theta*180/3.14159265<5) || (use_true ? vec_parts[i].truetheta*180/3.14159265>35 : vec_parts[i].theta*180/3.14159265>35))) passDihadron=false;
    if(pid_h2!=111 && ((use_true ? vec_parts[j].truetheta*180/3.14159265<5 : vec_parts[j].theta*180/3.14159265<5) || (use_true ? vec_parts[j].truetheta*180/3.14159265>35 : vec_parts[j].theta*180/3.14159265>35))) passDihadron=false;

    return passDihadron;
}


int hipo2tree_multiplicity_mc(
              const char * hipoFile = "/cache/clas12/rg-a/production/montecarlo/clasdis/fall2018/torus+1/v1/bkg50nA_10604MeV/50nA_OB_job_3313_0.hipo",
              const char * outputFile = "data.root",
              const double _electron_beam_energy = 10.6041,
              const int pid_h1=211,
              const int pid_h2=-211,
              const int maxEvents = 5000)
{
    
    
  // Open TTree and declare branches
  // -------------------------------------
  TFile *fOut = new TFile(outputFile,"RECREATE");
  TTree* tree = new TTree("dihadron_cuts","dihadron_cuts");
  EventTree *gen_tree = new EventTree("mc_particles");
  EventTree *reco_tree = new EventTree("reco_particles");
    
  int run,rec_passDihadron, MCmatch, MC_2h_match;
  double gen_x, gen_y, gen_z, gen_nu, gen_Q2, gen_W, gen_xF1, gen_xF2, gen_xF, gen_pT1, gen_pT2, gen_pTtot, gen_phi_h, gen_phi_R0, gen_phi_R1, gen_th, gen_E_e, gen_theta_e, gen_phi_e, gen_Mh, gen_z1, gen_z2,gen_Mx;
  double rec_x, rec_y, rec_z, rec_nu, rec_Q2, rec_W, rec_xF1, rec_xF2, rec_xF, rec_pT1, rec_pT2, rec_pTtot, rec_phi_h, rec_phi_R0, rec_phi_R1, rec_th, rec_E_e, rec_theta_e, rec_phi_e, rec_Mh, rec_z1, rec_z2,rec_Mx;

  tree->Branch("MCmatch", &MCmatch, "MCmatch/I");
  tree->Branch("MC_2h_match", &MC_2h_match, "MC_2h_match/I");
    // Generated variables
    tree->Branch("gen_x", &gen_x, "gen_x/D");
    tree->Branch("gen_y", &gen_y, "gen_y/D");
    tree->Branch("gen_z", &gen_z, "gen_z/D");
    tree->Branch("gen_nu", &gen_nu, "gen_nu/D");
    tree->Branch("gen_Q2", &gen_Q2, "gen_Q2/D");
    tree->Branch("gen_W", &gen_W, "gen_W/D");
    tree->Branch("gen_xF1", &gen_xF1, "gen_xF1/D");
    tree->Branch("gen_xF2", &gen_xF2, "gen_xF2/D");
    tree->Branch("gen_xF", &gen_xF, "gen_xF/D");
    tree->Branch("gen_pT1", &gen_pT1, "gen_pT1/D");
    tree->Branch("gen_pT2", &gen_pT2, "gen_pT2/D");
    tree->Branch("gen_pTtot", &gen_pTtot, "gen_pTtot/D");
    tree->Branch("gen_phi_h", &gen_phi_h, "gen_phi_h/D");
    tree->Branch("gen_phi_R0", &gen_phi_R0, "gen_phi_R0/D");
    tree->Branch("gen_phi_R1", &gen_phi_R1, "gen_phi_R1/D");
    tree->Branch("gen_th", &gen_th, "gen_th/D");
    tree->Branch("gen_E_e", &gen_E_e, "gen_E_e/D");
    tree->Branch("gen_theta_e", &gen_theta_e, "gen_theta_e/D");
    tree->Branch("gen_phi_e", &gen_phi_e, "gen_phi_e/D");
    tree->Branch("gen_Mh", &gen_Mh, "gen_Mh/D");
    tree->Branch("gen_Mx", &gen_Mx, "gen_Mx/D");
    tree->Branch("gen_z1", &gen_z1, "gen_z1/D");
    tree->Branch("gen_z2", &gen_z2, "gen_z2/D");
    
    // Reconstructed variables
    tree->Branch("rec_passDihadron",&rec_passDihadron,"rec_passDihadron/I");
    tree->Branch("rec_x", &rec_x, "rec_x/D");
    tree->Branch("rec_y", &rec_y, "rec_y/D");
    tree->Branch("rec_z", &rec_z, "rec_z/D");
    tree->Branch("rec_nu", &rec_nu, "rec_nu/D");
    tree->Branch("rec_Q2", &rec_Q2, "rec_Q2/D");
    tree->Branch("rec_W", &rec_W, "rec_W/D");
    tree->Branch("rec_xF1", &rec_xF1, "rec_xF1/D");
    tree->Branch("rec_xF2", &rec_xF2, "rec_xF2/D");
    tree->Branch("rec_xF", &rec_xF, "rec_xF/D");
    tree->Branch("rec_pT1", &rec_pT1, "rec_pT1/D");
    tree->Branch("rec_pT2", &rec_pT2, "rec_pT2/D");
    tree->Branch("rec_pTtot", &rec_pTtot, "rec_pTtot/D");
    tree->Branch("rec_phi_h", &rec_phi_h, "rec_phi_h/D");
    tree->Branch("rec_phi_R0", &rec_phi_R0, "rec_phi_R0/D");
    tree->Branch("rec_phi_R1", &rec_phi_R1, "rec_phi_R1/D");
    tree->Branch("rec_th", &rec_th, "rec_th/D");
    tree->Branch("rec_E_e", &rec_E_e, "rec_E_e/D");
    tree->Branch("rec_theta_e", &rec_theta_e, "rec_theta_e/D");
    tree->Branch("rec_phi_e", &rec_phi_e, "rec_phi_e/D");
    tree->Branch("rec_Mh", &rec_Mh, "rec_Mh/D");
    tree->Branch("rec_Mx", &rec_Mx, "rec_Mx/D");
    tree->Branch("rec_z1", &rec_z1, "rec_z1/D");
    tree->Branch("rec_z2", &rec_z2, "rec_z2/D");
    

  // Configure CLAS12 Reader and HipoChain
  // -------------------------------------
  clas12root::HipoChain _chain;
  clas12::clas12reader *_config_c12{nullptr};

  _chain.Add(hipoFile);
  _config_c12=_chain.GetC12Reader();

  // If not monte carlo, enforce QADB
  // -------------------------------------
  _config_c12->db()->turnOffQADB();


  // Configure PIDs for final state
  // -------------------------------------
  FS fs = get_FS(pid_h1,pid_h2);
  _config_c12->addAtLeastPid(11,1);     // At least 1 electron


  // Establish CLAS12 event parser
  // -------------------------------------
  auto &_c12=_chain.C12ref();
    
  // Create RCDB Connection
  // -------------------------------------
  clas12::clas12databases::SetRCDBRootConnection("/work/clas12/users/gmat/clas12/clas12_dihadrons/utils/rcdb.root"); 
  clas12::clas12databases db;
  
  // Add Analysis Objects
  // -------------------------------------
  CutManager _cm = CutManager();
  CLAS12Analysis clas12ana = CLAS12Analysis(_c12,_electron_beam_energy);
  clas12ana.set_run_config(_c12);
  Kinematics kin;
    
  // Add Analysis Structs
  // -------------------------------------  
  std::vector<part> vec_particles;
  std::vector<part> vec_mcparticles;
  EVENT_INFO event_info;
  EVENT event;

  // Declare TLorentzVectors
  // -------------------------------------
  TLorentzVector init_electron, electron,trueelectron,q,trueq,h1,trueh1,h2,trueh2,dihadron,truedihadron;
  init_electron.SetPxPyPzE(0,0,sqrt(_electron_beam_energy*_electron_beam_energy-Me*Me),_electron_beam_energy);
  TLorentzVector init_target(0,0,0,Mp);    
  int whileidx=0;
  int _ievent=0;
  int badAsym=0;
    
  while(_chain.Next()==true && (whileidx < maxEvents || maxEvents < 0)){
      
    if(whileidx%10000==0 && whileidx!=0){
      std::cout << whileidx << " events read | " << _ievent*100.0/whileidx << "% passed event selection | " << badAsym << " events skipped from QADB" << std::endl;
    }
      
    event_info = clas12ana.get_event_info(_c12);
    event_info.uID = whileidx;

          
    if(whileidx==0){
        init_electron.SetE(runBeamEnergy(event_info.run));
        init_electron.SetPz(sqrt(init_electron.E()*init_electron.E()-Me*Me));
    }
      
    whileidx++;
    // Set run specific information
    // -------------------------------------
    _cm.set_run(event_info.run);
    _cm.set_run_period(std::string(hipoFile));


    // *******************************************************************
    //     Reconstructed Particles
    //

    vec_particles = clas12ana.load_reco_particles(_c12);
    int idx_scattered_ele = clas12ana.find_reco_scattered_electron(vec_particles);
    if(idx_scattered_ele==-1)
      continue; // No scattered electron found
    vec_particles[idx_scattered_ele].is_scattered_electron=1;
    clas12ana.fill_reco_event_variables(event, vec_particles);
    if(event.y > 0.8)
      continue; // Maximum y cut
    vec_particles = _cm.filter_particles(vec_particles); // Apply Cuts
    if(clas12ana.reco_event_contains_scattered_electron(vec_particles)==false)
      continue; // REC::Particle must have found the scattered electron

    //
    //
    // *******************************************************************


    // *******************************************************************
    //     Monte Carlo Generated Particles
    //
      
    vec_mcparticles = clas12ana.load_mc_particles(_c12);
    int idx_scattered_ele_mc = clas12ana.find_mc_scattered_electron(vec_mcparticles);
    clas12ana.fill_mc_event_variables(event, vec_mcparticles);
    clas12ana.match_mc_to_reco(vec_particles, vec_mcparticles);

      
    // 
    //
    // *******************************************************************
      
    // Fill the ParticleTrees
    gen_tree->FillTree(vec_mcparticles,event,event_info);
    reco_tree->FillTree(vec_particles,event,event_info);
      
    // Set the event variables to 0
    gen_x = event.truex;
    gen_y = event.truey;
    gen_nu = event.truenu;
    gen_Q2 = event.trueQ2;
    gen_W = event.trueW;

    rec_x = event.x;
    rec_y = event.y;
    rec_nu = event.nu;
    rec_Q2 = event.Q2;
    rec_W = event.W;


    //Loop over all particles in the event to find electron
      
    electron.SetPxPyPzE(vec_particles[idx_scattered_ele].px,
                        vec_particles[idx_scattered_ele].py,
                        vec_particles[idx_scattered_ele].pz,
                        vec_particles[idx_scattered_ele].E);
      
    trueelectron.SetPxPyPzE(vec_mcparticles[idx_scattered_ele_mc].truepx,
                        vec_mcparticles[idx_scattered_ele_mc].truepy,
                        vec_mcparticles[idx_scattered_ele_mc].truepz,
                        vec_mcparticles[idx_scattered_ele_mc].trueE);

    rec_E_e = electron.E();
    rec_theta_e = electron.Theta();
    rec_phi_e = electron.Phi();
    
    gen_E_e = trueelectron.E();
    gen_theta_e = trueelectron.Theta();
    gen_phi_e = trueelectron.Phi();
      
    q=init_electron-electron;
    trueq=init_electron-trueelectron;
    
    // Get the indices of dihadrons from the MC Particles
    auto mc_dihadron_idxs = clas12ana.dihadron_idxs(pid_h1,pid_h2,vec_mcparticles);

    // Now loop over all MONTE CARLO dihadrons
    for(int a = 0 ; a < mc_dihadron_idxs.size() ; a++){
        std::vector<int> dihadron_idx = mc_dihadron_idxs.at(a);
        
        // Look at one specific MC dihadron pair
        // Set the event true dihadron kinematics
        // If those MC particles correspond to a reconstructed dihadron, set those kinematics
        
        clas12ana.fill_mc_reco_dihadron_variables(event, q, trueq, electron, trueelectron, vec_mcparticles, dihadron_idx, pid_h1, pid_h2);
            
        gen_z = event.truez;
        gen_xF1 = event.truexF1;
        gen_xF2 = event.truexF2;
        gen_xF = event.truexF;
        gen_pT1 = event.truepT1;
        gen_pT2 = event.truepT2;
        gen_pTtot = event.truepTtot;
        gen_phi_h = event.truephi_h;
        gen_phi_R0 = event.truephi_R0;
        gen_phi_R1 = event.truephi_R1;
        gen_th = event.trueth;
        gen_Mh = event.trueMh;
        gen_Mx = event.trueMx;
        gen_z1 = event.truez1;
        gen_z2 = event.truez2;
            
        rec_z = event.z;
        rec_xF1 = event.xF1;
        rec_xF2 = event.xF2;
        rec_xF = event.xF;
        rec_pT1 = event.pT1;
        rec_pT2 = event.pT2;
        rec_pTtot = event.pTtot;
        rec_phi_h = event.phi_h;
        rec_phi_R0 = event.phi_R0;
        rec_phi_R1 = event.phi_R1;
        rec_th = event.th;
        rec_Mh = event.Mh;
        rec_Mx = event.Mx;
        rec_z1 = event.z1;
        rec_z2 = event.z2;
            
        MCmatch = event.MCmatch;
        MC_2h_match = event.MC_2h_match;
        
        rec_passDihadron = checkDihadronCuts(pid_h1,pid_h2,rec_Mx,rec_z,rec_xF1,rec_xF2,idx_scattered_ele_mc, vec_mcparticles, event.i, event.ii, event.j, event.jj,false);
        
        // Before filling, make sure that the gen_Mx is greater than 1.5 to step away from exclusives
        // Missing mass cut for pi+pi- and pi+pi0 final states
        if((pid_h1==211 && pid_h2==-211) || (pid_h1==211 && pid_h2==111)){
            if(gen_Mx<1.5) continue;
        }
        
        tree->Fill();
        
    }
      
      
    // Get the indices of dihadrons from the RECONSTRUCTED Particles
    auto rec_dihadron_idxs = clas12ana.dihadron_idxs(pid_h1,pid_h2,vec_particles);

    // Now loop over all RECONSTRUCTED dihadrons
    for(int a = 0 ; a < rec_dihadron_idxs.size() ; a++){
        std::vector<int> dihadron_idx = rec_dihadron_idxs.at(a);
        
        // Look at one specific RECONSTRUCTED dihadron pair
        // If it has a Monte Carlo match --> SKIP (we already filled it in the previous for loop)
        // If it does not have a Monte Carlo match --> FILL 
        
        clas12ana.fill_mc_reco_dihadron_variables(event, q, trueq, electron, trueelectron, vec_particles, dihadron_idx, pid_h1, pid_h2);
            
        if(event.MC_2h_match==1)
            continue;
            
        gen_z = event.truez;
        gen_xF1 = event.truexF1;
        gen_xF2 = event.truexF2;
        gen_xF = event.truexF;
        gen_pT1 = event.truepT1;
        gen_pT2 = event.truepT2;
        gen_pTtot = event.truepTtot;
        gen_phi_h = event.truephi_h;
        gen_phi_R0 = event.truephi_R0;
        gen_phi_R1 = event.truephi_R1;
        gen_th = event.trueth;
        gen_Mh = event.trueMh;
        gen_Mx = event.trueMx;
        gen_z1 = event.truez1;
        gen_z2 = event.truez2;
            
        rec_z = event.z;
        rec_xF1 = event.xF1;
        rec_xF2 = event.xF2;
        rec_xF = event.xF;
        rec_pT1 = event.pT1;
        rec_pT2 = event.pT2;
        rec_pTtot = event.pTtot;
        rec_phi_h = event.phi_h;
        rec_phi_R0 = event.phi_R0;
        rec_phi_R1 = event.phi_R1;
        rec_th = event.th;
        rec_Mh = event.Mh;
        rec_Mx = event.Mx;
        rec_z1 = event.z1;
        rec_z2 = event.z2;
            
        MCmatch = event.MCmatch;
        MC_2h_match = event.MC_2h_match;
        
        rec_passDihadron = checkDihadronCuts(pid_h1,pid_h2,rec_Mx,rec_z,rec_xF1,rec_xF2,idx_scattered_ele_mc, vec_particles, event.i, event.ii, event.j, event.jj,false);
        
        tree->Fill();
        
    }
    
    _ievent++;
  }
  fOut->cd();
  gen_tree->Write();
  reco_tree->Write();
  tree->Write();
  fOut->Close(); 
  return 0;
}


