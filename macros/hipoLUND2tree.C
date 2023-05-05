#include "../src/Constants.h"
#include "../src/Structs.h"
#include "../src/Kinematics.C"

void generate_combinations(std::vector<int>& input, int num, int start_idx, std::vector<int>& curr_combination, std::vector<std::vector<int>>& result) {
    if (num == 0) {
        result.push_back(curr_combination);
        return;
    }
    for (int i = start_idx; i <= input.size() - num; i++) {
        curr_combination.push_back(input[i]);
        generate_combinations(input, num - 1, i + 1, curr_combination, result);
        curr_combination.pop_back();
    }
}

std::vector<std::vector<int>> unique_combinations(std::vector<int> input, int num) {
    std::vector<std::vector<int>> result;
    std::vector<int> curr_combination;
    std::sort(input.begin(), input.end());

    generate_combinations(input, num, 0, curr_combination, result);
    return result;
}

std::vector<std::vector<int>> remove_duplicates(std::vector<std::vector<int>> input) {
    std::vector<std::vector<int>> output;

    // Store the original indices and sort each inner vector based on the values
    std::vector<std::vector<size_t>> indices(input.size());
    for (size_t i = 0; i < input.size(); ++i) {
        indices[i].resize(input[i].size());
        std::iota(indices[i].begin(), indices[i].end(), 0);

        std::sort(indices[i].begin(), indices[i].end(),
                  [&](size_t a, size_t b) { return input[i][a] < input[i][b]; });
        std::sort(input[i].begin(), input[i].end());
    }

    // Sort and remove duplicates from the outer vector
    std::sort(input.begin(), input.end());
    input.erase(std::unique(input.begin(), input.end()), input.end());

    // Restore the original order of the inner vectors
    for (size_t i = 0; i < input.size(); ++i) {
        std::vector<int> temp(input[i].size());
        for (size_t j = 0; j < input[i].size(); ++j) {
            temp[indices[i][j]] = input[i][j];
        }
        input[i] = temp;
    }

    return input;
}

// This program parses through a Monte Carlo hipo file, ONLY reading the MC::Lund banks
// For each event, the scattered electron is found
// Then, the dihadrons are formed
// The important event variables are then saved to a TTree
int hipoLUND2tree(const char * input_file = "/cache/clas12/rg-a/production/montecarlo/clasdis/fall2018/torus-1/v1/bkg45nA_10604MeV/45nA_job_3051_0.hipo",
                  const char * outputFile = "hipoLUND2tree.root",
                  const int pid_h1=211,
                  const int pid_h2=-211,
                  const int maxEvents=50000000){
    
  // Open TTree and declare branches
  // -------------------------------------
  TFile *fOut = new TFile(outputFile,"RECREATE");
  TTree* tree = new TTree("dihadron","dihadron");
  double _electron_beam_energy;
  int run;
  double x,y,z,nu,s,Q2,W,xF1,xF2,xF,P1,P2,z1,z2,Mh,Mx,pT1,pT2,pTtot,phi_h,phi_R0,phi_R1,th;
    

  tree->Branch("x", &x, "x/D");
  tree->Branch("y", &y, "y/D");
  tree->Branch("z", &z, "z/D");
  tree->Branch("nu", &nu, "nu/D");
  tree->Branch("s", &s, "s/D");
  tree->Branch("Q2", &Q2, "Q2/D");
  tree->Branch("W", &W, "W/D");
  tree->Branch("xF1", &xF1, "xF1/D");
  tree->Branch("xF2", &xF2, "xF2/D");
  tree->Branch("xF", &xF, "xF/D");
  tree->Branch("z1", &z1, "z1/D");
  tree->Branch("z2", &z2, "z2/D");
  tree->Branch("Mh", &Mh, "Mh/D");
  tree->Branch("Mx", &Mx, "Mx/D");
  tree->Branch("pT1", &pT1, "pT1/D");
  tree->Branch("pT2", &pT2, "pT2/D");
  tree->Branch("pTtot", &pTtot, "pTtot/D");
  tree->Branch("phi_h", &phi_h, "phi_h/D");
  tree->Branch("phi_R0", &phi_R0, "phi_R0/D");
  tree->Branch("phi_R1", &phi_R1, "phi_R1/D");
  tree->Branch("th", &th, "th/D");
    
  // Configure CLAS12 Reader and HipoChain
  // -------------------------------------
  clas12root::HipoChain _chain;
  clas12::clas12reader *_config_c12{nullptr};

  _chain.Add(input_file);
  _config_c12=_chain.GetC12Reader();  
    
  // Add RUN::config bank
  // -------------------------------------
  int _idx_RUNconfig = _config_c12->addBank("RUN::config");
  int _irun = _config_c12->getBankOrder(_idx_RUNconfig,"run");
  
  // Establish CLAS12 event parser
  // -------------------------------------
  auto &_c12=_chain.C12ref();
    
  // Configure PIDs for final state
  // -------------------------------------
  FS fs = get_FS(pid_h1,pid_h2);  
    
  // Create CLAS12Analysis Objects
  // -------------------------------------
  Kinematics _kin;
    
  // Create particle struct vector
  // -------------------------------------
  std::vector<part> mcparts;
    
  // Declare TLorentzVectors
  // -------------------------------------
  TLorentzVector vec_eIn, vec_eOut,q,h1,h2,dihadron;
  TLorentzVector init_target(0,0,0,Mp);  
    
  // Declare indecies for tracking dihadron particles
  // ------------------------------------------------
  std::vector<int> h1_idxs;
  std::vector<int> h2_idxs;
  std::vector<std::vector<int>> dihadron_idxs;    

  // For loop over events in hipo file
  // -------------------------------------    
  int whileidx=0;
    
  while(_chain.Next()==true && (whileidx < maxEvents || maxEvents < 0)){
    if(whileidx%10000==0 && whileidx!=0){
      std::cout << whileidx << " events read" << std::endl;
    }
    whileidx++;
      
    auto event = _c12->event();
    
    // Clear vectors
    mcparts.clear();
    h1_idxs.clear();
    h2_idxs.clear();
    dihadron_idxs.clear();
      
    // Get run specific information
    // -------------------------------------
    run = _c12->getBank(_idx_RUNconfig)->getInt(_irun,0);
      
    if(run!=11){
        std::cout << "File " << input_file << " is not Monte Carlo...Aborting..." << std::endl;
        return -1;
    }
      
    auto mcparticles=_c12->mcparts();
    for(int idx = 0 ; idx < mcparticles->getRows() ; idx++){
      int pid = mcparticles->getPid(idx);
      
      
      if(pid!=11 && pid!=fs.pid_h1 && pid!=fs.pid_h2) // Skip over unimportant particles for building dihadron
          continue;
      
      part partstruct;
      partstruct.pid = mcparticles->getPid(idx);
      partstruct.px = mcparticles->getPx(idx);
      partstruct.py = mcparticles->getPy(idx);
      partstruct.pz = mcparticles->getPz(idx);
      partstruct.m = mcparticles->getMass(idx);
      partstruct.pt = _kin.Pt(partstruct.px,partstruct.py);
      partstruct.p  = _kin.P(partstruct.px,partstruct.py,partstruct.pz);
      partstruct.E  = _kin.E(partstruct.m,partstruct.p);
      if(idx==0) {
          vec_eIn.SetPxPyPzE(0,0,partstruct.pz,partstruct.E);
          _electron_beam_energy=partstruct.E;
      }
      
      partstruct.theta = _kin.th(partstruct.pt,partstruct.pz);
      partstruct.eta = _kin.eta(partstruct.theta);
      partstruct.phi   = _kin.phi(partstruct.px,partstruct.py);

      partstruct.trueparentid = mcparticles->getParent(idx)-1;
      partstruct.trueparentpid = mcparticles->getPid(partstruct.trueparentid);
        
      if(partstruct.pid==11 && partstruct.trueparentid==0){ // scattered electron
        s=pow(Mp,2)+pow(Me,2)+2*Mp*_electron_beam_energy;
        Q2=_kin.Q2(_electron_beam_energy,partstruct.E,_kin.cth(partstruct.px,partstruct.py,partstruct.pz));
        y=_kin.y(_electron_beam_energy,partstruct.E);
        nu=_kin.nu(_electron_beam_energy,partstruct.E);
        W=_kin.W(Q2,Mp,nu);
        x=_kin.x(Q2,s,y);
        partstruct.is_scattered_electron=1;
        vec_eOut.SetPxPyPzE(partstruct.px,partstruct.py,partstruct.pz,partstruct.E);
      }
        
      // Only add dihadron worthy particles to the list
      if(pid==11) continue;
      if(pid==22 && partstruct.trueparentpid!=111) continue; // skip photons that don't come from pi0's
      
      mcparts.push_back(partstruct);
    }
    

      
    q=vec_eIn-vec_eOut;
    
    int Nmax = mcparts.size();
    //Loop over all particles in the event to determine hadron indecies

        for(int i = 0; i<Nmax; i++){
            if(mcparts[i].pid==pid_h1 || (mcparts[i].pid==22 && pid_h1==111)) h1_idxs.push_back(i);
            if(mcparts[i].pid==pid_h2 || (mcparts[i].pid==22 && pid_h2==111)) h2_idxs.push_back(i);
        }
        //Now form all possible dihadron index pairs
        if(pid_h1==pid_h2 && pid_h1!=111){dihadron_idxs=unique_combinations(h1_idxs,2); }
        else if(pid_h1==pid_h2 && pid_h1==111){dihadron_idxs=unique_combinations(h1_idxs,4);}
        else if(pid_h1!=pid_h2 && pid_h1==111 && pid_h2 != 111){
            for(int i = 0 ; i < h2_idxs.size(); i++){
                for(int j = 0 ; j < h1_idxs.size(); j++){
                    for(int k = j+1 ; k < h1_idxs.size(); k++){
                        std::vector<int> dihadron_idx = {h1_idxs.at(j),h1_idxs.at(k),h2_idxs.at(i)}; // 2 photons at start
                        dihadron_idxs.push_back(dihadron_idx);
                    }
                }
            }
        }
        else if(pid_h1!=pid_h2 && pid_h1!=111 && pid_h2 == 111){
            for(int i = 0 ; i < h1_idxs.size(); i++){
                for(int j = 0 ; j < h2_idxs.size(); j++){
                    for(int k = j+1 ; k < h2_idxs.size(); k++){
                        std::vector<int> dihadron_idx = {h1_idxs.at(i),h2_idxs.at(j),h2_idxs.at(k)}; // 2 photons at end
                        dihadron_idxs.push_back(dihadron_idx);
                    }
                }
            }
        }
        else{
            for(int i = 0 ; i < h1_idxs.size(); i++){
                for(int j = 0 ; j < h2_idxs.size(); j++){
                    std::vector<int> dihadron_idx = {h1_idxs.at(i), h2_idxs.at(j)};
                    dihadron_idxs.push_back(dihadron_idx);
                }
            }
        }
    
        // Remove any instance of duplicate dihadrons
        dihadron_idxs = remove_duplicates(dihadron_idxs);
        // Now loop over all dihadrons
        for(int a = 0 ; a < dihadron_idxs.size() ; a++){
            std::vector<int> dihadron_idx = dihadron_idxs.at(a);
            int i=0;
            int ii=0;
            int j=0;
            int jj=0;
            if(pid_h1==111){
                i=dihadron_idx.at(0);
                ii=dihadron_idx.at(1);
            }else{
                i=dihadron_idx.at(0);
            }
            if(pid_h2==111&&pid_h1!=111){
                j=dihadron_idx.at(1);
                jj=dihadron_idx.at(2);
            }else if(pid_h2==111&&pid_h1==111){
                j=dihadron_idx.at(2);
                jj=dihadron_idx.at(3);
            }else if(pid_h1==111){
                j=dihadron_idx.at(2);
            }else{
                j=dihadron_idx.at(1);
            }
  
            if(pid_h1==111){
                h1.SetPxPyPzE(mcparts[i].px+mcparts[ii].px,mcparts[i].py+mcparts[ii].py,mcparts[i].pz+mcparts[ii].pz,mcparts[i].E+mcparts[ii].E);
            }else{
                h1.SetPxPyPzE(mcparts[i].px,mcparts[i].py,mcparts[i].pz,mcparts[i].E);
            }
            if(pid_h2==111){
                h2.SetPxPyPzE(mcparts[j].px+mcparts[jj].px,mcparts[j].py+mcparts[jj].py,mcparts[j].pz+mcparts[jj].pz,mcparts[j].E+mcparts[jj].E);
            }
            else{
                h2.SetPxPyPzE(mcparts[j].px,mcparts[j].py,mcparts[j].pz,mcparts[j].E);
            }
            
            if(pid_h1==pid_h2){
                z1 = _kin.z(init_target,h1,q);
                z2 = _kin.z(init_target,h2,q);
                if(z1<z2){
                    TLorentzVector temp = h1;
                    h2=h1;
                    h1=temp;
                }
            }
            
            // Build the dihadron
            dihadron = h1+h2;
            // fill results
            Mh = dihadron.M();
            pT1 = _kin.Pt(q,h1,init_target);
            pT2 = _kin.Pt(q,h2,init_target);
            pTtot = _kin.Pt(q,dihadron,init_target);
            phi_h = _kin.phi_h(q,vec_eIn,h1,h2);
            phi_R0 = _kin.phi_R(q,vec_eIn,h1,h2,0);
            phi_R1 = _kin.phi_R(q,vec_eIn,h1,h2,1);
            th     = _kin.com_th(h1,h2);
            xF1 = _kin.xF(q,h1,init_target,W);
            xF2 = _kin.xF(q,h2,init_target,W);
            xF     = _kin.xF(q,dihadron,init_target,W);
            z1 = _kin.z(init_target,h1,q);
            z2 = _kin.z(init_target,h2,q);
            z = z1+z2;
            Mx = (vec_eIn+init_target-vec_eOut-dihadron).M();

            P1=h1.P();
            P2=h2.P();
            if(y>0.8) continue;
            if(W<2) continue;
            if(Q2<1) continue;
            if(pid_h1!=111&&P1<1.25) continue;
            if(pid_h2!=111&&P2<1.25) continue;
            if((pid_h1==211&&pid_h2==-211)||(pid_h1==211&&pid_h2==111)){
                if(Mx<1.5) continue;
            }
            if(z>0.95) continue;
            if(xF1<0||xF2<0) continue;
            tree->Fill();
        } // end dihadron loop
    }// end event loop
    tree->Write();
    fOut->Close();
    return 0;
}