#include "CLAS12Analysis.h"

CLAS12Analysis::CLAS12Analysis(){}

CLAS12Analysis::CLAS12Analysis(const std::unique_ptr<clas12::clas12reader>& _c12, TLorentzVector eIn, TLorentzVector pIn){
    hipoBankInterface = HipoBankInterface(_c12);
    init_electron = eIn;
    target = pIn;
    _electron_beam_energy = eIn.E();
    s = init_electron.M2()+target.M2()+2*target.M()*_electron_beam_energy;
}

CLAS12Analysis::CLAS12Analysis(const std::unique_ptr<clas12::clas12reader>& _c12, double beamE){
    hipoBankInterface = HipoBankInterface(_c12);
    init_electron.SetPxPyPzE(0,0,sqrt(beamE*beamE-Me*Me),beamE);
    target.SetPxPyPzE(0,0,0,Mp);
    _electron_beam_energy = beamE;
    s = init_electron.M2()+target.M2()+2*target.M()*_electron_beam_energy;
}


std::vector<part> CLAS12Analysis::load_reco_particles(const std::unique_ptr<clas12::clas12reader>& _c12){
    
    std::vector<part> vec_particles;
    
    // Loop over reconstructed particles
    // -------------------------------------------------------
    auto particles=_c12->getDetParticles();
    for(unsigned int idx = 0 ; idx < particles.size() ; idx++){
      // Create new part struct
      part partstruct;
      // Extract each particle from event one-at-a-time
      // -------------------------------------------------------
      auto particle = particles.at(idx);
      partstruct.pid = particle->getPid(); 

      partstruct.chi2 = particle->getChi2Pid();
      partstruct.theta = particle->getTheta();
      partstruct.eta = _kin.eta(partstruct.theta);
      partstruct.phi = particle->getPhi();
      partstruct.p = particle->getP();
      partstruct.px = _kin.Px(partstruct.p,partstruct.theta,partstruct.phi);
      partstruct.py = _kin.Py(partstruct.p,partstruct.theta,partstruct.phi);
      partstruct.pz = _kin.Pz(partstruct.p,partstruct.theta,partstruct.phi);
      partstruct.pt = _kin.Pt(partstruct.px,partstruct.py);
      if(partstruct.pid!=22)
        partstruct.m = particle->getPdgMass();
      else
        partstruct.m = 0;
      partstruct.beta = particle->getBeta();
      partstruct.pindex = particle->getIndex();
      partstruct.vx = particle->par()->getVx();
      partstruct.vy = particle->par()->getVy();
      partstruct.vz = particle->par()->getVz();
      partstruct.status = particle->getStatus();
      partstruct.E = _kin.E(partstruct.m,partstruct.p);
    
      // Ensure hadrons are not in CD
      if (partstruct.pid == 2212 || partstruct.pid == -2212 ||
        partstruct.pid == 2112 ||
        partstruct.pid == -321 || partstruct.pid == -211 ||
        partstruct.pid == 211 || partstruct.pid == 321) {
          if(partstruct.status>=4000 && partstruct.status<5000)
              continue;
      }
    
      hipoBankInterface.loadBankData(_c12,partstruct);  
      vec_particles.push_back(partstruct);

    }
    
    return vec_particles;
}



std::vector<part> CLAS12Analysis::load_mc_particles(const std::unique_ptr<clas12::clas12reader>& _c12){
    
    
    std::vector<part> vec_mcparticles;
    
      
    // Loop over all Monte Carlo particles
    // -------------------------------------
    auto mcparticles=_c12->mcparts();
    for(int idx = 0 ; idx < mcparticles->getRows(); idx++){
      part partstruct;
      if(mcparticles->getType(idx)!=1) // Reject non-final state
        {continue;} 
      partstruct.truepid = mcparticles->getPid(idx);
      partstruct.truepx = mcparticles->getPx(idx);
      partstruct.truepy = mcparticles->getPy(idx);
      partstruct.truepz = mcparticles->getPz(idx);
      partstruct.truem = mcparticles->getMass(idx);
      partstruct.truept = _kin.Pt(partstruct.truepx,partstruct.truepy);
      partstruct.truep  = _kin.P(partstruct.truepx,partstruct.truepy,partstruct.truepz);
      partstruct.trueE  = _kin.E(partstruct.truem,partstruct.truep);

      partstruct.truetheta = _kin.th(partstruct.truept,partstruct.truepz);
      partstruct.trueeta = _kin.eta(partstruct.truetheta);
      partstruct.truephi   = _kin.phi(partstruct.truepx,partstruct.truepy);

      partstruct.truevx = mcparticles->getVx(idx);
      partstruct.truevy = mcparticles->getVy(idx);
      partstruct.truevz = mcparticles->getVz(idx);

      partstruct.trueparentid = mcparticles->getParent(idx)-1;
      partstruct.trueparentpid = mcparticles->getPid(partstruct.trueparentid);
      partstruct.trueparentparentid = mcparticles->getParent(partstruct.trueparentid)-1;
      if(partstruct.trueparentparentid==-1){
          partstruct.trueparentparentpid = -999;
      }else{
          partstruct.trueparentparentpid = mcparticles->getPid(partstruct.trueparentparentid);
      }
      // for loop over the idxs until we find if this particle came from CFR
      int parent_idx = mcparticles->getParent(idx)-1;
      int parent_pid = 0;
      while(parent_idx>=0){
          parent_pid = mcparticles->getPid(parent_idx);
          if(parent_pid==0){break;}
          if(6-abs(parent_pid)>=0){
              partstruct.is_CFR=1;
              break;
          }
          parent_idx = mcparticles->getParent(parent_idx)-1;
      }
      if(partstruct.is_CFR!=1) partstruct.is_CFR=0;
      if(partstruct.truepid==11 && partstruct.trueparentid==0){ // scattered electro
        partstruct.is_scattered_electron=1;
      }
      
      // Add particle to list
      vec_mcparticles.push_back(partstruct);
    }
    
    return vec_mcparticles;
    
}

bool CLAS12Analysis::reco_event_contains_final_state(std::vector<part> vec_particles, FS fs){
    int num_e=0;
    int num_h1=0;
    int num_h2=0;
    for(part particle : vec_particles){
      if(particle.pid==11 && particle.is_scattered_electron==1) num_e++;
      else if(particle.pid==fs.pid_h1) num_h1++;
      else if(particle.pid==fs.pid_h2) num_h2++;
    }
    if(num_e<1 || num_h1 < fs.num_h1 || num_h2 < fs.num_h2) {
      return false;
    }
    
    return true;

}

int CLAS12Analysis::find_reco_scattered_electron(std::vector<part>& vec_particles){
    // Code for determine the scattered electron from REC::Particle
    // --> Find pid==11 particle with largest energy
    // -->   If no electron is found, skip
    // --> Check if the status of the maximum energy electron is in FD
    // -->   Skip if not (i.e. always skip events if the max energy electron was not in FD)
    // --> Set that particle as the scattered electron
    int idx_e=-1;
    double max_energy = -1; 
    for (int i = 0; i < vec_particles.size(); i++) {
      part partstruct = vec_particles[i];
      // check if the particle is an electron
      if (partstruct.pid == 11) {
        // compare energy with the current maximum and update if necessary
        if (partstruct.E > max_energy) {
          max_energy = partstruct.E;
          idx_e=i;
         }
       }
    }
    
    
    
    if(idx_e==-1) return -1;
    
    
    // Toss events where the REC::Particle scattered electron was not in the FD
    
    if((vec_particles[idx_e].status <= -3000 || vec_particles[idx_e].status > -2000))
            return -1;
    
    return idx_e;
}

DIS_EVENT CLAS12Analysis::calc_reco_event_variables(std::vector<part> parts){
    
    part scattered_electron;
    for(auto p: parts){
        if(p.is_scattered_electron){
            scattered_electron = p;
            break;
        }
    }
    
    // If none was found, print an error message
    if(scattered_electron.is_scattered_electron != 1){
        cout << "ERROR: No rec::scattered electron found..." << endl;
    }
    
    DIS_EVENT event;
    event.Q2=_kin.Q2(_electron_beam_energy,scattered_electron.E,
               _kin.cth(scattered_electron.px,scattered_electron.py,scattered_electron.pz));
    event.y=_kin.y(_electron_beam_energy,scattered_electron.E);     
    event.nu=_kin.nu(_electron_beam_energy,scattered_electron.E);
    event.W=_kin.W(event.Q2,Mp,event.nu);
    event.x=_kin.x(event.Q2,s,event.y);

    return event;
}

DIS_EVENT CLAS12Analysis::calc_mc_event_variables(std::vector<part> parts){
    
    part scattered_electron;
    for(auto p: parts){
        if(p.is_scattered_electron){
            scattered_electron = p;
            break;
        }
    }
    
    // If none was found, print an error message
    if(scattered_electron.is_scattered_electron != 1){
        cout << "ERROR: No mc::scattered electron found..." << endl;
    }
    
    DIS_EVENT event;
    event.trueQ2=_kin.Q2(_electron_beam_energy,scattered_electron.trueE,
               _kin.cth(scattered_electron.truepx,scattered_electron.truepy,scattered_electron.truepz));
    event.truey=_kin.y(_electron_beam_energy,scattered_electron.trueE);     
    event.truenu=_kin.nu(_electron_beam_energy,scattered_electron.trueE);
    event.trueW=_kin.W(event.trueQ2,Mp,event.truenu);
    event.truex=_kin.x(event.trueQ2,s,event.truey);
    
    return event;
}


void CLAS12Analysis::match_mc_to_reco(std::vector<part>& vec_particles,
                                      std::vector<part>& vec_mcparticles){
    
    for (int i=0; i < vec_particles.size(); i++){
      for (int j=0; j < vec_mcparticles.size(); j++){
        float dth = abs(vec_particles[i].theta - vec_mcparticles[j].truetheta)*180/PI;
        float dphi = abs(vec_particles[i].phi - vec_mcparticles[j].truephi)*180/PI;
        float dE = abs(vec_particles[i].E - vec_mcparticles[j].trueE);
        
        if (dth<2 && (dphi<4 || abs(dphi-2*PI)<4) && dE<1){
	  // Perform Pairing
	  vec_particles[i].truepx = vec_mcparticles[j].truepx;
	  vec_particles[i].truepy = vec_mcparticles[j].truepy;
	  vec_particles[i].truepz = vec_mcparticles[j].truepz;
	  vec_particles[i].truep = vec_mcparticles[j].truep;
	  vec_particles[i].truept = vec_mcparticles[j].truept;
	  vec_particles[i].trueE = vec_mcparticles[j].trueE;
	  vec_particles[i].truem = vec_mcparticles[j].truem;
	  vec_particles[i].truetheta = vec_mcparticles[j].truetheta;
	  vec_particles[i].trueeta = vec_mcparticles[j].trueeta;
	  vec_particles[i].truephi = vec_mcparticles[j].truephi;
	  vec_particles[i].truevx = vec_mcparticles[j].truevx;
	  vec_particles[i].truevy = vec_mcparticles[j].truevy;
	  vec_particles[i].truevz = vec_mcparticles[j].truevz;
      vec_particles[i].is_CFR = vec_mcparticles[j].is_CFR;
	  vec_particles[i].truepid = vec_mcparticles[j].truepid;
	  vec_particles[i].trueparentid = vec_mcparticles[j].trueparentid;
	  vec_particles[i].trueparentpid = vec_mcparticles[j].trueparentpid;
      vec_particles[i].trueparentparentid = vec_mcparticles[j].trueparentparentid;
	  vec_particles[i].trueparentparentpid = vec_mcparticles[j].trueparentparentpid;
	  break;
	}
      }
    }
    
    
}

