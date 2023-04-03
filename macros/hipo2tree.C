#include "../src/CutManager.C"
#include "../src/HipoBankInterface.C"
#include "../src/Constants.h"
#include "../src/Structs.h"
#include "../src/Kinematics.C"
#include "../src/ParseBinYAML.C"
#include "../src/ParseText.C"


int hipo2tree(//const char * hipoFile = "/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v1/dst/train/nSidis/nSidis_005032.hipo",
              const char * hipoFile = "/cache/clas12/rg-a/production/montecarlo/clasdis/fall2018/torus-1/v1/bkg45nA_10604MeV/45nA_job_3051_0.hipo",
              const char * outputFile = "hipo2tree.root",
              const double _electron_beam_energy = 10.6,
              const int pid_h1=211,
              const int pid_h2=111,
              const int maxEvents = 10000,
              bool hipo_is_mc = true){
 
  // Open TTree and declare branches
  // -------------------------------------
  TFile *fOut = new TFile(outputFile,"RECREATE");
  TTree* tree = new TTree("EventTree","EventTree");
    
  int run;
  // Reconstructed variables
  double x, y, W, nu,Pol,Q2;
  int hel;
  double s=pow(Mp,2)+pow(Me,2)+2*Mp*_electron_beam_energy;

  // Monte Carlo Info
  double truex, truey, trueW, truenu,trueQ2;

  // Particle Info
  int Nmax = 100; // maximum number of particles
  int pindex[Nmax], status[Nmax], pid[Nmax], truepid[Nmax], trueparentid[Nmax], trueparentpid[Nmax], trueparentparentid[Nmax], trueparentparentpid[Nmax];
  double px[Nmax], py[Nmax], pz[Nmax], p[Nmax], E[Nmax],m[Nmax];
  double vx[Nmax], vy[Nmax], vz[Nmax], chi2[Nmax], beta[Nmax];
  double truepx[Nmax], truepy[Nmax], truepz[Nmax], truep[Nmax], trueE[Nmax];
  int is_CFR[Nmax];
  double theta[Nmax], eta[Nmax], phi[Nmax], truept[Nmax], truem[Nmax], truetheta[Nmax], trueeta[Nmax], truephi[Nmax], truevx[Nmax], truevy[Nmax], truevz[Nmax];
  int pcal_sector[Nmax], ecin_sector[Nmax], ecout_sector[Nmax];
  double pcal_x[Nmax], pcal_y[Nmax], pcal_z[Nmax];
  double ecin_x[Nmax], ecin_y[Nmax], ecin_z[Nmax];
  double ecout_x[Nmax], ecout_y[Nmax], ecout_z[Nmax];
  double pcal_e[Nmax], pcal_lu[Nmax], pcal_lv[Nmax], pcal_lw[Nmax], pcal_m2u[Nmax], pcal_m2v[Nmax], pcal_m2w[Nmax];
  double ecin_e[Nmax], ecin_lu[Nmax], ecin_lv[Nmax], ecin_lw[Nmax], ecin_m2u[Nmax], ecin_m2v[Nmax], ecin_m2w[Nmax];
  double ecout_e[Nmax], ecout_lu[Nmax], ecout_lv[Nmax], ecout_lw[Nmax], ecout_m2u[Nmax], ecout_m2v[Nmax], ecout_m2w[Nmax];
  int sector[Nmax];
  double traj_x1[Nmax], traj_y1[Nmax], traj_z1[Nmax], traj_x2[Nmax], traj_y2[Nmax], traj_z2[Nmax], traj_x3[Nmax], traj_y3[Nmax], traj_z3[Nmax];
 
  // Set branches
  tree->Branch("run",&run,"run/I");
  tree->Branch("Pol",&Pol,"Pol/D");
  tree->Branch("Nmax",&Nmax,"Nmax/I");
  tree->Branch("x", &x, "x/D");
  tree->Branch("y", &y, "y/D");
  tree->Branch("s", &s, "s/D");
  tree->Branch("W", &W, "W/D");
  tree->Branch("Q2",&Q2,"Q2/D");
  tree->Branch("nu", &nu, "nu/D");
  tree->Branch("truex", &truex, "truex/D");
  tree->Branch("truey", &truey, "truey/D");
  tree->Branch("trueQ2",&trueQ2,"trueQ2/D");
  tree->Branch("trueW", &trueW, "trueW/D");
  tree->Branch("truenu", &truenu, "truenu/D");
  tree->Branch("hel", &hel, "hel/I");

  tree->Branch("truex", &truex, "truex/D");
  tree->Branch("truey", &truey, "truey/D");
  tree->Branch("trueW", &trueW, "trueW/D");
  tree->Branch("truenu", &truenu, "truenu/D");

  tree->Branch("pindex", pindex, "pindex[Nmax]/I");
  tree->Branch("status", status, "status[Nmax]/I");
  tree->Branch("px", px, "px[Nmax]/D");
  tree->Branch("py", py, "py[Nmax]/D");
  tree->Branch("pz", pz, "pz[Nmax]/D");
  tree->Branch("p", p, "p[Nmax]/D");
  tree->Branch("E", E, "E[Nmax]/D");
  tree->Branch("pid", pid, "pid[Nmax]/I");
  tree->Branch("vx", vx, "vx[Nmax]/D");
  tree->Branch("vy", vy, "vy[Nmax]/D");
  tree->Branch("vz", vz, "vz[Nmax]/D");
  tree->Branch("chi2", chi2, "chi2[Nmax]/D");
  tree->Branch("beta", beta, "beta[Nmax]/D");
  tree->Branch("m", m, "m[Nmax]/D");
  tree->Branch("theta", theta, "theta[Nmax]/D");
  tree->Branch("eta", eta, "eta[Nmax]/D");
  tree->Branch("phi", phi, "phi[Nmax]/D");
  tree->Branch("truepx", truepx, "truepx[Nmax]/D");
  tree->Branch("truepy", truepy, "truepy[Nmax]/D");
  tree->Branch("truepz", truepz, "truepz[Nmax]/D");
  tree->Branch("truep", truep, "truep[Nmax]/D");
  tree->Branch("truept", truept, "truept[Nmax]/D");
  tree->Branch("truem", truem, "truem[Nmax]/D");
  tree->Branch("truetheta", truetheta, "truetheta[Nmax]/D");
  tree->Branch("trueeta", trueeta, "trueeta[Nmax]/D");
  tree->Branch("truephi", truephi, "truephi[Nmax]/D");
  tree->Branch("truevx", truevx, "truevx[Nmax]/D");
  tree->Branch("truevy", truevy, "truevy[Nmax]/D");
  tree->Branch("truevz", truevz, "truevz[Nmax]/D");
  tree->Branch("trueE", trueE, "trueE[Nmax]/D");
  tree->Branch("is_CFR", is_CFR, "is_CFR[Nmax]/I");
  tree->Branch("truepid", truepid, "truepid[Nmax]/I");
  tree->Branch("trueparentid", trueparentid, "trueparentid[Nmax]/I");
  tree->Branch("trueparentpid", trueparentpid, "trueparentpid[Nmax]/I");
  tree->Branch("trueparentparentid", trueparentparentid, "trueparentparentid[Nmax]/I");
  tree->Branch("trueparentparentpid", trueparentparentpid, "trueparentparentpid[Nmax]/I");
    
  tree->Branch("pcal_sector", pcal_sector, "pcal_sector[Nmax]/I");
  tree->Branch("pcal_e", pcal_e, "pcal_e[Nmax]/D");
  tree->Branch("pcal_x", pcal_x, "pcal_x[Nmax]/D");
  tree->Branch("pcal_y", pcal_y, "pcal_y[Nmax]/D");
  tree->Branch("pcal_z", pcal_z, "pcal_z[Nmax]/D");
  tree->Branch("pcal_lu", pcal_lu, "pcal_lu[Nmax]/D");
  tree->Branch("pcal_lv", pcal_lv, "pcal_lv[Nmax]/D");
  tree->Branch("pcal_lw", pcal_lw, "pcal_lw[Nmax]/D");
  tree->Branch("pcal_m2u", pcal_m2u, "pcal_m2u[Nmax]/D");
  tree->Branch("pcal_m2v", pcal_m2v, "pcal_m2v[Nmax]/D");
  tree->Branch("pcal_m2w", pcal_m2w, "pcal_m2w[Nmax]/D");

  tree->Branch("ecin_sector", ecin_sector, "ecin_sector[Nmax]/I");
  tree->Branch("ecin_e", ecin_e, "ecin_e[Nmax]/D");
  tree->Branch("ecin_x", ecin_x, "ecin_x[Nmax]/D");
  tree->Branch("ecin_y", ecin_y, "ecin_y[Nmax]/D");
  tree->Branch("ecin_z", ecin_z, "ecin_z[Nmax]/D");
  tree->Branch("ecin_lu", ecin_lu, "ecin_lu[Nmax]/D");
  tree->Branch("ecin_lv", ecin_lv, "ecin_lv[Nmax]/D");
  tree->Branch("ecin_lw", ecin_lw, "ecin_lw[Nmax]/D");
  tree->Branch("ecin_m2u", ecin_m2u, "ecin_m2u[Nmax]/D");
  tree->Branch("ecin_m2v", ecin_m2v, "ecin_m2v[Nmax]/D");
  tree->Branch("ecin_m2w", ecin_m2w, "ecin_m2w[Nmax]/D");

  tree->Branch("ecout_sector", ecout_sector, "ecout_sector[Nmax]/I");
  tree->Branch("ecout_e", ecout_e, "ecout_e[Nmax]/D");
  tree->Branch("ecout_x", ecout_x, "ecout_x[Nmax]/D");
  tree->Branch("ecout_y", ecout_y, "ecout_y[Nmax]/D");
  tree->Branch("ecout_z", ecout_z, "ecout_z[Nmax]/D");
  tree->Branch("ecout_lu", ecout_lu, "ecout_lu[Nmax]/D");
  tree->Branch("ecout_lv", ecout_lv, "ecout_lv[Nmax]/D");
  tree->Branch("ecout_lw", ecout_lw, "ecout_lw[Nmax]/D");
  tree->Branch("ecout_m2u", ecout_m2u, "ecout_m2u[Nmax]/D");
  tree->Branch("ecout_m2v", ecout_m2v, "ecout_m2v[Nmax]/D");
  tree->Branch("ecout_m2w", ecout_m2w, "ecout_m2w[Nmax]/D");

  tree->Branch("sector", sector, "sector[Nmax]/I");
  tree->Branch("traj_x1", traj_x1, "traj_x1[Nmax]/D");
  tree->Branch("traj_y1", traj_y1, "traj_y1[Nmax]/D");
  tree->Branch("traj_z1", traj_z1, "traj_z1[Nmax]/D");
  tree->Branch("traj_x2", traj_x2, "traj_x2[Nmax]/D");
  tree->Branch("traj_y2", traj_y2, "traj_y2[Nmax]/D");
  tree->Branch("traj_z2", traj_z2, "traj_z2[Nmax]/D");
  tree->Branch("traj_x3", traj_x3, "traj_x3[Nmax]/D");
  tree->Branch("traj_y3", traj_y3, "traj_y3[Nmax]/D");
  tree->Branch("traj_z3", traj_z3, "traj_z3[Nmax]/D");
    
    
  // Configure CLAS12 Reader and HipoChain
  // -------------------------------------
  clas12root::HipoChain _chain;
  clas12::clas12reader *_config_c12{nullptr};

  _chain.Add(hipoFile);
  _config_c12=_chain.GetC12Reader();
  // If not monte carlo, enforce QADB
  // -------------------------------------
  if(hipo_is_mc==true)
    _config_c12->db()->turnOffQADB();
  

  // Configure PIDs for final state
  // -------------------------------------
  FS fs = get_FS(pid_h1,pid_h2);
  _config_c12->addExactPid(11,1);     // Exactly 1 electron
  _config_c12->addAtLeastPid(fs.pid_h1,fs.num_h1);
  if(fs.pid_h2!=0) _config_c12->addAtLeastPid(fs.pid_h2,fs.num_h2); // Doesn't run if duplicate final state
  
  // Add RUN::config bank
  // -------------------------------------
  int _idx_RUNconfig = _config_c12->addBank("RUN::config");
  int _irun = _config_c12->getBankOrder(_idx_RUNconfig,"run");
  int _ievnum = _config_c12->getBankOrder(_idx_RUNconfig,"event");
  int _itorus = _config_c12->getBankOrder(_idx_RUNconfig,"torus");

  // Establish CLAS12 event parser
  // -------------------------------------
  auto &_c12=_chain.C12ref();
  
  // Create CLAS12Analysis Objects
  // -------------------------------------
  HipoBankInterface _hipoInterface = HipoBankInterface(_c12);
  CutManager _cm = CutManager();
  Kinematics _kin;

  // Create particle struct vector
  // -------------------------------------
  std::vector<part> vec_particles;
  std::vector<part> vec_mcparticles;

  int whileidx=0;
  int _ievent=0;
  while(_chain.Next()==true && (whileidx < maxEvents || maxEvents < 0)){
    if(whileidx%10000==0 && whileidx!=0){
      std::cout << whileidx << " events read | " << _ievent*100.0/whileidx << "% passed event selection" << std::endl;
    }

    auto event = _c12->event();
      
    whileidx++;
    
    // Clear vectors
    vec_particles.clear();
    vec_mcparticles.clear();
      
    // Get run specific information
    // -------------------------------------
    run = _c12->getBank(_idx_RUNconfig)->getInt(_irun,0);
    if(hipo_is_mc)
      run *= _c12->getBank(_idx_RUNconfig)->getFloat(_itorus,0); // Multiply run number by torus bending
    _cm.set_run(run);
      
    // Get helicity
    // -------------------------------------
    
    if(runHelicityFlip(run)){
      hel = -event->getHelicity();
    }
    else
      {
        hel = event->getHelicity();
      }
      
    
    // Get polarization
    // -------------------------------------
    Pol = runPolarization(run);
      
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
      if(partstruct.theta*180/PI<5 || partstruct.theta*180/PI>35) continue;
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
      
      if(partstruct.pid==11){//only electron is assumed to be the scattered e
	Q2=_kin.Q2(_electron_beam_energy,partstruct.E,_kin.cth(partstruct.px,partstruct.py,partstruct.pz));
	y=_kin.y(_electron_beam_energy,partstruct.E);
        if(y>0.8)
	  continue;
	nu=_kin.nu(_electron_beam_energy,partstruct.E);
	W=_kin.W(Q2,Mp,nu);
	x=_kin.x(Q2,s,y);
      }
      _hipoInterface.loadBankData(_c12,partstruct);
      vec_particles.push_back(partstruct);
          
    }

    // Apply the cuts to make a new vec_particles
    // Calls the CutManager which parses through the vector
    // and makes relevant cuts for each particle
    // 
    // Return a list of filtered particles
    // -------------------------------------------------
      
    vec_particles = _cm.filter_particles(vec_particles);
      
    // Determine if this filtered list still contains the final state
    // particles we are interested in. If not, next event
    // -------------------------------------------------
    int num_e=0;
    int num_h1=0;
    int num_h2=0;
    for(part particle : vec_particles){
      if(particle.pid==11) num_e++;
      else if(particle.pid==fs.pid_h1) num_h1++;
      else if(particle.pid==fs.pid_h2) num_h2++;
    }
    if(num_e<1 || num_h1 < fs.num_h1 || num_h2 < fs.num_h2) {
      continue;
    }
      
      
      
    // Loop over all Monte Carlo particles
    // -------------------------------------
    auto mcparticles=_c12->mcparts();
    for(int idx = 0 ; idx < mcparticles->getRows() && hipo_is_mc ; idx++){
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
      if(partstruct.truepid==11 && partstruct.trueparentid==0){ // scattered electron
        trueQ2=_kin.Q2(_electron_beam_energy,partstruct.trueE,_kin.cth(partstruct.truepx,partstruct.truepy,partstruct.truepz));
        truey=_kin.y(_electron_beam_energy,partstruct.trueE);
        truenu=_kin.nu(_electron_beam_energy,partstruct.trueE);
        trueW=_kin.W(trueQ2,Mp,truenu);
        truex=_kin.x(trueQ2,s,truey);
      }
      
      // Add particle to list
      vec_mcparticles.push_back(partstruct);
    }
    // Perform MC<-->Reco particle matching
    // ---------------------------------------------
    for (int i=0; i < vec_particles.size(); i++){
      for (int j=0; j < vec_mcparticles.size(); j++){
        float dth = abs(vec_particles[i].theta - vec_mcparticles[j].truetheta)*180/PI;
        float dphi = abs(vec_particles[i].phi - vec_mcparticles[j].truephi)*180/PI;
        float dE = abs(vec_particles[i].E - vec_mcparticles[j].trueE);
        
        if (dth<2 && (dphi<4 || abs(dphi-2*PI)<4) && dE<0.5){
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
      
    // Loop over all particles and fill variables for the TTree
    // --------------------------------------------------------
    Nmax=vec_particles.size();
    for(int i = 0; i < vec_particles.size(); i++) {
      part par = vec_particles[i];
        
      pindex[i] = par.pindex;
      status[i] = par.status;
      px[i] = par.px;
      py[i] = par.py;
      pz[i] = par.pz;
      p[i] = par.p;
      E[i] = par.E;
      pid[i] = par.pid;
      vx[i] = par.vx;
      vy[i] = par.vy;
      vz[i] = par.vz;
      chi2[i] = par.chi2;
      beta[i] = par.beta;
      m[i] = par.m;
      theta[i] = par.theta;
      eta[i] = par.eta;
      phi[i] = par.phi;
      truepx[i] = par.truepx;
      truepy[i] = par.truepy;
      truepz[i] = par.truepz;
      truep[i] = par.truep;
      truept[i] = par.truept;
      trueE[i] = par.trueE;
      truem[i] = par.truem;
      truetheta[i] = par.truetheta;
      trueeta[i] = par.trueeta;
      truephi[i] = par.truephi;
      truevx[i] = par.truevx;
      truevy[i] = par.truevy;
      truevz[i] = par.truevz;
      is_CFR[i] = par.is_CFR;
      truepid[i] = par.truepid;
      trueparentid[i] = par.trueparentid;
      trueparentpid[i] = par.trueparentpid;
      trueparentparentid[i] = par.trueparentparentid;
      trueparentparentpid[i] = par.trueparentparentpid;
      pcal_sector[i] = par.pcal_sector;
      pcal_e[i] = par.pcal_e;
      pcal_x[i] = par.pcal_x;
      pcal_y[i] = par.pcal_y;
      pcal_z[i] = par.pcal_z;
      pcal_lu[i] = par.pcal_lu;
      pcal_lv[i] = par.pcal_lv;
      pcal_lw[i] = par.pcal_lw;
      pcal_m2u[i] = par.pcal_m2u;
      pcal_m2v[i] = par.pcal_m2v;
      pcal_m2w[i] = par.pcal_m2w;
      ecin_sector[i] = par.ecin_sector;
      ecin_e[i] = par.ecin_e;
      ecin_x[i] = par.ecin_x;
      ecin_y[i] = par.ecin_y;
      ecin_z[i] = par.ecin_z;
      ecin_lu[i] = par.ecin_lu;
      ecin_lv[i] = par.ecin_lv;
      ecin_lw[i] = par.ecin_lw;
      ecin_m2u[i] = par.ecin_m2u;
      ecin_m2v[i] = par.ecin_m2v;
      ecin_m2w[i] = par.ecin_m2w;
      ecout_sector[i] = par.ecout_sector;
      ecout_e[i] = par.ecout_e;
      ecout_x[i] = par.ecout_x;
      ecout_y[i] = par.ecout_y;
      ecout_z[i] = par.ecout_z;
      ecout_lu[i] = par.ecout_lu;
      ecout_lv[i] = par.ecout_lv;
      ecout_lw[i] = par.ecout_lw;
      ecout_m2u[i] = par.ecout_m2u;
      ecout_m2v[i] = par.ecout_m2v;
      ecout_m2w[i] = par.ecout_m2w;
      sector[i] = par.sector;
      traj_x1[i] = par.traj_x1;
      traj_y1[i] = par.traj_y1;
      traj_z1[i] = par.traj_z1;
      traj_x2[i] = par.traj_x2;
      traj_y2[i] = par.traj_y2;
      traj_z2[i] = par.traj_z2;
      traj_x3[i] = par.traj_x3;
      traj_y3[i] = par.traj_y3;
      traj_z3[i] = par.traj_z3;
    }
    tree->Fill();
    _ievent++;
  }
  tree->Write();
  fOut->Close(); 
  return 0;
}
