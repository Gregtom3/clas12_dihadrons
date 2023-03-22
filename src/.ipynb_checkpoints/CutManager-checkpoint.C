#include "CutManager.h"
#include "Structs.h"

// Constructors
CutManager::CutManager(){
    _run = 0;
    _torusBending = 0;
}

CutManager::CutManager(int run){
    set_run(run);
}

// Public member functions
// Set the run number and torus for the CutManager
void CutManager::set_run(int run){
    _run=run;
    _torusBending=runTorusBending(run);
}

// Return a vector of particles that passes the cuts
std::vector<part> CutManager::filter_particles(std::vector<part> particles){
    
    std::vector<part> filtered_particles;
    // Get scattered electron
    part electron;
    for (auto particle : particles) {
        if(particle.pid==11) electron=particle;
    }
    
    for (auto particle : particles) {
        bool pass = false;
        int pid = particle.pid;
        
        switch(pid){
            case 11:
                pass = apply_electron_cuts(particle);
                break;
            case 211:
                pass = apply_pion_cuts(particle);
                break;
            case -211:
                pass = apply_pion_cuts(particle);
                break;
            case 22:
                pass = apply_photon_cuts(particle,electron);
                break;
            default:
                pass = true;
                break;
        }
        
        if (pass==true)
            filtered_particles.push_back(particle);
    }

    return filtered_particles;
}


// Protected member functions
// Apply all relevant electron cuts
bool CutManager::apply_electron_cuts(part particle){
    if(DC_fiducial_cut(particle)==false) return false;
    if(VzCut(particle)==false) return false;
    if(EleSampFrac(particle)==false) return false;
    if(minEpcal(particle)==false) return false;
    if(caloEdges(particle,0)==false) return false;
    if(Ele3calo(particle)==false) return false;
    return true;
}

// Apply all relevant pion cuts
bool CutManager::apply_pion_cuts(part particle){
    if(DC_fiducial_cut(particle)==false) return false;
    //if(chi2pid(particle,0)==false) return false;
    //if(minPiMomentum(particle)==false) return false;
    //if(caloEdges(particle,0)==false) return false;
    return true;
}

// Apply all relevant photon cuts
bool CutManager::apply_photon_cuts(part particle, part electron){
    if(minEpcal(particle)==false) return false;
    if(photonMinEtot(particle)==false) return false;
    if(photonElectronAngle(particle,electron)==false) return false;
    if(photonBetaCut(particle)==false) return false;
    if(caloEdges(particle,1)==false) return false;
    return true;
}

// Place a fiducial cut on particle in drift chamber
// Based on particle, torus bending, etc. determine which cut to place
bool CutManager::DC_fiducial_cut(part particle){
    int pid = particle.pid;
    if(pid==11) return DC_fiducial_cut_XY(particle);
    if(pid==211 || pid==-211 || pid==2212){
        if(_torusBending==-1) return DC_fiducial_cut_theta_phi(particle);
        else if(_torusBending==1) return DC_fiducial_cut_XY(particle);
        else return true;
    }
    return true;
}

// 
bool CutManager::DC_fiducial_cut_theta_phi(part particle){
    const auto minparams = ((_torusBending==-1) ? minparams_in_theta_phi : minparams_out_theta_phi);
  const auto maxparams = ((_torusBending==-1) ? maxparams_in_theta_phi : maxparams_out_theta_phi);
  double theta_DCr = 5000;
  double phi_DCr_raw = 5000;
  double x=0;
  double y=0;
  double z=0;
  for(int r = 0; r<3; r++){  
      x=0;y=0;z=0;
      switch(r){
          case 0:
              x=particle.traj_x1;
              y=particle.traj_y1;
              z=particle.traj_z1;
              break;
          case 1:
              x=particle.traj_x2;
              y=particle.traj_y2;
              z=particle.traj_z2;
              break;
          case 2:
              x=particle.traj_x3;
              y=particle.traj_y3;
              z=particle.traj_z3;
              break;
      }
  int sector = particle.sector;

  theta_DCr = 180 / PI * acos(z / sqrt(pow(x,2) + pow(y,2) + pow(z,2)));
  phi_DCr_raw = 180 / PI * atan2(y / sqrt(pow(x,2) + pow(y,2) + pow(z,2)), 
					x /sqrt(pow(x,2) + pow(y,2) + pow(z,2)));

  double phi_DCr = 5000;

  if (sector == 1) phi_DCr = phi_DCr_raw;
  if (sector == 2) phi_DCr = phi_DCr_raw - 60;
  if (sector == 3) phi_DCr = phi_DCr_raw - 120;
  if (sector == 4 && phi_DCr_raw > 0) phi_DCr = phi_DCr_raw - 180;
  if (sector == 4 && phi_DCr_raw < 0) phi_DCr = phi_DCr_raw + 180;
  if (sector == 5) phi_DCr = phi_DCr_raw + 120;
  if (sector == 6) phi_DCr = phi_DCr_raw + 60;

  int pid = 0;

  switch (particle.pid)
    {
    case 11: pid = 0; break;
    case 2212: pid = 1; break;
    case 211: pid = 2; break;
    case -211: pid = 3; break;
    case 321: pid = 4; break;
    case -321: pid = 5; break;
    default: return false; break;
    }
  
  

  double calc_phi_min = minparams[pid][sector - 1][r][0] + minparams[pid][sector - 1][r][1] * std::log(theta_DCr) 
    + minparams[pid][sector - 1][r][2] * theta_DCr + minparams[pid][sector - 1][r][3] * theta_DCr * theta_DCr;

  double calc_phi_max = maxparams[pid][sector - 1][r][0] + maxparams[pid][sector - 1][r][1] * std::log(theta_DCr)
    + maxparams[pid][sector - 1][r][2] * theta_DCr + maxparams[pid][sector - 1][r][3] * theta_DCr * theta_DCr;
  
  if(isnan(calc_phi_min)||isnan(calc_phi_max)) return false;
  if((phi_DCr < calc_phi_min) || (phi_DCr > calc_phi_max)) return false;
  }
  return true;
}

bool CutManager::DC_fiducial_cut_XY(part particle){
    
   const auto minparams = ((_torusBending==-1) ? minparams_in_XY : minparams_out_XY);
   const auto maxparams = ((_torusBending==-1) ? maxparams_in_XY : maxparams_out_XY);
   double X=0;
   double Y=0;
   for(int r = 0 ; r < 3; r++){
       X=0;
       Y=0;
       switch(r){
           case 0:
               X = particle.traj_x1;
               Y = particle.traj_y1;
               break;
           case 1:
               X = particle.traj_x2;
               Y = particle.traj_y2;
               break;
           case 2:
               X = particle.traj_x3;
               Y = particle.traj_y3;
               break;
       }

       int sector = particle.sector;

       if(sector == 2)
         {
           const double X_new = X * std::cos(-60 * PI / 180) - Y * std::sin(-60 * PI / 180);
           Y = X * std::sin(-60 * PI / 180) + Y * std::cos(-60 * PI / 180);
           X = X_new;
         }

       if(sector == 3)
         {
           const double X_new = X * std::cos(-120 * PI / 180) - Y * std::sin(-120 * PI / 180);
           Y = X * std::sin(-120 * PI / 180) + Y * std::cos(-120 * PI / 180);
           X = X_new;
         }

       if(sector == 4)
         {
           const double X_new = X * std::cos(-180 * PI / 180) - Y * std::sin(-180 * PI / 180);
           Y = X * std::sin(-180 * PI / 180) + Y * std::cos(-180 * PI / 180);
           X = X_new;
         }

       if(sector == 5)
         {
           const double X_new = X * std::cos(120 * PI / 180) - Y * std::sin(120 * PI / 180);
           Y = X * std::sin(120 * PI / 180) + Y * std::cos(120 * PI / 180);
           X = X_new;
         }

       if(sector == 6)
         {
           const double X_new = X * std::cos(60 * PI / 180) - Y * std::sin(60 * PI / 180);
           Y = X * std::sin(60 * PI / 180) + Y * std::cos(60 * PI / 180);
           X = X_new;
         }


       int pid = 0;

       switch (particle.pid)
         {
         case 11: 
           pid = 0; 
           break;
         case 2212: 
           pid = 1; 
           break;
         case 211: 
           pid = 2; 
           break;
         case -211: 
           pid = 3; 
           break;
         case 321: 
           pid = 4; 
           break;
         case -321: 
           pid = 5; 
           break;
         default: 
           return false; 
           break;
         }
       double calc_min = minparams[pid][sector - 1][r][0] + minparams[pid][sector - 1][r][1] * X;
       double calc_max = maxparams[pid][sector - 1][r][0] + maxparams[pid][sector - 1][r][1] * X;
       if(isnan(calc_min)||isnan(calc_max)) return false;
       if((Y<calc_min) || (Y>calc_max)) {  return false;}
   }
    return true;
}

// Perform cut on the chi2pid as a function of particle momentum
bool CutManager::chi2pid(part particle,int isStrict){
    bool passChargedPionChi2=false;
    int pid=particle.pid;
    double chi2=particle.chi2;
    double p=particle.p;
    if(pid==211||pid==-211){
    // Determine pion charge dependent C value
    float C = 0.0;
    (pid==211 ? C=0.88 : C=0.93);
    // 2 different pion chi2pid regions
    // standard
    // strict
    if(isStrict==0){
      if(p<2.44)
	passChargedPionChi2=chi2<C*3;
      else
	passChargedPionChi2=chi2<C*(0.00869 + 14.98587 * exp(-p/1.18236) + 1.81751 * exp(-p/4.86394));
    } 
    else if(isStrict==1){
      if(p<2.44)
	passChargedPionChi2=chi2<C*3;
      else if(p<4.6)
	passChargedPionChi2=chi2< C * (0.00869 + 14.98587 * exp(-p/1.18236) + 1.81751 * exp(-p/4.86394));
      else
	passChargedPionChi2=chi2< C * (-1.14099 + 24.14992 * exp(-p/1.36554) + 2.66876 * exp(-p/6.80552));

    }
  }
  return passChargedPionChi2;
    return true;
}


// Perform Electron Sampling Fraction cut
bool CutManager::EleSampFrac(part particle){
    double p = particle.p;
    double Ele_ECIN_e = particle.ecin_e;
    double Ele_PCAL_e = particle.pcal_e;
    double Ele_ECOUT_e = particle.ecout_e;
    int Ele_sector = particle.pcal_sector;
    
    bool sfcutDiag, sfcutSigma;
  // calorimeter diagonal cut, on PCAL and ECIN SF correlation
  if(p<4.5) sfcutDiag=true; // only applies above HTCC threshold
  else sfcutDiag = Ele_ECIN_e/p > 0.2 - Ele_PCAL_e/p; 

  // compute SF
  float eleSampFrac = (Ele_PCAL_e + Ele_ECIN_e + Ele_ECOUT_e) / p;

  // need defined EC sector to apply SF cut
  if(Ele_sector>=1 && Ele_sector<=6) {

    // parameters for SF mu(p) and sigma(p) functions ///////   from F.X. and RGA common analysis note
    Double_t sfMu[3][6];
    Double_t sfSigma[3][6];

    // From the runNumber, determine the sampling fraction quantities
    sampFracInfo(_run,sfMu,sfSigma);

    // calculate mu(p) and sigma(p), where p=`p`
    Double_t mu    = sfMu[0][Ele_sector-1] + (sfMu[1][Ele_sector-1]/1000) * pow(p-sfMu[2][Ele_sector-1],2);
    Double_t sigma = sfSigma[0][Ele_sector-1] + sfSigma[1][Ele_sector-1] / (10 * (p-sfSigma[2][Ele_sector-1]));
    
    // SF must be within 3.5 sigma of mean; here SF is from PCAL+ECAL+ECIN/p
    sfcutSigma = TMath::Abs(eleSampFrac-mu) < 3.5*sigma;

  } else {
    sfcutSigma = false;
  }

  // return full SF cut result
  return sfcutDiag && sfcutSigma;
}


// Minimum pion momentum cut
bool CutManager::minPiMomentum(part particle){
    if(particle.p < 1.25) return false;
    return true;
}


// Vz in target cell cut
bool CutManager::VzCut(part particle){
    if(particle.pid==11){
        if(particle.vz<-8||particle.vz>3) return false;
    }
    return true;
}

// Minimum energy deposited in PCAL cut
bool CutManager::minEpcal(part particle){
    int pid=particle.pid;
    if(pid==11) return particle.pcal_e>0.07;
    else if(pid==22) return particle.pcal_e>0.00;
    return true;
}

// Fiducial cut on calorimeter edges
bool CutManager::caloEdges(part particle,int strength){
    double min_lu,max_lu,min_lv,max_lv;
    if(strength==0){
        min_lu=9;
        max_lu=400;
        min_lv=9;
        max_lv=400;
    }
    else if(strength==1){
        min_lu=14;
        max_lu=400;
        min_lv=14;
        max_lv=400;
    }
    else
    {
        min_lu=9;
        max_lu=400;
        min_lv=9;
        max_lv=400;
    }
    
    if(particle.pcal_lu < min_lu || particle.pcal_lu > max_lu) return false;
    if(particle.pcal_lv < min_lv || particle.pcal_lv > max_lv) return false;
//     if(particle.ecin_lu < min_lu || particle.ecin_lu > max_lu) return false;
//     if(particle.ecin_lv < min_lv || particle.ecin_lv > max_lv) return false;
//     if(particle.ecout_lu < min_lu || particle.ecout_lu > max_lu) return false;
//     if(particle.ecout_lv < min_lv || particle.ecout_lv > max_lv) return false;
    return true;
}

bool CutManager::Ele3calo(part particle){
    if(particle.pcal_e<=0 || particle.ecin_e<=0 || particle.ecout_e<=0) return false;
    return true;
}

bool CutManager::photonMinEtot(part particle){
    if(particle.E<0.2) return false;
    return true;
}

bool CutManager::photonElectronAngle(part particle_A, part particle_B){
    TLorentzVector lv_A, lv_B;
    lv_A.SetPx(particle_A.px);
    lv_A.SetPy(particle_A.py);
    lv_A.SetPz(particle_A.pz);
    lv_A.SetE(particle_A.E);

    lv_B.SetPx(particle_B.px);
    lv_B.SetPy(particle_B.py);
    lv_B.SetPz(particle_B.pz);
    lv_B.SetE(particle_B.E);

    double angle = lv_A.Angle(lv_B.Vect());

    if (angle < 8 * TMath::DegToRad())
        return false;
    else
        return true;
}

bool CutManager::photonBetaCut(part particle){
    if(particle.beta<0.9||particle.beta>1.1) return false;
    return true;
}