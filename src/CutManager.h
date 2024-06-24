#ifndef CutManager_h
#define CutManager_h
#include "Constants.h"

enum RUN_PERIOD {
    RGA,
    RGB,
    RGC,
    RGK
};

class CutManager{
 public:
  // Constructors
  CutManager();
  CutManager(int); // Run Number
  // Public member variables and functions
  void set_run(int);
  void set_run_period(std::string);
  RUN_PERIOD get_run_period();
  std::vector<part> filter_particles(std::vector<part>);
    
 protected:
  // Protected member functions
  bool apply_electron_cuts(part);
  bool apply_pion_cuts(part,part);
  bool apply_photon_cuts(part,part);
  bool DC_fiducial_cut(part);
  bool DC_fiducial_cut_theta_phi(part);
  bool DC_fiducial_cut_XY(part);
  bool chi2pid(part,int);
  bool EleSampFrac(part);
  bool minPiMomentum(part);
  bool hadronStatus(part);
  bool VzCut(part);
  bool minEpcal(part);
  bool caloEdges(part,int);
  bool Ele3calo(part);
  bool minEleMomentum(part);
  bool photonMinEtot(part);
  bool photonElectronAngle(part,part);
  bool photonBetaCut(part);
  bool vz_e_pi(part,part);
 private:
        
  // Private member variables
  int _run=0;
  int _torusBending=0;
  bool _is_MC=false;
  RUN_PERIOD _run_period=RGA;
        
};
#endif
