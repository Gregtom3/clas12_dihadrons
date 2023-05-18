#ifndef CLAS12Analysis_H
#define CLAS12Analysis_H

#include "Constants.h"
#include "Structs.h"
#include "HipoBankInterface.h"
#include "Kinematics.h"

class CLAS12Analysis {
public:
    // Constructor and destructor
    CLAS12Analysis();
    CLAS12Analysis(const std::unique_ptr<clas12::clas12reader>&, TLorentzVector, TLorentzVector);
    CLAS12Analysis(const std::unique_ptr<clas12::clas12reader>&, double);
    ~CLAS12Analysis();
    
    // Member functions
    std::vector<part> load_reco_particles(const std::unique_ptr<clas12::clas12reader>&);
    std::vector<part> load_mc_particles(const std::unique_ptr<clas12::clas12reader>&);
    int find_reco_scattered_electron(std::vector<part>&);
    bool reco_event_contains_final_state(std::vector<part>, FS);
    
    DIS_EVENT calc_reco_event_variables(std::vector<part>);
    DIS_EVENT calc_mc_event_variables(std::vector<part>);
    
    void match_mc_to_reco(std::vector<part>&, std::vector<part>&);
    
private:
    // Private member variables
    TLorentzVector init_electron, target;
    HipoBankInterface hipoBankInterface;
    Kinematics _kin;
    double _electron_beam_energy;
    double s; // com energy
};

#endif  // CLAS12ANALYSIS_H