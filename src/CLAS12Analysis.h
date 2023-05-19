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
    void set_run_config(const std::unique_ptr<clas12::clas12reader>&);
    void set_beams(TLorentzVector, TLorentzVector);
    std::vector<part> load_reco_particles(const std::unique_ptr<clas12::clas12reader>&);
    std::vector<part> load_mc_particles(const std::unique_ptr<clas12::clas12reader>&);
    int find_reco_scattered_electron(std::vector<part>&);
    bool reco_event_contains_final_state(std::vector<part>, FS);
    
    EVENT_INFO get_event_info(const std::unique_ptr<clas12::clas12reader>&);
    EVENT calc_reco_event_variables(std::vector<part>);
    EVENT calc_mc_event_variables(std::vector<part>);
    
    void match_mc_to_reco(std::vector<part>&, std::vector<part>&);
    
    std::vector<std::vector<int>> dihadron_idxs(int,int,int[],int);
    
protected:
    // Dihadron indexing code
    void generate_combinations(std::vector<int>& input, int num, int start_idx, std::vector<int>& curr_combination, std::vector<std::vector<int>>& result);
    std::vector<std::vector<int>> unique_combinations(std::vector<int> input, int num);
    std::vector<std::vector<int>> remove_duplicates(std::vector<std::vector<int>> input);

private:
    // Private member variables
    TLorentzVector init_electron, target;
    HipoBankInterface hipoBankInterface;
    Kinematics _kin;
    double _electron_beam_energy;
    double s; // com energy
    int _idx_RUNconfig;
    int _irun;
    int _ievnum;
    int _itorus;
};

#endif  // CLAS12ANALYSIS_H
