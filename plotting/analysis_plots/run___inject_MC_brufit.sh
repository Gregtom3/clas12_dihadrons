#!/bin/bash


# Bins for different binning schemes
# They are defined as strings where each number is separated by a comma
z_bins="0,0.35,0.43,0.49,0.550,0.620,0.700,0.8,1"
Mh_bins="0.3,0.48,0.6,0.7,0.85,1.0,1.15,1.3,2"
x_bins="0,0.1,0.13,0.16,0.19,0.235,0.3,0.5,1"
pTtot_bins="0,0.15,0.25,0.35,0.45,0.55,0.7,0.9,2"
Mx_bins="-0.5,0.2,0.6,1,1.6,2,2.4,2.8,3.4,4"

# Array of binning schemes and corresponding bins
declare -A bins
bins=(["z"]=$z_bins ["Mh"]=$Mh_bins ["x"]=$x_bins ["pTtot"]=$pTtot_bins) #["Mx"]=$Mx_bins)

# Define asymmetry functions for 'all_zeros', 'all_simple', and 'all_simple_tight_sideband' cases
asymmetry_functions_zeros='["0." ,"0.", "sin(:phi_R0:)"],
                            ["0." ,"0.", "sin(:phi_h:)"],
                            ["0." ,"0.", "sin(-:phi_h:+2*:phi_R0:)"],
                            ["0." ,"0.", "sin(:phi_h:-:phi_R0:)"],
                            ["0." ,"0.", "sin(2*:phi_h:-2*:phi_R0:)"],
                            ["0." ,"0.", "sin(2*:phi_h:-:phi_R0:)"],
                            ["0." ,"0.", "sin(3*:phi_h:-2*:phi_R0:)"]'

asymmetry_functions_simple='["0.05" ,"-0.02", "sin(:phi_R0:)"],
                             ["0.05" ,"-0.02", "sin(:phi_h:)"],
                             ["0.05" ,"-0.02", "sin(-:phi_h:+2*:phi_R0:)"],
                             ["0.05" ,"-0.02", "sin(:phi_h:-:phi_R0:)"],
                             ["0.05" ,"-0.02", "sin(2*:phi_h:-2*:phi_R0:)"],
                             ["0.05" ,"-0.02", "sin(2*:phi_h:-:phi_R0:)"],
                             ["0.05" ,"-0.02", "sin(3*:phi_h:-2*:phi_R0:)"]'
                             
asymmetry_functions_simple_phiR_only='["0.05" ,"-0.02", "sin(:phi_R0:)"],
                                      ["0." ,"0.", "sin(:phi_h:)"],
                                      ["0." ,"0.", "sin(-:phi_h:+2*:phi_R0:)"],
                                      ["0." ,"0.", "sin(:phi_h:-:phi_R0:)"],
                                      ["0." ,"0.", "sin(2*:phi_h:-2*:phi_R0:)"],
                                      ["0." ,"0.", "sin(2*:phi_h:-:phi_R0:)"],
                                      ["0." ,"0.", "sin(3*:phi_h:-2*:phi_R0:)"]'
                             
asymmetry_functions_Mh_incline='["0.05+0.05*:Mh:" ,"-0.02", "sin(:phi_R0:)"],
                                ["0.05+0.05*:Mh:" ,"-0.02", "sin(:phi_h:)"],
                                ["0.05+0.05*:Mh:" ,"-0.02", "sin(-:phi_h:+2*:phi_R0:)"],
                                ["0.05+0.05*:Mh:" ,"-0.02", "sin(:phi_h:-:phi_R0:)"],
                                ["0.05+0.05*:Mh:" ,"-0.02", "sin(2*:phi_h:-2*:phi_R0:)"],
                                ["0.05+0.05*:Mh:" ,"-0.02", "sin(2*:phi_h:-:phi_R0:)"],
                                ["0.05+0.05*:Mh:" ,"-0.02", "sin(3*:phi_h:-2*:phi_R0:)"]'
                                
asymmetry_functions_M2_incline='["0.05" ,"-0.02+0.05*:M2:", "sin(:phi_R0:)"],
                                ["0.05" ,"-0.02+0.05*:M2:", "sin(:phi_h:)"],
                                ["0.05" ,"-0.02+0.05*:M2:", "sin(-:phi_h:+2*:phi_R0:)"],
                                ["0.05" ,"-0.02+0.05*:M2:", "sin(:phi_h:-:phi_R0:)"],
                                ["0.05" ,"-0.02+0.05*:M2:", "sin(2*:phi_h:-2*:phi_R0:)"],
                                ["0.05" ,"-0.02+0.05*:M2:", "sin(2*:phi_h:-:phi_R0:)"],
                                ["0.05" ,"-0.02+0.05*:M2:", "sin(3*:phi_h:-2*:phi_R0:)"]'

asymmetry_functions_M2_incline_phiR_only='["0.03" ,"-0.03+0.067*:M2:", "sin(:phi_R0:)"],
                                          ["0." ,"0.", "sin(:phi_h:)"],
                                          ["0." ,"0.", "sin(-:phi_h:+2*:phi_R0:)"],
                                          ["0." ,"0.", "sin(:phi_h:-:phi_R0:)"],
                                          ["0." ,"0.", "sin(2*:phi_h:-2*:phi_R0:)"],
                                          ["0." ,"0.", "sin(2*:phi_h:-:phi_R0:)"],
                                          ["0." ,"0.", "sin(3*:phi_h:-2*:phi_R0:)"]'
# Loop over different binning schemes
for binning_scheme in "${!bins[@]}"; do
    
    # Define filenames based on binning scheme
    FILENAME_ZEROS="bruinjection_files/all_zeros_${binning_scheme}_binned.txt"
    FILENAME_SIMPLE="bruinjection_files/all_simple_${binning_scheme}_binned.txt"
    FILENAME_SIMPLE_PHIR_ONLY="bruinjection_files/phiR_only_simple_${binning_scheme}_binned.txt"
    FILENAME_MH_INCLINE="bruinjection_files/all_Mh_incline_${binning_scheme}_binned.txt"
    FILENAME_M2_INCLINE="bruinjection_files/all_M2_incline_${binning_scheme}_binned.txt"
    FILENAME_M2_INCLINE_PHIR_ONLY="bruinjection_files/M2_incline_phiR_only_${binning_scheme}_binned.txt"
    FILENAME_SELF_INCLINE="bruinjection_files/all_${binning_scheme}_incline_${binning_scheme}_binned.txt"
    
    # Clear the existing content of the file
    > $FILENAME_ZEROS
    > $FILENAME_SIMPLE
    > $FILENAME_SIMPLE_PHIR_ONLY
    > $FILENAME_MH_INCLINE
    > $FILENAME_M2_INCLINE
    > $FILENAME_M2_INCLINE_PHIR_ONLY
    > $FILENAME_SELF_INCLINE
    
    # Write new content to the files

    # All zeros case
    echo "self.injector.set_beam_polarization(1)" >> $FILENAME_ZEROS
    echo "self.injector.load_weight_funcs([$asymmetry_functions_zeros])" >> $FILENAME_ZEROS
    echo "self.injector.create_bins(\"$binning_scheme\", [${bins[$binning_scheme]}])" >> $FILENAME_ZEROS

    # All simple case
    echo "self.injector.set_beam_polarization(1)" >> $FILENAME_SIMPLE
    echo "self.injector.load_weight_funcs([$asymmetry_functions_simple])" >> $FILENAME_SIMPLE
    echo "self.injector.create_bins(\"$binning_scheme\", [${bins[$binning_scheme]}])" >> $FILENAME_SIMPLE

    # phiR only simple case
    echo "self.injector.set_beam_polarization(1)" >> $FILENAME_SIMPLE_PHIR_ONLY
    echo "self.injector.load_weight_funcs([$asymmetry_functions_simple_phiR_only])" >> $FILENAME_SIMPLE_PHIR_ONLY
    echo "self.injector.create_bins(\"$binning_scheme\", [${bins[$binning_scheme]}])" >> $FILENAME_SIMPLE_PHIR_ONLY
    
    # All Mh incline
    echo "self.injector.set_beam_polarization(1)" >> $FILENAME_MH_INCLINE
    echo "self.injector.load_weight_funcs([$asymmetry_functions_Mh_incline])" >> $FILENAME_MH_INCLINE
    echo "self.injector.create_bins(\"$binning_scheme\", [${bins[$binning_scheme]}])" >> $FILENAME_MH_INCLINE
    
    # All M2 incline
    echo "self.injector.set_beam_polarization(1)" >> $FILENAME_M2_INCLINE
    echo "self.injector.load_weight_funcs([$asymmetry_functions_M2_incline])" >> $FILENAME_M2_INCLINE
    echo "self.injector.create_bins(\"$binning_scheme\", [${bins[$binning_scheme]}])" >> $FILENAME_M2_INCLINE
    
    # PhiR only M2 incline
    echo "self.injector.set_beam_polarization(1)" >> $FILENAME_M2_INCLINE_PHIR_ONLY
    echo "self.injector.load_weight_funcs([$asymmetry_functions_M2_incline_phiR_only])" >> $FILENAME_M2_INCLINE_PHIR_ONLY
    echo "self.injector.create_bins(\"$binning_scheme\", [${bins[$binning_scheme]}])" >> $FILENAME_M2_INCLINE_PHIR_ONLY
    
    # Self binning incline
    asymmetry_functions_self_incline='["0.05+0.05*:${binning_scheme}:" ,"-0.02", "sin(:phi_R0:)"],
                                      ["0.05+0.05*:${binning_scheme}:" ,"-0.02", "sin(:phi_h:)"],
                                      ["0.05+0.05*:${binning_scheme}:" ,"-0.02", "sin(-:phi_h:+2*:phi_R0:)"],
                                      ["0.05+0.05*:${binning_scheme}:" ,"-0.02", "sin(:phi_h:-:phi_R0:)"],
                                      ["0.05+0.05*:${binning_scheme}:" ,"-0.02", "sin(2*:phi_h:-2*:phi_R0:)"],
                                      ["0.05+0.05*:${binning_scheme}:" ,"-0.02", "sin(2*:phi_h:-:phi_R0:)"],
                                      ["0.05+0.05*:${binning_scheme}:" ,"-0.02", "sin(3*:phi_h:-2*:phi_R0:)"]'
    echo "self.injector.set_beam_polarization(1)" >> $FILENAME_SELF_INCLINE
    echo "self.injector.load_weight_funcs([$asymmetry_functions_self_incline])" >> $FILENAME_SELF_INCLINE
    echo "self.injector.create_bins(\"$binning_scheme\", [${bins[$binning_scheme]}])" >> $FILENAME_SELF_INCLINE

done

# Define path to directory with injection files
INJECTION_DIR="./bruinjection_files"

# Loop over all files in the injection directory
for FILE in ${INJECTION_DIR}/*
do
  # Formulate the command string
  CMD="python3 asymmetry___inject_MC.py --infile \"/work/clas12/users/gmat/clas12/clas12_dihadrons/projects/pipi0_paper_RGA_only/volatile/data/piplus_pi0/MC_RGA_3051_0.root\" --project_name \"pipi0_paper_RGA_only\" --n_trials 5 --n_cpus 10 --program \"${FILE}\""

  # Create a SLURM batch script
  BATCH_SCRIPT="sbatch.sh"

  # Write SLURM options and command to the batch script
  echo "#!/bin/bash
#SBATCH --job-name=jobname
#SBATCH --output=/farm_out/%u/%x-%j-%N.out
#SBATCH --error=/farm_out/%u/%x-%j-%N.err
#SBATCH --partition=production
#SBATCH --account=clas12
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=2000
#SBATCH --time=24:00:00

${CMD}
" > ${BATCH_SCRIPT}

  # Submit the batch script to SLURM
  sbatch ${BATCH_SCRIPT}

  # Remove the batch script
  rm ${BATCH_SCRIPT}
  exit
done
