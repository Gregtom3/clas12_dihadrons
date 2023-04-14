#!/bin/bash
# This script sends in multiple types of jobs to calculate the BSAs for dihadrons at CLAS12
# The program uses the "merged" ROOT TTrees from the last step of the analysis pipeline
# Currently, 3 methods of calculating BSAs are used
# Each of these methods call "macros/calc_asymmetry.C" 

# 1. standard
#   This method is reserved for non-pi0 dihadrons where background corrections are not needed
#   Therefore, 


volatile=/volatile/clas12/users/gmat/clas12analysis.sidis.data
NOW=$( date '+%F_%H_%M_%S' )
farmout=/farm_out/gmat/clas12analysis.sidis.data/clas12_dihadrons/$NOW

# Define the function printred, which prints the argument with the color red
printred () {
  echo -e "\e[31m$1\e[0m"
}

# Define the function printblue, which prints the argument with the color blue
printblue () {
  echo -e "\e[34m$1\e[0m"
}

# Define the function printgreen, which prints the argument with the color green
printgreen () {
  echo -e "\e[32m$1\e[0m"
}

mkdir_green () {
    mkdir -p "$1"
    if [ -n "$2" ]; then
        printf "\t%.s" $(seq 1 "$2")
    fi
    printgreen "mkdir $1"
}

# Assign the variable PWD the value of the current working directory
PWD=`pwd`

# If the current working directory's base name is not "clas12_dihadrons"
if [ $(basename $(pwd)) != "clas12_dihadrons" ]; then
  # Print an error message
  printred "This script must be run from the clas12_dihadrons directory"
  # Exit
  exit
fi

# Prompt the user for a project name
PROJECT_DIR=$PWD/projects
echo "Available projects"
ls $PROJECT_DIR

if [[ -n "$1" ]]; then
    PROJECT_NAME=$1
else
    read -p "Please enter a project name: " PROJECT_NAME
fi

# Check if the project exists in the project directory
if [ -d "$PROJECT_DIR/$PROJECT_NAME" ]; then
  echo "Project exists"
  # Set the project directory to the selected project
  PROJECT_DIR=$PROJECT_DIR/$PROJECT_NAME
else
  # Print an error message and exit the program
  echo "Project does not exist...Aborting..."
  exit 1
fi

VOLATILE_DIR="$volatile/clas12_dihadrons/projects/$PROJECT_NAME"
DATA_DIR="$VOLATILE_DIR/data"
BRU_DIR=$VOLATILE_DIR/asym
FARMOUT_DIR="$farmout"

# Create folders within FARMOUT_DIR
mkdir_green "$FARMOUT_DIR/log"
mkdir_green "$FARMOUT_DIR/err"
mkdir_green "$FARMOUT_DIR/slurm"
mkdir_green "$FARMOUT_DIR/shell"

BINNING_DIR=$PWD/utils/binning_files
echo "Available binning schemes"
ls $BINNING_DIR
if [[ -n "$2" ]]; then
    BINNING_FILE=$2
else
    read -p "Please enter a bin file: " BINNING_FILE
fi

BINNING_FILE=$BINNING_DIR/$BINNING_FILE


# Read the contents of the YAML file into a variable
binyaml=$(cat $BINNING_FILE)

# Count the number of times ".root" appears in the YAML file
schemes=$(echo "$binyaml" | grep -o "numDimensions:" | wc -l)

CUT_TITLES=("v6" "v4" "v4" "v6" "v1" "v3")
#CUT_TITLES=("v0" "v0" "v0" "v0" "v0" "v0")

# Create pion pid pairs
pion_pairs=("piplus_piplus" "piplus_pi0" "piminus_pi0" "piminus_piminus" "pi0_pi0" "piplus_piminus")
pion_pairs=("piplus_pi0" )
# Create list of unique datasets
#datasets=("MC_RGA_inbending" "MC_RGA_outbending" "Fall2018_RGA_inbending" "Fall2018_RGA_outbending" "Spring2019_RGA_inbending")
datasets=("MC_RGA_inbending")

# First write the scripts for the fitting code without ML (i.e. without pi0's)
for ((i=0; i<${#pion_pairs[@]}; i++)); do
  pion_pair=${pion_pairs[$i]}
  
  if [[ $pion_pair == *"pi0"* ]]; then
      continue
  fi

  CUT_TITLE=${CUT_TITLES[$i]}
  for dataset in "${datasets[@]}"; do
    for ((binscheme=0; binscheme<$schemes; binscheme++)); do
        nbins=$(python $PWD/utils/read_bin_nums.py $BINNING_FILE $binscheme)
        for ((binnum=0; binnum<$nbins; binnum++)); do
        
            FILE=$DATA_DIR/$pion_pair/${dataset}_merged.root
            
            file=${pion_pair}_${dataset}_${binscheme}_${binnum}


            slurmshell=$FARMOUT_DIR/shell/brudihadron_$file.sh
            slurmslurm=$FARMOUT_DIR/slurm/brudihadron_$file.slurm

            touch -f $slurmshell
            touch -f $slurmslurm

            chmod +x $slurmshell
            cat >> $slurmslurm << EOF
#!/bin/bash
#SBATCH --account=clas12
#SBATCH --partition=production
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=job_brudihadron_$file
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --output=$FARMOUT_DIR/log/brudihadron_$file.out
#SBATCH --error=$FARMOUT_DIR/err/brudihadron_$file.err
$slurmshell
EOF

            echo "#!/bin/tcsh" >> $slurmshell
            echo "module unload root" >> $slurmshell
            echo "source /group/clas12/packages/setup.csh" >> $slurmshell
            echo "module load clas12/pro" >> $slurmshell

            echo "/u/site/12gev_phys/2.4/Linux_CentOS7.7.1908-gcc9.2.0/root/6.20.04/bin/root $BRUFIT/macros/LoadBru.C -b -q -l $PWD/macros/calc_asymmetry.C\(\\\"${FILE}\\\",\\\"${BINNING_FILE}\\\",\\\"${BRU_DIR}\\\",$binscheme,$binnum,\\\"${CUT_TITLE}\\\",0\)" >> $slurmshell

            echo "STANDARD: Submitting ${pion_pair} slurm job for $(basename "$FILE"), binning scheme $((binscheme+1)) of $((schemes)), bin $((binnum+1)) of $((nbins))"
            sbatch --quiet $slurmslurm
        done
      done
    done
done









# Next, write the scripts for the sPlot + Sideband fitting code with ML (i.e. with pi0's)
for ((i=0; i<${#pion_pairs[@]}; i++)); do
  pion_pair=${pion_pairs[$i]}
  
  if ! [[ $pion_pair == *"pi0"* ]]; then
      continue
  fi

  if [[ $pion_pair == "pi0_pi0" ]]; then
      continue
  fi
  CUT_TITLE=${CUT_TITLES[$i]}
  for dataset in "${datasets[@]}"; do
    for ((ML=0; ML<2; ML++)); do
        for ((binscheme=0; binscheme<$schemes; binscheme++)); do
            nbins=$(python $PWD/utils/read_bin_nums.py $BINNING_FILE $binscheme)
            FILE=$DATA_DIR/$pion_pair/${dataset}_merged.root

            file=splot_${pion_pair}_${dataset}_${binscheme}_isML_$ML

            slurmshell=$FARMOUT_DIR/shell/brudihadron_$file.sh
            slurmslurm=$FARMOUT_DIR/slurm/brudihadron_$file.slurm

            touch -f $slurmshell
            touch -f $slurmslurm

            chmod +x $slurmshell
            cat >> $slurmslurm << EOF
#!/bin/bash
#SBATCH --account=clas12
#SBATCH --partition=production
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=job_brudihadron_$file
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --output=$FARMOUT_DIR/log/brudihadron_$file.out
#SBATCH --error=$FARMOUT_DIR/err/brudihadron_$file.err
$slurmshell
EOF

            echo "#!/bin/tcsh" >> $slurmshell
            echo "module unload root" >> $slurmshell
            echo "source /group/clas12/packages/setup.csh" >> $slurmshell
            echo "module load clas12/pro" >> $slurmshell

            echo "/u/site/12gev_phys/2.4/Linux_CentOS7.7.1908-gcc9.2.0/root/6.20.04/bin/root $BRUFIT/macros/LoadBru.C -b -q -l $PWD/macros/calc_asymmetry.C\(\\\"${FILE}\\\",\\\"${BINNING_FILE}\\\",\\\"${BRU_DIR}\\\",$binscheme,0,\\\"${CUT_TITLE}\\\",$ML,1,0,0\)" >> $slurmshell

            echo "SPLOT (ML=$ML): Submitting ${pion_pair} slurm job for $(basename "$FILE"), binning scheme $((binscheme+1)) of $((schemes))"
            # Now write the functions that submit the individual binned fits
            for ((binnum=0; binnum<$nbins; binnum++)); do
            
                FILE=$DATA_DIR/$pion_pair/${dataset}_merged.root
                # ---------------------------------------------------------------------------------------------
                # SWEIGHTED FIT HERE
                # ---------------------------------------------------------------------------------------------
                file=sweight_${pion_pair}_${dataset}_${binscheme}_${binnum}_isML_$ML

                slurmshell2=$FARMOUT_DIR/shell/brudihadron_$file.sh
                slurmslurm2=$FARMOUT_DIR/slurm/brudihadron_$file.slurm

                touch -f $slurmshell2
                touch -f $slurmslurm2

                chmod +x $slurmshell2
                cat >> $slurmslurm2 << EOF
#!/bin/bash
#SBATCH --account=clas12
#SBATCH --partition=production
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=job_brudihadron_$file
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --output=$FARMOUT_DIR/log/brudihadron_$file.out
#SBATCH --error=$FARMOUT_DIR/err/brudihadron_$file.err
$slurmshell2
EOF
            
                echo "#!/bin/tcsh" >> $slurmshell2
                echo "module unload root" >> $slurmshell2
                echo "source /group/clas12/packages/setup.csh" >> $slurmshell2
                echo "module load clas12/pro" >> $slurmshell2

                echo "/u/site/12gev_phys/2.4/Linux_CentOS7.7.1908-gcc9.2.0/root/6.20.04/bin/root $BRUFIT/macros/LoadBru.C -b -q -l $PWD/macros/calc_asymmetry.C\(\\\"${FILE}\\\",\\\"${BINNING_FILE}\\\",\\\"${BRU_DIR}\\\",$binscheme,$binnum,\\\"${CUT_TITLE}\\\",$ML,0,1,0\)" >> $slurmshell2

                echo "sbatch $slurmslurm2" >> $slurmshell
                
                # ---------------------------------------------------------------------------------------------
                # SIDEBAND FIT HERE
                # ---------------------------------------------------------------------------------------------
                file=sideband_${pion_pair}_${dataset}_${binscheme}_${binnum}_isML_$ML


                slurmshell2=$FARMOUT_DIR/shell/brudihadron_$file.sh
                slurmslurm2=$FARMOUT_DIR/slurm/brudihadron_$file.slurm

                touch -f $slurmshell2
                touch -f $slurmslurm2

                chmod +x $slurmshell2
                cat >> $slurmslurm2 << EOF
#!/bin/bash
#SBATCH --account=clas12
#SBATCH --partition=production
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=job_brudihadron_$file
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --output=$FARMOUT_DIR/log/brudihadron_$file.out
#SBATCH --error=$FARMOUT_DIR/err/brudihadron_$file.err
$slurmshell2
EOF
            
                echo "#!/bin/tcsh" >> $slurmshell2
                echo "module unload root" >> $slurmshell2
                echo "source /group/clas12/packages/setup.csh" >> $slurmshell2
                echo "module load clas12/pro" >> $slurmshell2

                echo "/u/site/12gev_phys/2.4/Linux_CentOS7.7.1908-gcc9.2.0/root/6.20.04/bin/root $BRUFIT/macros/LoadBru.C -b -q -l $PWD/macros/calc_asymmetry.C\(\\\"${FILE}\\\",\\\"${BINNING_FILE}\\\",\\\"${BRU_DIR}\\\",$binscheme,$binnum,\\\"${CUT_TITLE}\\\",$ML,0,0,1\)" >> $slurmshell2
                
                echo "    SWEIGHT+SIDEBAND (ML=$ML): Submitting ${pion_pair} slurm job for $(basename "$FILE"), binning scheme $((binscheme+1)) of $((schemes)), bin $((binnum+1)) of $((nbins))"
                echo "sbatch $slurmslurm2" >> $slurmshell
                
            done

            echo "bash $PWD/scripts/wait_to_delete_inject_file.sh job_brudihadron_sweight" >> $slurmshell
            echo "bash $PWD/scripts/wait_to_delete_inject_file.sh job_brudihadron_sideband" >> $slurmshell
            echo "/u/site/12gev_phys/2.4/Linux_CentOS7.7.1908-gcc9.2.0/root/6.20.04/bin/root $BRUFIT/macros/LoadBru.C -b -q -l $PWD/macros/calc_asymmetry.C\(\\\"${FILE}\\\",\\\"${BINNING_FILE}\\\",\\\"${BRU_DIR}\\\",$binscheme,0,\\\"${CUT_TITLE}\\\",$ML,1,0,0,1\)" >> $slurmshell
            

            sbatch --quiet $slurmslurm
          done
        done
    done
done






# Next, write the scripts for the sideband fitting code with ML (i.e. with pi0's)
# for ((i=0; i<${#pion_pairs[@]}; i++)); do
#   pion_pair=${pion_pairs[$i]}
  
#   if ! [[ $pion_pair == *"pi0"* ]]; then
#       continue
#   fi

#   if [[ $pion_pair == "pi0_pi0" ]]; then
#       continue
#   fi
  
#   CUT_TITLE=${CUT_TITLES[$i]}
#   for dataset in "${datasets[@]}"; do
#     for ((ML=0; ML<2; ML++)); do
#         for ((binscheme=0; binscheme<$schemes; binscheme++)); do
#             nbins=$(python $PWD/utils/read_bin_nums.py $BINNING_FILE $binscheme)
#             for ((binnum=0; binnum<$nbins; binnum++)); do
#                 FILE=$DATA_DIR/$pion_pair/${dataset}_merged.root

#                 file=sideband_${pion_pair}_${dataset}_${binscheme}_${binnum}_isML_$ML


#                 slurmshell=$FARMOUT_DIR/shell/brudihadron_$file.sh
#                 slurmslurm=$FARMOUT_DIR/slurm/brudihadron_$file.slurm

#                 touch -f $slurmshell
#                 touch -f $slurmslurm

#                 chmod +x $slurmshell
#                 cat >> $slurmslurm << EOF
# #!/bin/bash
# #SBATCH --account=clas12
# #SBATCH --partition=production
# #SBATCH --mem-per-cpu=4000
# #SBATCH --job-name=job_brudihadron_$file
# #SBATCH --cpus-per-task=8
# #SBATCH --time=24:00:00
# #SBATCH --output=$FARMOUT_DIR/log/brudihadron_$file.out
# #SBATCH --error=$FARMOUT_DIR/err/brudihadron_$file.err
# $slurmshell
# EOF
            
#                 echo "#!/bin/tcsh" >> $slurmshell
#                 echo "module unload root" >> $slurmshell
#                 echo "source /group/clas12/packages/setup.csh" >> $slurmshell
#                 echo "module load clas12/pro" >> $slurmshell

#                 echo "/u/site/12gev_phys/2.4/Linux_CentOS7.7.1908-gcc9.2.0/root/6.20.04/bin/root $BRUFIT/macros/LoadBru.C -b -q -l $PWD/macros/calc_asymmetry.C\(\\\"${FILE}\\\",\\\"${BINNING_FILE}\\\",\\\"${BRU_DIR}\\\",$binscheme,$binnum,\\\"${CUT_TITLE}\\\",$ML,0,0,1\)" >> $slurmshell
#                 echo "SIDEBAND (ML=$ML): Submitting ${pion_pair} slurm job for $(basename "$FILE"), binning scheme $((binscheme+1)) of $((schemes)), bin $((binnum+1)) of $((nbins))"
#                 sbatch --quiet $slurmslurm
#             done
#           done
#         done
#     done
# done
