#!/bin/bash
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
FARMOUT_DIR="$farmout"


# Create pion pid pairs
pion_pairs=("piplus_piplus" "piplus_pi0" "piminus_pi0" "piminus_piminus" "pi0_pi0" "piplus_piminus")
declare -A pion_pairs_pids
pion_pairs_pids[0,0]=211
pion_pairs_pids[0,1]=211

pion_pairs_pids[1,0]=211
pion_pairs_pids[1,1]=111

pion_pairs_pids[2,0]=-211
pion_pairs_pids[2,1]=111

pion_pairs_pids[3,0]=-211
pion_pairs_pids[3,1]=-211

pion_pairs_pids[4,0]=111
pion_pairs_pids[4,1]=111

pion_pairs_pids[5,0]=211
pion_pairs_pids[5,1]=-211

# Create folders within FARMOUT_DIR
mkdir_green "$FARMOUT_DIR/log"
mkdir_green "$FARMOUT_DIR/err"
mkdir_green "$FARMOUT_DIR/slurm"
mkdir_green "$FARMOUT_DIR/shell"


# Create list of unique datasets
datasets=("Fall2018_RGA_inbending" "Fall2018_RGA_outbending" "Spring2019_RGA_inbending" "MC_RGA_inbending" "MC_RGA_outbending" "Data_RGC" "MC_RGC")

# Loop over pion pairs
for pion_pair in "${pion_pairs[@]}"; do
    # Loop over datasets
    for dataset in "${datasets[@]}"; do
        
        rootdir=${DATA_DIR}/${pion_pair}
        
        slurmshell=$FARMOUT_DIR/shell/merge_${pion_pair}_${dataset}.sh
        slurmslurm=$FARMOUT_DIR/slurm/merge_${pion_pair}_${dataset}.slurm

        touch -f $slurmshell
        touch -f $slurmslurm

        chmod +x $slurmshell
        cat >> $slurmslurm << EOF
#!/bin/bash
#SBATCH --account=clas12
#SBATCH --partition=production
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=job_merge_${pion_pair}_${dataset}
#SBATCH --cpus-per-task=2
#SBATCH --time=24:00:00
#SBATCH --output=$FARMOUT_DIR/log/merge_${pion_pair}_${dataset}.out
#SBATCH --error=$FARMOUT_DIR/err/merge_${pion_pair}_${dataset}.err
$slurmshell
EOF

        echo "#!/bin/tcsh" >> $slurmshell
        echo "source /group/clas12/packages/setup.csh" >> $slurmshell
        echo "module load clas12/pro" >> $slurmshell

        echo "clas12root -b -q $PWD/macros/merge_dihadrons.C\(\\\"${rootdir}\\\",\\\"${dataset}\\\"\)" >> $slurmshell
                
        echo "Submitting slurm job for $pion_pair $dataset"
        sbatch --quiet $slurmslurm
    done
done
