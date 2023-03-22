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

read -p "Please enter a project name: " PROJECT_NAME

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
DATA_DIR="$VOLATILE_DIR/data/raw"
FARMOUT_DIR="$farmout"

# Create pion pid pairs for the ML portion
pion_pairs=("piplus_pi0" "piminus_pi0" "pi0_pi0")
declare -A pion_pairs_pids
pion_pairs_pids[0,0]=211
pion_pairs_pids[0,1]=111

pion_pairs_pids[1,0]=-211
pion_pairs_pids[1,1]=111

pion_pairs_pids[2,0]=111
pion_pairs_pids[2,1]=111

# Create folders within FARMOUT_DIR
mkdir_green "$FARMOUT_DIR/log"
mkdir_green "$FARMOUT_DIR/err"
mkdir_green "$FARMOUT_DIR/slurm"
mkdir_green "$FARMOUT_DIR/shell"

# Determine the file names
files=$(ls $DATA_DIR/${pion_pairs[0]} | sort -V)

printblue "Analyzing $(echo $files | wc -w)"
for file in $files; do
    slurmshell=$FARMOUT_DIR/shell/photonML_$file.sh
    slurmslurm=$FARMOUT_DIR/slurm/photonML_$file.slurm

    touch -f $slurmshell
    touch -f $slurmslurm

    chmod +x $slurmshell
    cat >> $slurmslurm << EOF
#!/bin/bash
#SBATCH --account=clas12
#SBATCH --partition=production
#SBATCH --mem-per-cpu=2000
#SBATCH --job-name=job_photonML_$ana_$runNumber
#SBATCH --cpus-per-task=2
#SBATCH --time=24:00:00
#SBATCH --output=$FARMOUT_DIR/log/photonML_$ana_$runNumber.out
#SBATCH --error=$FARMOUT_DIR/err/photonML_$ana_$runNumber.err
$slurmshell
EOF

    echo "#!/bin/tcsh" >> $slurmshell
    echo "source /group/clas12/packages/setup.csh" >> $slurmshell
    echo "module load clas12/pro" >> $slurmshell

    # For loop over each dihadron pair
    j=0
    for pair in "${pion_pairs[@]}"; do
        FILE="$DATA_DIR/$pair/$file"
        echo "clas12root -b -q $PWD/macros/photonML.C\(\\\"${FILE}\\\",0\)" >> $slurmshell
        echo "clas12root -b -q $PWD/macros/photonML.C\(\\\"${FILE}\\\",1\)" >> $slurmshell
        j=$((j+1))
    done
    
    echo "Submitting slurm job for $file"
    sbatch --quiet $slurmslurm
done
