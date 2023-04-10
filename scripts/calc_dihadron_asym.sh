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


# CUT_LIBRARY=$PWD/utils/cut_library.yaml
# echo "Available cut schemes (green options)"

# while IFS= read -r line; do
#   if [[ $line =~ ^[[:space:]] ]]; then
#     echo -e "$line"
#   else
#     echo -e "\e[32m$line\e[0m"
#   fi
# done < $CUT_LIBRARY


# if [[ -n "$3" ]]; then
#     CUT_TITLE=$3
# else
#     read -p "Please enter a cut style: " CUT_TITLE
# fi

CUT_TITLES=("v3" "v4" "v6" "v1" "v5" "v6")

# Create pion pid pairs
pion_pairs=("piplus_piplus" "piplus_pi0" "piminus_pi0" "piminus_piminus" "pi0_pi0" "piplus_piminus")
# Create list of unique datasets
datasets=("MC_RGA_inbending" "MC_RGA_outbending" "Fall2018_RGA_inbending" "Fall2018_RGA_outbending" "Spring2019_RGA_inbending")

for ((i=0; i<${#pion_pairs[@]}; i++)); do
  pion_pair=${pion_pairs[$i]}
  CUT_TITLE=${CUT_TITLES[$i]}
  for dataset in "${datasets[@]}"; do
    for ((binnum=0; binnum<$schemes; binnum++)); do

        FILE=$DATA_DIR/$pion_pair/${dataset}_merged.root
        file=${pion_pair}_${dataset}_${binnum}
        
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
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --output=$FARMOUT_DIR/log/brudihadron_$file.out
#SBATCH --error=$FARMOUT_DIR/err/brudihadron_$file.err
$slurmshell
EOF

        echo "#!/bin/tcsh" >> $slurmshell
        echo "module unload root" >> $slurmshell
        echo "source /group/clas12/packages/setup.csh" >> $slurmshell
        echo "module load clas12/pro" >> $slurmshell

        echo "/u/site/12gev_phys/2.4/Linux_CentOS7.7.1908-gcc9.2.0/root/6.20.04/bin/root $BRUFIT/macros/LoadBru.C -b -q -l $PWD/macros/calc_asymmetry.C\(\\\"${FILE}\\\",\\\"${BINNING_FILE}\\\",\\\"${BRU_DIR}\\\",$binnum,\\\"${CUT_TITLE}\\\",0\)" >> $slurmshell
        
        # Needs ML?
        if [[ $pion_pair == *"pi0"* ]]; then
            echo "/u/site/12gev_phys/2.4/Linux_CentOS7.7.1908-gcc9.2.0/root/6.20.04/bin/root $BRUFIT/macros/LoadBru.C -b -q -l $PWD/macros/calc_asymmetry.C\(\\\"${FILE}\\\",\\\"${BINNING_FILE}\\\",\\\"${BRU_DIR}\\\",$binnum,\\\"${CUT_TITLE}\\\",1\)" >> $slurmshell
        fi

        echo "Submitting slurm job for $(basename "$FILE"), binning scheme $((binnum+1)) of $((schemes))"
        sbatch --quiet $slurmslurm
    done
  done
done
