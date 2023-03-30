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

# Check if python3 version is 3.9.7

PYTHON_VERSION=$(python3 --version) 

if [[ "$PYTHON_VERSION" != "Python 3.9.7" ]]; then
    printred "Your python3 version is not 3.9.7."
    echo -e "Please run...\n\t module unload root\n\t module unload python\n\t module load root"
    echo "This program requires the latest version of ROOT on ifarm (6.26.10)"
    echo "Exiting"
    exit 0
fi

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
    PROJECT_NAME="$1"
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

# Create folders within FARMOUT_DIR
mkdir_green "$FARMOUT_DIR/log"
mkdir_green "$FARMOUT_DIR/err"
mkdir_green "$FARMOUT_DIR/slurm"
# --------------------------------------------------------------------------------
#               Determine which model to use for the predictions
# --------------------------------------------------------------------------------
#Function to get the basename and substring from the filename containing "_model_" 
getBasenameSubstring(){
    PROJECT_DIR=$1
    FILE=$2
    BASENAME=$(basename $FILE)
    SUBSTRING=${BASENAME%_model_*}
    echo $SUBSTRING
}

LIST=()
CALO_LIST=()
TRACK_LIST=()

for file in $(find $PROJECT_DIR -name "model_params.txt"); do
    FIRST_LINE=$(head -n1 "$file")
    if [[ ! ${LIST[*]} =~ $FIRST_LINE ]]; then
        LIST+=(${FIRST_LINE})
        CALO_LIST+=(${FIRST_LINE}_calo)
        TRACK_LIST+=(${FIRST_LINE}_track)
    fi
done

#Print out the list of unique substrings
echo "Unique substrings found:"
echo ${CALO_LIST[@]} | tr " " "\n"
echo ${TRACK_LIST[@]} | tr " " "\n"

#Prompt user to pick a substring
echo -n "Which substring would you like to use? "
if [[ -n "$2" ]]; then
    CHOICE="$2"
else
    read CHOICE
fi

# Find the index of CHOICE in either CALO_LIST or TRACK_LIST
if [[ "${CALO_LIST[*]}" =~ "${CHOICE}" ]]; then
    INDEX=$(echo ${CALO_LIST[@]} | tr " " "\n" | grep -n "${CHOICE}" | cut -d ":" -f 1)
    ARRAY="CALO_LIST"
elif [[ "${TRACK_LIST[*]}" =~ "${CHOICE}" ]]; then
    INDEX=$(echo ${TRACK_LIST[@]} | tr " " "\n" | grep -n "${CHOICE}" | cut -d ":" -f 1)
    ARRAY="TRACK_LIST"
else
    echo "Error: CHOICE not found in CALO_LIST or TRACK_LIST."
    exit 1
fi

echo "Index of CHOICE in $ARRAY: $((INDEX-1))"
echo "${LIST[$(($((INDEX-1))))]}"
FINAL_LIST=()
for file in $(find $PROJECT_DIR -name "*_model_*"); do
    if [[ $file == */${LIST[$(($((INDEX-1))))]}/* ]]; then
        FINAL_LIST+=($file)
    fi
done

# LIST=()

# for file in $(find $PROJECT_DIR -name "model_params.txt"); do
#     FIRST_LINE=$(head -n1 "$file")
#     if [[ ! ${LIST[*]} =~ $FIRST_LINE ]]; then
#         LIST+=(${SUBSTRING}_calo)
#         LIST+=(${SUBSTRING}_track)
#     fi
# done


# #Search for files in project directory containing "_model_" in their filename
# for file in $(find $PROJECT_DIR -name "*_model_*"); do
#   #Get the substring from the filename
#   SUBSTRING=$(getBasenameSubstring $PROJECT_DIR $file)
#   #Check if substring is already in list
#   if [[ ! ${LIST[*]} =~ $SUBSTRING ]]; then
#     #If not, add it to the list
#     LIST+=($SUBSTRING)
#   fi
# done

# #Print out the list of unique substrings
# echo "Unique substrings found:"
# echo ${LIST[@]}


# #Prompt user to pick a substring
# echo -n "Which substring would you like to use? "
# if [[ -n "$2" ]]; then
#     CHOICE="$2"
# else
#     read CHOICE
# fi


# #Prompt user to pick a method
# echo -n "Would you like to predict with \"calo\" or \"track\"? "
# if [[ -n "$3" ]]; then
#     METHOD="$3"
# else
#     read METHOD
# fi

# #Create sublist of files with the chosen substring 
# FINAL_LIST=()
# for file in $(find $PROJECT_DIR -name "*_model_*"); do
#   #Get the substring from the filename
#   SUBSTRING=$(getBasenameSubstring $PROJECT_DIR $file)
#   if [[ $SUBSTRING == $CHOICE ]]; then
#     #Check for calo or track substring
#     if [[ $file == *"/${METHOD}/"* ]]; then
#       FINAL_LIST+=($file)
#     fi
#   fi
# done

# Create pion pid pairs for the ML portion
pion_pairs=("piplus_pi0" "piminus_pi0" "pi0_pi0")
# Create list of unique datasets
datasets=("Fall2018_RGA_inbending" "Fall2018_RGA_outbending" "Spring2019_RGA_inbending" "MC_RGA_inbending" "MC_RGA_outbending")

# Returns the monte carlo dataset used to train the appropriate model
function dataset_to_model() {
    case "$1" in
        "Fall2018_RGA_inbending")
            echo "MC_RGA_inbending"
            ;;
        "Fall2018_RGA_outbending")
            echo "MC_RGA_outbending"
            ;;
        "Spring2019_RGA_inbending")
            echo "MC_RGA_inbending"
            ;;
        "MC_RGA_inbending")
            echo "MC_RGA_inbending"
            ;;
        "MC_RGA_outbending")
            echo "MC_RGA_outbending"
            ;;
    esac
}

# Nested for loop over pion pairs and datasets
for pion_pair in ${pion_pairs[@]}; do
    for dataset in ${datasets[@]}; do
        # Get the suffix to look for in FINAL_LIST
        SUFFIX=$(dataset_to_model $dataset)
        # Look for the corresponding element in FINAL_LIST
        for model in ${FINAL_LIST[@]}; do
            if [[ $model == *"$SUFFIX"* ]] && [[ $model == *"$pion_pair"* ]]; then
                  slurmslurm=$FARMOUT_DIR/slurm/predict_photonML_${pion_pair}_${dataset}.slurm

                  touch -f $slurmslurm

                  cat >> $slurmslurm << EOF
#!/bin/bash
#SBATCH --account=clas12
#SBATCH --partition production
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=job_photonMLpredict_${pion_pair}_${dataset}
#SBATCH --cpus-per-task=2
#SBATCH --time=24:00:00
#SBATCH --output=$FARMOUT_DIR/log/predict_photonML_${pion_pair}_${dataset}.out
#SBATCH --error=$FARMOUT_DIR/err/predict_photonML_${pion_pair}_${dataset}.err
/u/apps/python3/3.9.7/bin/python3 $PWD/machine_learning/photonID/predict.py "${DATA_DIR}/${pion_pair}" "$dataset" "$model"
EOF
                echo "Submitting slurm job for ${pion_pair} , ${dataset}"
                sbatch --quiet $slurmslurm
            fi
        done
    done
done
