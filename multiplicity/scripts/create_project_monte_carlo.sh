#!/bin/bash
# Assign the variable PWD the value of the current working directory
PWD=`pwd`

NOW=$( date '+%F_%H_%M_%S' )
FARMOUT_DIR=/farm_out/gmat/clas12analysis.sidis.data/clas12_dihadrons/$NOW

# Check if the current directory is "multiplicity"
if [[ ! $(basename "$(pwd)") == "multiplicity" ]]; then
    echo "Please navigate to the 'multiplicity' directory before running this script."
    exit 1
fi

# Print the names of existing projects with their creation dates
echo -e "Existing projects:\n------------------"
for project_path in ./projects/*/; do
    project=$(basename "$project_path")
    creation_date=$(stat -c %y "$project_path" | awk -F" " '{print $1}')
    echo -e "$project\t\t$creation_date"
done

# Prompt the user for a project name
read -p "Enter the new project name: " project_name

# Check if the project directory already exists
project_directory="./projects/$project_name"
# Create a link to the volatile directory within the project directory
volatile_directory="/volatile/clas12/users/gmat/clas12analysis.sidis.data/dihadron_multiplicities_2023/$project_name"

if [ -d "$project_directory" ]; then
    read -p "The project '$project_name' already exists. Do you want to overwrite it? (y/n): " overwrite_choice
    if [ "$overwrite_choice" != "y" ]; then
        echo "Operation cancelled. The project was not created."
        exit 0
    else
        rm -r $project_directory
        rm -r $volatile_directory
    fi
fi

# Create pion pid pairs
pion_pairs=("piplus_pi0" "piplus_piplus" "piminus_pi0" "piminus_piminus" "pi0_pi0" "piplus_piminus")
declare -A pion_pairs_pids
pion_pairs_pids[0,0]=211
pion_pairs_pids[0,1]=111

pion_pairs_pids[1,0]=211
pion_pairs_pids[1,1]=211

pion_pairs_pids[2,0]=-211
pion_pairs_pids[2,1]=111

pion_pairs_pids[3,0]=-211
pion_pairs_pids[3,1]=-211

pion_pairs_pids[4,0]=111
pion_pairs_pids[4,1]=111

pion_pairs_pids[5,0]=211
pion_pairs_pids[5,1]=-211

# Create subdirectories for each pion pair in the data directory
for pair in "${pion_pairs[@]}"; do
    mkdir -p "$volatile_directory/data/$pair"
done

# Create the project directory
mkdir -p "$project_directory"
mkdir -p "$volatile_directory"

# Create folders within FARMOUT_DIR
mkdir -p "$FARMOUT_DIR/log"
mkdir -p "$FARMOUT_DIR/err"
mkdir -p "$FARMOUT_DIR/slurm"
mkdir -p "$FARMOUT_DIR/shell"

# Create links to the farmout and volatile directorys
ln -s "$volatile_directory" "$project_directory/volatile"
ln -s "$FARMOUT_DIR" "$project_directory/farmout"

# Declare the hipo files from Monte Carlo
declare -a hipodirs=("/cache/clas12/rg-a/production/montecarlo/clasdis/fall2018/torus-1/v1/bkg45nA_10604MeV/")

# Prompt the user for the number of files to read
read -p "Enter the number of files you would like to read: " file_count

# Prompt the user for the value of maxEvents
read -p "Enter the value of maxEvents: " maxEvents

# Print the specified number of files from the directory
for dir in "${hipodirs[@]}"; do
    echo "Directory: $dir"
    files=("$dir"*)
    for ((i=0; i<$file_count && i<${#files[@]}; i++)); do
        file_path=${files[i]}
        read runNumber <<< $(basename $file_path | grep -oP '(?<=_job_).*(?=.hipo)')
        
        hipo_file=$file_path
        root_file="MC_RGA_${runNumber}.root"
        
        beamE=10.6041
        pid_h1=211
        pid_h2=-211
        
        slurmshell=$FARMOUT_DIR/shell/hipo2tree_${runNumber}.sh
        slurmslurm=$FARMOUT_DIR/slurm/hipo2tree_${runNumber}.slurm

        touch -f $slurmshell
        touch -f $slurmslurm

        chmod +x $slurmshell
        cat >> $slurmslurm << EOF
#!/bin/bash
#SBATCH --account=clas12
#SBATCH --partition=production
#SBATCH --mem-per-cpu=2000
#SBATCH --job-name=job_hipo2tree_$runNumber
#SBATCH --cpus-per-task=2
#SBATCH --time=24:00:00
#SBATCH --output=$FARMOUT_DIR/log/hipo2tree_$runNumber.out
#SBATCH --error=$FARMOUT_DIR/err/hipo2tree_$runNumber.err
$slurmshell
EOF
        j=0
        for pair in "${pion_pairs[@]}"; do
            if [[ "$pair" != "piplus_piminus" ]]; then
		j=$((j+1))
                continue
            fi
            outfile="$volatile_directory/data/$pair/$root_file"
            pid1=${pion_pairs_pids[${j},0]}
            pid2=${pion_pairs_pids[${j},1]}
            echo "clas12root -b -q $PWD/macros/hipo2tree_multiplicity_mc.C\(\\\"${hipo_file}\\\",\\\"${outfile}\\\",$beamE,$pid1,$pid2,$maxEvents\)" >> $slurmshell
            j=$((j+1))
        done
        
        echo "Submitting slurm job for $hipo_file"
        sbatch --quiet $slurmslurm
    done
done

echo "Project created successfully!"
