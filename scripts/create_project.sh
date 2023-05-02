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


rmdir_red () {
    if [ -n "$2" ]; then
        printf "\t%.s" $(seq 1 "$2")
    fi
    printred "Waiting 10 seconds before...rm -r $1"
    sleep 10
    rm -r $1
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
if [[ -n "$2" ]]; then
    PROJECT_NAME=$2
else
    read -p "Please enter a project name: " PROJECT_NAME
fi

# Check if the -o flag is not used
if [[ $* != *-o* ]]; then
    if [ -d "$PWD/projects/$PROJECT_NAME" ]; then
      # If it does, exit
      printred "$PROJECT_NAME already exists. Use -o to overwrite it"
      exit
    fi
fi

if [[ -z "$PROJECT_NAME" ]]; then
  echo "Error: PROJECT_NAME is empty. Please provide a project name."
  exit 1
fi

# Make the directory and call the location to the directory PROJECT_DIR
PROJECT_DIR="$PWD/projects/$PROJECT_NAME"
VOLATILE_DIR="$volatile/clas12_dihadrons/projects/$PROJECT_NAME"

if [ -d "$PROJECT_DIR" ]; then
    rmdir_red "$PROJECT_DIR"
fi

if [ -d "$VOLATILE_DIR" ]; then
    rmdir_red "$VOLATILE_DIR"
fi

FARMOUT_DIR="$farmout"
mkdir_green "$PROJECT_DIR"
mkdir_green "$VOLATILE_DIR"
mkdir_green "$FARMOUT_DIR"

# Create folders within PROJECT_DIR
mkdir_green "$PROJECT_DIR/models"
mkdir_green "$PROJECT_DIR/plots"


# Create folders within VOLATILE_DIR
mkdir_green "$VOLATILE_DIR/data"
mkdir_green "$VOLATILE_DIR/asym"

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
    mkdir_green "$VOLATILE_DIR/data/$pair" 1
done

# Create folders within FARMOUT_DIR
mkdir_green "$FARMOUT_DIR/log"
mkdir_green "$FARMOUT_DIR/err"
mkdir_green "$FARMOUT_DIR/slurm"
mkdir_green "$FARMOUT_DIR/shell"

# Create links to volatile and farmout
ln -snf $VOLATILE_DIR $PROJECT_DIR/volatile
ln -snf $FARMOUT_DIR $PROJECT_DIR/farmout



# Take the user's input
if [[ -n "$3" ]]; then
    nFiles=$3
else
    read -p "Please enter the number of files: " nFiles
fi

# Take the user's input
if [[ -n "$4" ]]; then
    nEvents=$4
else
    read -p "Please enter the number of events per file: " nEvents
fi

# Define a function that returns the hipo files for analysis
get_hipo_dirs()
{
    #declare -a hipofiles
    local version=$1
    local rungroup=$2
    if [ $rungroup == "rg-a" ]; then
        if [ $version == "MC" ]; then
            declare -a hipodirs=("/cache/clas12/$rungroup/production/montecarlo/clasdis/fall2018/torus-1/v1/bkg45nA_10604MeV/" "/cache/clas12/$rungroup/production/montecarlo/clasdis/fall2018/torus+1/v1/bkg50nA_10604MeV/")
        elif [ $version == "nSidis" ]; then
            declare -a hipodirs=("/cache/clas12/$rungroup/production/recon/fall2018/torus-1/pass1/v1/dst/train/nSidis/" "/cache/clas12/$rungroup/production/recon/fall2018/torus+1/pass1/v1/dst/train/nSidis/" "/cache/clas12/$rungroup/production/recon/spring2019/torus-1/pass1/v1/dst/train/nSidis/")
        fi
    fi
    if [ $rungroup == "rg-b" ]; then
        if [ $version == "MC" ]; then
            declare -a hipodirs=("/cache/clas12/rg-a/production/montecarlo/clasdis/fall2018/torus-1/v1/bkg45nA_10604MeV/" "/cache/clas12/rg-a/production/montecarlo/clasdis/fall2018/torus+1/v1/bkg50nA_10604MeV/")
        elif [ $version == "sidisdvcs" ]; then
            declare -a hipodirs=("/cache/clas12/rg-b/production/recon/spring2019/torus-1/pass1/v0/dst/train/sidisdvcs/" "/cache/clas12/rg-b/production/recon/fall2019/torus+1/pass1/v1/dst/train/sidisdvcs/" "/cache/clas12/rg-b/production/recon/spring2020/torus-1/pass1/v1/dst/train/sidisdvcs/")
        fi
    fi
    if [ $rungroup == "rg-c" ]; then
        if [ $version == "MC" ]; then
            declare -a hipodirs=("/work/cebaf24gev/sidis/reconstructed/polarized-plus-10.5GeV-proton/hipo/") #"/work/cebaf24gev/sidis/reconstructed/polarized-plus-10.5GeV-neutron/hipo/")
        elif [ $version == "sidisdvcs" ]; then
            declare -a hipodirs=("/volatile/clas12/$rungroup/production/dst/8.7.0_TBT/dst/train/sidisdvcs/")
        fi
    fi
    echo ${hipodirs[@]}
}


# Define a function that gets the beamE and runNumber from the hipofile
function hipo_beamE_runNumber {
    #Read the input
    local ana=$2
    local hipo=$1
    local rungroup=$3
    
    #determine runNumber
    if [ $rungroup == "rg-a" ]; then
        if [ $ana == "MC" ]; then
            read runNumber <<< $(basename $hipo | grep -oP '(?<=_job_).*(?=.hipo)')
        else
            read runNumber <<< $(basename $hipo | grep -oP '(?<='$ana'_00).*(?=.hipo)')
        fi
    elif [ $rungroup == "rg-b" ]; then
        if [ $ana == "MC" ]; then
            read runNumber <<< $(basename $hipo | grep -oP '(?<=_job_).*(?=.hipo)')
        else
            read runNumber <<< $(basename $hipo | sed -E 's/sidisdvcs_0*([0-9]+)\.hipo/\1/')
        fi
    elif [ $rungroup == "rg-c" ]; then
        if [ $ana == "MC" ]; then
            read runNumber <<< $(basename $hipo | grep -oP '\d+(?=.hipo)')
        else
            read runNumber <<< $(basename $hipo | grep -oP '(?<='$ana'_0).*(?=.hipo)')
        fi
    fi
    
    #determine beamE
    if [ $ana != "MC" ]; then
        if (($runNumber >= 5032 && runNumber <= 5666)); then
            beamE=10.6041
        elif (($runNumber >= 6616 && runNumber <= 6783)); then
            beamE=10.1998
        elif (($runNumber >= 6120 && runNumber <= 6399)); then
            beamE=10.5986
        elif (($runNumber >= 6409 && runNumber <= 6604)); then
            beamE=10.1998
        elif (($runNumber >= 11093 && runNumber <= 11283)); then
            beamE=10.4096
        elif (($runNumber >= 11284 && runNumber <= 11300)); then
            beamE=4.17179
        elif (($runNumber >= 11323 && runNumber <= 11571)); then
            beamE=10.3894
        else # Assume RGC
            beamE=10.559
        fi
    else
        if [[ $hipo=="*rg-a*/fall2018/*" ]]; then
            beamE=10.604
        fi
        if [[ $hipo==*"10.5GeV"* ]]; then
            beamE=10.5
        fi
        beamE=10.604
    fi
    #Return the two values
    echo "$runNumber $beamE"
}

#Call the function

rungroups=("rg-b")
#rungroups=("rg-a")
#rungroups=("rg-c")
versions=("data" "MC")
#versions=("MC")

for ana in "${versions[@]}"
do
    printblue "VERSION=$ana"
    for rungroup in "${rungroups[@]}"
    do
        if [ $rungroup == "rg-a" ]; then
            rg="RGA"
        elif [ $rungroup == "rg-b" ]; then
            rg="RGB"
        elif [ $rungroup == "rg-c" ]; then
            rg="RGC"
        fi
        
        if [ $rungroup == "rg-a" ] && [ "$ana" == "data" ]; then
            ana="nSidis"
        elif [ $rungroup == "rg-b" ] && [ "$ana" == "data" ]; then
            ana="sidisdvcs"
        elif [ $rungroup == "rg-c" ] && [ "$ana" == "data" ]; then
            ana="sidisdvcs"
        fi
        
        printgreen "\tRUNGROUP=$rungroup"
        if [ "$ana" == "nSidis" ] || [ "$ana" == "sidisdvcs" ]; then
            hipo_is_mc=0
        else
            hipo_is_mc=1
        fi
        hipodirs=$(get_hipo_dirs $ana "$rungroup")
        for hipodir in ${hipodirs[@]}
        do
            i=0
            for hipo in "$hipodir"*.hipo
            do
                if [ $i -eq $nFiles ]; then
                    break
                else
                    i=$((i+1))
                fi

                read runNumber beamE <<< $(hipo_beamE_runNumber $hipo $ana $rungroup)
                echo $hipo $runNumber
                slurmshell=$FARMOUT_DIR/shell/hipo2tree_${ana}_${runNumber}.sh
                slurmslurm=$FARMOUT_DIR/slurm/hipo2tree_${ana}_${runNumber}.slurm

                touch -f $slurmshell
                touch -f $slurmslurm

                chmod +x $slurmshell
                cat >> $slurmslurm << EOF
#!/bin/bash
#SBATCH --account=clas12
#SBATCH --partition=production
#SBATCH --mem-per-cpu=2000
#SBATCH --job-name=job_hipo2tree_$ana_$runNumber
#SBATCH --cpus-per-task=2
#SBATCH --time=24:00:00
#SBATCH --output=$FARMOUT_DIR/log/hipo2tree_$ana_$runNumber.out
#SBATCH --error=$FARMOUT_DIR/err/hipo2tree_$ana_$runNumber.err
$slurmshell
EOF

                echo "#!/bin/tcsh" >> $slurmshell
                echo "source /group/clas12/packages/setup.csh" >> $slurmshell
                echo "module load clas12/pro" >> $slurmshell

                # For loop over each dihadron pair
                j=0
                for pair in "${pion_pairs[@]}"; do
                    outfile="$VOLATILE_DIR/data/$pair/${ana}_${rg}_${runNumber}.root"
                    pid1=${pion_pairs_pids[${j},0]}
                    pid2=${pion_pairs_pids[${j},1]}
                    echo "clas12root -b -q $PWD/macros/hipo2tree.C\(\\\"${hipo}\\\",\\\"${outfile}\\\",$beamE,$pid1,$pid2,$nEvents,$hipo_is_mc\)" >> $slurmshell
                    j=$((j+1))
                done
                
                echo "Submitting slurm job for $hipo"
                sbatch --quiet $slurmslurm
            done
        done
    done
done

    
