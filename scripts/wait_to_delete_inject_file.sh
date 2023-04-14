#!/bin/bash

function wait_for_jobs() {
    local job_name=$1
    local jobsLeft=-999
    
    while [ $jobsLeft -ne 0 ]
    do
        read jobsLeft <<< $(echo "$(squeue -u gmat --format="%.18i %.9P %.30j %.8u %.8T %.10M %.9l %.6D %R")" | grep $job_name | awk 'END{print NR}')
        echo "Jobs remaining: " $jobsLeft
        sleep 30
    done
}

wait_for_jobs $1