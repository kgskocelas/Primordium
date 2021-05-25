#!/bin/bash

if [ "$#" != "1" ]; then
    echo "Error! Expecting one argument: the name of the experiment!"
    echo "Recommended format: YYYY_MM_DD_name_here"
else
    echo "Creating experiment: $1!"
    # Grab some useful variables
    # Experiment name -> command line argument
    SR_EXP_NAME=$1
    # Root directory -> The root level of the repo, should be directory just above 'experiments'
    SR_ROOT_DIR=$(pwd | grep -oP ".+/(?=experiments)")
    # Experiment directory -> current directory
    SR_EXP_DIR=${SR_ROOT_DIR}/experiments/${SR_EXP_NAME}

    mkdir ${SR_EXP_DIR}
    cp ${SR_ROOT_DIR}/experiments/scripts/templates/experiment/* ${SR_EXP_DIR}
    echo "Experiment created! It can be found at: ${SR_EXP_DIR}"
    echo "To get started, navigate to that dir and edit config.sh"
    echo "After that, start running the bash scripts, starting with 00_generate_timing_distributions.sh"
    echo "Don't forget to edit the README.md in the experiment directory!"
fi
