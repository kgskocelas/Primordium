#!/bin/bash

# Load in config options
source ./config.sh

# Handle command line args
DO_RUN=0
while getopts ":r" opt; do
  case $opt in
    r)
     DO_RUN=1 
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      ;;
  esac
done

# Create necessary directories and clear old data if present
mkdir -p ${SR_TIMING_OUTPUT_DIR}
mkdir -p ${SR_TIMING_JOB_DIR}
rm ${SR_TIMING_JOB_DIR}/* -r > /dev/null 2> /dev/null

# Prep variables for local format (in this case, flags)
SR_ONE_CHECK_LOCAL=@@one_check
if [ "${SR_ONE_CHECK,,}" == "false" -o "${SR_ONE_CHECK,,}" == "f" -o "${SR_ONE_CHECK,,}" == "0" ]; 
    then
    SR_ONE_CHECK_LOCAL=@@multi_check
fi
SR_INFINITE_LOCAL=@@infinite
if [ "${SR_INFINITE,,}" == "false" -o "${SR_INFINITE,,}" == "f" -o "${SR_INFINITE,,}" == "0" ]; 
    then
    SR_INFINITE_LOCAL=@@finite
fi

# Generate the jobs!
python3 ${SR_ROOT_DIR}/experiments/scripts/timing_job_prep.py @@executable_path ${SR_ROOT_DIR}/application/bin/SpatialRestraint @@output_dir ${SR_TIMING_OUTPUT_DIR} @@job_dir ${SR_TIMING_JOB_DIR} @@ones ${SR_TIMING_ONES} @@cost ${SR_COST} @@mc_size ${SR_MC_SIZE} @@mut_rate ${SR_TIMING_MUT_RATE} @@threshold ${SR_THRESHOLD} @@samples ${SR_SAMPLES} @@num_tasks ${SR_TIMING_TASKS} @@seed_offset ${SR_TIMING_SEED_OFFSET} @@time ${SR_TIMING_TIME} @@memory ${SR_TIMING_MEMORY} ${SR_ONE_CHECK_LOCAL} ${SR_INFINITE_LOCAL}

# Create instance of roll_q for timing jobs
cp ${SR_ROOT_DIR}/experiments/scripts/third_party/roll_q ${SR_TIMING_DIR}/roll_q -r
# Setup roll_q files for timing jobs
find ${SR_TIMING_JOB_DIR} -name *.sb | sort > ${SR_TIMING_DIR}/roll_q/roll_q_job_array.txt
echo "0" > ${SR_TIMING_DIR}/roll_q/roll_q_idx.txt

# Launch timing jobs via roll_q!
if [ "${DO_RUN}" == "1" ]; then
    echo "Launching jobs!"
    ${SR_TIMING_DIR}/roll_q/roll_q.sh
fi
if [ "${DO_RUN}" == "0" ]; then
    echo "Run with -r to launch jobs!"
fi
  


# Shift-J in Vim to gather selected lines into one line
#python3 ${SR_ROOT_DIR}/experiments/scripts/timing_job_prep.py
#    --executable_path ${ROOT_DIR}/application/bin/SpatialRestraint 
#    --output_dir ${SR_TIMING_OUTPUT_DIR}
#    --job_dir ${SR_TIMING_JOB_DIR}
#    --ones ${SR_TIMING_ONES}
#    --cost ${SR_COST}
#    --mc_size ${SR_MC_SIZE}
#    --mut_rate ${SR_TIMING_MUT_RATE}
#    --threshold ${SR_THRESHOLD}
#    --samples ${SR_SAMPLES}
#    --num_tasks ${SR_TIMING_TASKS}
#    --seed_offset ${SR_TIMING_SEED_OFFSET}
#    --time ${SR_TIMING_TIME}
#    --memory ${SR_TIMING_MEMORY}
#    ${SR_ONE_CHECK_LOCAL}
#    ${SR_INFINITE_LOCAL}

