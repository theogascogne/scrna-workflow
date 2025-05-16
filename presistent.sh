#!/bin/bash

# =============================
#
# Prototype to bypass Snakemake's default usage of 'Rscript' for each job's shell/run.
#
# Purpose:
# This experiment aims to forward Rscript commands from Snakemake to a persistent R session,
# with the goal of improving runtime by avoiding repeated Rscript startup overhead.
# This prototype has since been deprecated as it is unstable, may cause crashe due to maxing out RAM or have issues communicating with global.r
#
# Functionalities:
# - Error handling: acts as an intermediary to prompt retry, skip, or quit on command failure.
#   (This was an initial idea and may not be fully needed now.)
# - Logging: tracks CPU and RAM usage during command execution, but this feature is
#   currently not required anymore.
#
# How the prototype works:
# - Takes Rscript commands forwarded from Snakemake.
# - Uses three configurable variables (errorhandle, logging, background) to alter behavior.
#   If all are FALSE, commands run normally.
#
# For the 'background' mode (prototype feature):
# - If the trigger file ('trigger.txt') does not exist, launches 'global.r' and waits 10+ seconds for it to initialize.
# - Writes commands as key-value pairs into the trigger file, which global.r monitors.
# - Waits for global.r to update the trigger file with SUCCESS or FAILURE messages.
# - Relays this result back to Snakemake, allowing it to queue the next job accordingly.
# - If the trigger file already exists from a previous session, the script writes immediately to it.
#   (This can cause problems if the file persists across sessions; remember to delete 'trigger.txt'
#    after crash, graceful exit, or completion.)
# - The script's waiting behavior also forces Snakemake to wait for job completion.
#
# =============================


# Initialization of variables and constants
LOGFILE="workflow_metrics.log"
TRIGGER_FILE="$HOME/miniconda3/envs/test/lib/python3.9/site-packages/cellsnake/scrna/workflow/trigger.txt"

# Workflow-wide metrics initialization
WORKFLOW_START=$(date +%s)
TOTAL_CPU_SUM=0.0
TOTAL_CPU_SAMPLES=0
WORKFLOW_MAX_RAM_KB=0

# Flags for error handling, logging, and background execution
errorhandle=false  # Set true to enable retry on failure
logging=false      # Set true to enable resource usage logging
background=true    # Run global.r in background if true

# =============================
# Function: handle_failure
# 
# Prompt user for action on command failure:
# Retry (y), Skip (n), or Quit (q).
# Auto quits after 60 seconds of no input.
# =============================

handle_failure() {
  while true; do
    echo "Command failed. Retry? (y/n) or quit (q): (auto quit in 60s)"
    read -t 60 ans
    if [ $? -ne 0 ]; then
      echo "No input received in 60 seconds. Quitting the workflow."
      exit 1
    fi
    case "$ans" in
      y|Y) echo "Retrying..."; return 1 ;;
      n|N) echo "Skipping this job."; return 0 ;;
      q|Q) echo "Quitting the workflow."; exit 1 ;;
      *) echo "Invalid input. Please enter y, n, or q." ;;
    esac
  done
}

# =============================
# Function: run_with_metrics
#
# Runs a command in the background, tracks its CPU and RAM usage over time.
# Logs:
# - Command run
# - Duration (seconds)
# - Peak RAM usage (MB)
# - Average CPU usage (%)
# Updates workflow-wide aggregate metrics.
# =============================
run_with_metrics() {
  echo "Running: $*" | tee -a "$LOGFILE"
  start_time=$(date +%s)
  start_hr=$(date "+%Y-%m-%d %H:%M:%S")

  "$@" &
  CMD_PID=$!

  CPU_SAMPLES=()
  RAM_SAMPLES=()

  # Poll CPU and RAM usage every second while command runs
  while kill -0 $CMD_PID 2>/dev/null; do
    CPU=$(ps -p $CMD_PID -o %cpu= | awk '{print $1}')
    RAM_KB=$(ps -p $CMD_PID -o rss= | awk '{print $1}')

    if [[ "$CPU" =~ ^[0-9.]+$ ]]; then
      CPU_SAMPLES+=("$CPU")
    fi
    if [[ "$RAM_KB" =~ ^[0-9]+$ ]]; then
      RAM_SAMPLES+=("$RAM_KB")
    fi

    sleep 1
  done

  wait $CMD_PID
  STATUS=$?

  end_time=$(date +%s)
  duration=$((end_time - start_time))

  max_ram_kb=0
  cpu_sum=0.0

  # Calculate peak RAM usage
  for ram in "${RAM_SAMPLES[@]}"; do
    if [ "$ram" -gt "$max_ram_kb" ]; then max_ram_kb=$ram; fi
  done

  # Calculate total CPU usage sum
  for cpu in "${CPU_SAMPLES[@]}"; do
    cpu_sum=$(awk -v a="$cpu_sum" -v b="$cpu" 'BEGIN {print a + b}')
  done

  sample_count=${#CPU_SAMPLES[@]}
  cpu_avg=$(awk -v sum="$cpu_sum" -v count="$sample_count" 'BEGIN {if (count > 0) printf "%.2f", sum / count; else print 0}')
  ram_mb=$(awk -v kb="$max_ram_kb" 'BEGIN {printf "%.2f", kb / 1024}')

  # Update workflow-wide aggregated metrics
  TOTAL_CPU_SUM=$(awk -v a="$TOTAL_CPU_SUM" -v b="$cpu_sum" 'BEGIN {print a + b}')
  TOTAL_CPU_SAMPLES=$((TOTAL_CPU_SAMPLES + sample_count))
  if [ "$max_ram_kb" -gt "$WORKFLOW_MAX_RAM_KB" ]; then
    WORKFLOW_MAX_RAM_KB=$max_ram_kb
  fi

  # Log metrics to file
  {
    echo "---------------------------------------------"
    echo "Command         : $*"
    echo "Duration        : $duration seconds"
    echo ""
    echo "Peak RAM Usage  : ${ram_mb} MB"
    echo "Avg CPU Usage   : ${cpu_avg} %"
    echo "---------------------------------------------"
    echo ""
  } >> "$LOGFILE"

  return $STATUS
}

# Writes the command to a trigger file that global.r watches to start processing.
write_to_trigger_file() {
  echo "Writing to trigger file..."
  echo "--rscript $@" > "$TRIGGER_FILE"
}


# Waits until the trigger file exists (created/updated by global.r).
wait_for_trigger() {
  echo "Waiting for trigger file to be created/updated by global.r..."
  while [ ! -f "$TRIGGER_FILE" ]; do
    sleep 1
  done
}


# Monitors the trigger file for completion status:
# Detects "success" or "failure" strings to determine run result.
monitor_trigger_file() {
  echo "Waiting for global.r to finish..."
  while true; do
    if grep -q "success" "$TRIGGER_FILE"; then
      RESULT="success"
      break
    elif grep -q "failure" "$TRIGGER_FILE"; then
      RESULT="failure"
      break
    fi
    sleep 1
  done
}

# =============================
# Main execution block
# =============================
if [ "$background" = true ]; then
  # If trigger file doesn't exist, start global.r in background
  if [ ! -f "$TRIGGER_FILE" ]; then
    echo "Trigger file not found. Launching global.r..."
    Rscript ~/miniconda3/envs/test/lib/python3.9/site-packages/cellsnake/scrna/global.r &
    sleep 10  # Give global.r time to initialize
  fi

  # Write command to trigger file
  write_to_trigger_file "$@"

  # Wait for global.r to finish, monitor status
  monitor_trigger_file

  # Handle success or failure
  if [[ "$RESULT" == "success" ]]; then
    echo "Global.R completed successfully."
    echo "" > "$TRIGGER_FILE"
    exit 0
  else
    echo "Global.R failed."
    echo "" > "$TRIGGER_FILE"
    exit 1
  fi
else
  # If not background mode, run command normally
  if [ "$logging" = true ]; then
    run_with_metrics "$@"
    STATUS=$?
  else
    "$@"
    STATUS=$?
  fi

  # If error handling enabled and command failed, trigger failure handler
  if [ "$errorhandle" = true ] && [ $STATUS -ne 0 ]; then
    handle_failure
  fi
fi


