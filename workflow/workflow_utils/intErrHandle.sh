#!/bin/bash

# Function to handle retry/skip/quit behavior
handle_failure() {
  while true; do
    echo "Command failed. Retry? (y/n) or quit (q): "
    read ans
    if [[ "$ans" == "y" ]]; then
      echo "Retrying..."
      return 1  # Retry the command by returning non-zero status
    elif [[ "$ans" == "q" ]]; then
      echo "Quitting the workflow."
      exit 1  # Exit the entire workflow
    else
      echo "Invalid input. Please enter y, n, or q."
    fi
  done
}

# Run the command and keep it alive to handle failure
while true; do
  "$@" || handle_failure  # Run the command and call handle_failure if it fails
  if [ $? -eq 0 ]; then    # If command is successful, break out of the loop
    break
  fi
done