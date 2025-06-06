#!/usr/bin/env Rscript

# Start with a message that the script is starting
cat("Starting global.r script...\n")

# Define the trigger file path immediately
trigger_file <- path.expand("~/miniconda3/envs/test/lib/python3.9/site-packages/cellsnake/scrna/workflow/trigger.txt")

# Ensure the trigger file exists at startup
if (!file.exists(trigger_file)) {
  # Create the trigger file if it doesn't exist
  file.create(trigger_file)
  cat("Trigger.txt was created at the start.\n")
}

# Load required libraries
cat("Loading libraries...\n")
tryCatch({
  require(optparse)
  require(tidyverse)
  require(Seurat)
  require(patchwork)
  require(tools)
  require(data.table)
}, error = function(e) {
  cat("Error loading libraries: ", e$message, "\n")
  write("failure", file = trigger_file)
  stop("Libraries failed to load. Exiting...\n")
})

cat("Libraries loaded successfully.\n")

cat("Monitoring trigger file for commands...\n")

# Function to parse the trigger file and create an option list dynamically
parse_trigger_file <- function(trigger_file) {
  trigger_content <- readLines(trigger_file, warn = FALSE)
  
  if (length(trigger_content) > 0) {
    parts <- unlist(strsplit(trigger_content, " "))  # Split by space
    parts <- parts[parts != ""]  # Remove empty parts
    
    option_list <- list()
    
    i <- 1
    while (i <= length(parts)) {
      var_name <- parts[i]
      
      # If it's a standalone command (like exit_script), store it
      if (grepl("^--", var_name)) {
        if (i + 1 <= length(parts) && !grepl("^--", parts[i + 1])) {
          # If next part isn't a flag (i.e., it’s an argument), assign it
          var_value <- parts[i + 1]
          option_list[[var_name]] <- var_value
          i <- i + 2  # Move to the next argument pair
        } else {
          # If it's a flag without a value, set to TRUE (this is where flags like --doublet.filter will be set as TRUE)
          option_list[[var_name]] <- TRUE
          i <- i + 1
        }
      } else {
        option_list[[var_name]] <- TRUE  # Handle exit_script as a standalone command
        i <- i + 1
      }
    }
    
    return(option_list)
  }
  return(NULL)
}

run_r_script <- function(option_list) {
  if (!"--rscript" %in% names(option_list)) {
    cat("Error: No --rscript specified in trigger file.\n")
    write("failure", file = trigger_file)
    return()
  }
  
  # Clear trigger file *before* executing (prevents repeating the command)
  write("", file = trigger_file)
  
  r_script <- option_list[["--rscript"]]
  option_list[["--rscript"]] <- NULL
  
  # Construct the argument list to pass to the R script
  args_vector <- unlist(lapply(names(option_list), function(opt_name) {
    c(opt_name, option_list[[opt_name]])
  }))
  
  cat("Running", r_script, "with args:", paste(args_vector, collapse = " "), "\n")
  
  script_env <- new.env()
  script_env$argv <- args_vector
  script_env$parse_args <- function(parser) {
    parse_args(parser, args = args_vector)
  }
  
  tryCatch({
    # Start execution of the R script
    source(r_script, local = script_env)
    cat("Execution completed successfully.\n")
    
    # Write "success" to the trigger file
    write("success", file = trigger_file)
  }, error = function(e) {
    cat("Error during execution:", e$message, "\n")
    
    # Write "failure" to the trigger file if something goes wrong
    write("failure", file = trigger_file)
  })
}

exit_script <- function() {
  cat("Exit signal received. Stopping global.r...\n")
  
  # Remove the trigger file before exiting
  if (file.exists(trigger_file)) {
    file.remove(trigger_file)  # Delete the trigger file upon exit
    cat("trigger.txt was deleted upon exit.\n")
  }
  
  # Exit the script
  q("no")
}

# Ensure trigger.txt is cleared on script exit
on.exit({
  if (file.exists(trigger_file)) {
    file.remove(trigger_file)  # Delete the trigger file upon script exit
    cat("trigger.txt was deleted on exit.\n")
  }
})

# Main loop that waits for the trigger file and executes commands
while (TRUE) {
  if (file.exists(trigger_file)) {
    trigger_content <- readLines(trigger_file, warn = FALSE)
    
    # Proceed only if trigger_content is non-empty and valid
    if (length(trigger_content) > 0 && nchar(trimws(trigger_content)) > 0) {
      
      # Parse the trigger content and create an option list
      option_list <- parse_trigger_file(trigger_file)
      
      if (!is.null(option_list)) {
        # Check if the option list contains the script to run
        if ("--rscript" %in% names(option_list)) {
          run_r_script(option_list)  # Run the script with the dynamic option list
        } else if ("exit_script" %in% names(option_list)) {
          exit_script()
        }
      }
    }
  }
  
  # Sleep for 1 second before checking the trigger file again
  Sys.sleep(1)
}
