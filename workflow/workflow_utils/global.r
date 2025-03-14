#!/usr/bin/env Rscript

# Start with a message that the script is starting
cat("Starting global.r script...\n")

# Load required libraries
cat("Loading libraries...\n")
tryCatch({
  library(optparse)  # Ensure optparse is loaded
  library(Seurat)
}, error = function(e) {
  cat("Error loading libraries: ", e$message, "\n")
  stop("Libraries failed to load. Exiting...\n")
})

cat("Libraries loaded successfully.\n")

# Define the trigger file path
trigger_file <- path.expand("~/miniconda3/envs/test/lib/python3.9/site-packages/cellsnake/scrna/workflow/trigger.txt")
cat("Monitoring trigger file for commands...\n")

# Ensure the trigger file exists at startup
if (file.exists(trigger_file)) {
  # Clear the contents of the trigger file (without deleting it)
  write("", file = trigger_file)
  cat("Existing trigger.txt was cleared at the start.\n")
} else {
  # Create the trigger file if it doesn't exist
  file.create(trigger_file)
  cat("Trigger.txt was created at the start.\n")
}

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
        if (i + 1 <= length(parts)) {
          var_value <- parts[i + 1]
          option_list[[var_name]] <- var_value
          i <- i + 2  # Move to the next argument pair
        } else {
          option_list[[var_name]] <- TRUE  # Handle flags with no values
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
    return()
  }

  r_script <- option_list[["--rscript"]]  # Extract script path
  option_list[["--rscript"]] <- NULL  # Remove r_script path from options

  cat("Running", r_script, "with dynamic arguments...\n")

  # Construct optparse argument list
  args_vector <- unlist(lapply(names(option_list), function(opt_name) {
    c(opt_name, option_list[[opt_name]])  # Ensure this line has correct parentheses
  }))

  cat("Simulated command-line arguments:\n", paste(args_vector, collapse = " "), "\n")

  # Create an isolated environment for the script
  script_env <- new.env()

  # Pass arguments into the script environment
  script_env$argv <- args_vector  # Simulate command-line arguments

  # Override parse_args to force optparse to use our args
  script_env$parse_args <- function(parser) {
    parse_args(parser, args = args_vector)
  }

  # Debugging Output
  cat("Starting execution of", r_script, "at", format(Sys.time()), "\n")

  # Source the script within this environment
  source(r_script, local = script_env)

  cat(r_script, "execution completed at", format(Sys.time()), "\n")
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

          # Clear trigger file after processing
          write("", file = trigger_file)
        } else if ("exit_script" %in% names(option_list)) {
          exit_script()
        }
      }
    }
  }
  
  # Sleep for 1 second before checking again
  Sys.sleep(1)
}