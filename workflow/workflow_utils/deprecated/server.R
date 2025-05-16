#!/usr/bin/env Rscript

require(optparse)
require(tidyverse)
require(Seurat)
require(patchwork)
require(tools)
require(data.table)

host <- "127.0.0.1"
port <- 6011

cat("Server starting on", paste0(host, ":", port), "\n")

repeat {
  cat("Waiting for client...\n")
  con <- tryCatch(
    socketConnection(
      host = host,
      port = port,
      server = TRUE,
      blocking = TRUE,
      open = "a+b"
    ),
    error = function(e) {
      cat("Connection error:", e$message, "\n")
      return(NULL)
    }
  )
  
  if (is.null(con)) {
    Sys.sleep(1)
    next
  }
  
  cat("Client connected.\n")
  
  repeat {
    msg <- tryCatch(readLines(con, 1), error = function(e) "")
    
    if (length(msg) == 0 || msg == "") {
      cat("Client disconnected.\n")
      break
    }
    
    cat("Received:", msg, "\n")
    response <- ""
    
    if (msg == "quit") {
      response <- "Shutting down server."
      writeLines(response, con)
      close(con)
      quit(save = "no")
    } else {
      response <- tryCatch({
        # Split the message into parts (script path + args)
        parts <- strsplit(msg, " +")[[1]]
        script_path <- parts[1]
        args_vector <- parts[-1]
        
        if (!file.exists(script_path)) {
          stop("Script not found: ", script_path)
        }
        
        cat("Running R script with arguments:\n")
        cat("Script:", script_path, "\n")
        cat("Args:", paste(args_vector, collapse = " "), "\n")
        
        # Create an isolated environment for the script
        script_env <- new.env()
        
        # Simulate command-line arguments
        script_env$argv <- args_vector
        
        # Override parse_args in the environment
        script_env$parse_args <- function(parser) {
          optparse::parse_args(parser, args = args_vector)
        }
        
        # Source the script in the prepared environment
        source(script_path, local = script_env)
        
        "Script executed successfully."
      }, error = function(e) {
        err_msg <- paste("ERROR executing command:", e$message)
        cat(err_msg, "\n")
        err_msg
      })
    }
    
    writeLines(response, con)
  }
  
  close(con)
}