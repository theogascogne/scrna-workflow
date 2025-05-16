# client.R

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript client.R <message>")
}

msg <- args[1]

host <- "127.0.0.1"
port <- 6011

con <- tryCatch(
  socketConnection(
    host = host,
    port = port,
    blocking = TRUE,
    open = "r+b"
  ),
  error = function(e) {
    stop("Connection failed: ", e$message)
  }
)

writeLines(msg, con)

response <- tryCatch(readLines(con, 1), error = function(e) "No response")
cat("Response from server:", response, "\n")

close(con)