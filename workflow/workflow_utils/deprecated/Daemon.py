import subprocess
import os
import time
import sys
import daemon

# Paths
LOCK_FILE = "/tmp/global_r.lock"
LOG_FILE = os.path.expanduser('~/miniconda3/envs/test/lib/python3.9/site-packages/cellsnake/scrna/workflow/daemon.log')

# Function to log messages
def log_message(message):
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    with open(LOG_FILE, 'a') as log:  # Append mode
        log.write(f"{timestamp}: {message}\n")

# Function to start an R process and load Seurat
def start_r_process(script_path, log_output=True):
    if os.path.exists(LOCK_FILE):
        if log_output:
            log_message("global.r is already running. Exiting...")
        return  

    with open(LOCK_FILE, 'w') as lock:
        lock.write("locked")

    try:
        if log_output:
            log_message(f"Starting R process: {script_path}...")

        process = subprocess.Popen(
            ["Rscript", script_path], 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE
        )

        # Capture output and write to log (only for daemon mode)
        if log_output:
            for line in iter(process.stdout.readline, b""):
                log_message(f"stdout: {line.decode().strip()}")
            for line in iter(process.stderr.readline, b""):
                log_message(f"stderr: {line.decode().strip()}")

        process.stdout.close()
        process.stderr.close()
        process.wait()

    finally:
        os.remove(LOCK_FILE)
        if log_output:
            log_message("R process completed. Lock file removed.")

# Function to run as a daemon
def run_daemon(script_path, trigger_filename="trigger.txt"):
    open(LOG_FILE, 'w').close()  # Clear log at daemon start

    # Wait for the trigger file to be created
    trigger_file_path = os.path.join(os.path.dirname(script_path), trigger_filename)

    while not os.path.isfile(trigger_file_path):
        log_message(f"Waiting for {trigger_filename} to be created by global.r...")
        time.sleep(1)  # Check every 1 second

    log_message(f"{trigger_filename} detected. Proceeding with daemon...")

    with daemon.DaemonContext():
        log_message("Daemon running in background...")
        start_r_process(script_path, log_output=True)  # Log output in daemon mode

# Function to run in the background (&) without daemonizing
def run_in_background(script_path, trigger_filename="trigger.txt"):
    # Start R process
    process = subprocess.Popen(["Rscript", script_path])

    # Path to trigger file
    trigger_file_path = os.path.join(os.path.dirname(script_path), trigger_filename)

    try:
        while not os.path.isfile(trigger_file_path):
            print(f"Waiting for {trigger_filename} to be created by global.r...")
            time.sleep(10)  # Check every 10 seconds

        print(f"{trigger_filename} detected. Closing daemon.py.")

    finally:
        print("Daemon.py is exiting.")
        sys.exit(0)  # Gracefully exit the script



# Main function
def main():
    if len(sys.argv) < 3:
        print("Usage: python daemon.py <script_path> <mode (d/bg)>")
        sys.exit(1)

    script_path = sys.argv[1]
    mode = sys.argv[2].lower()  # Convert to lowercase

    if mode == "bg":
        run_in_background(script_path)
    elif mode == "d":
        run_daemon(script_path)
    else:
        print("Invalid mode. Use 'd' for daemon or 'bg' for background.")
        sys.exit(1)

if __name__ == "__main__":
    main()