# src/utils/utils.py

import subprocess
import shutil
from pathlib import Path
import logging
import sys
from datetime import datetime

def setup_logging(log_dir: Path, log_name_prefix: str):
    """
    Sets up logging to write to both a file and the console.

    Args:
        log_dir (Path): The directory where log files will be stored.
        log_name_prefix (str): A descriptive prefix for the log file name,
                               e.g., "yellow-fever_pipeline" or "rsv-a_fetch_data".
    """
    # Ensure the logs directory exists
    log_dir.mkdir(parents=True, exist_ok=True)
    
    # Create a timestamped and descriptive log file name
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filepath = log_dir / f"{timestamp}_{log_name_prefix}.log"

    # Get the root logger
    logger = logging.getLogger()
    logger.setLevel(logging.INFO) # Set the lowest level to capture

    # Prevent adding handlers multiple times if this function is called again
    if logger.hasHandlers():
        logger.handlers.clear()

    # Create file handler for detailed logging
    file_handler = logging.FileHandler(log_filepath)
    file_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(file_formatter)
    logger.addHandler(file_handler)

    # Create console handler for high-level status updates
    console_handler = logging.StreamHandler(sys.stdout)
    console_formatter = logging.Formatter('%(levelname)s: %(message)s')
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)
    
    logging.info(f"Logging initialized. Log file at: {log_filepath}")
    return log_filepath

def run_command(command: str, working_dir: Path):
    """Executes a shell command in a specified directory and checks for errors."""
    logging.info(f"Executing: {command}")
    try:
        result = subprocess.run(
            command,
            shell=True,
            check=True,
            capture_output=True,
            text=True,
            cwd=working_dir
        )
        if result.stdout:
            logging.debug(f"Stdout: {result.stdout}")
    except subprocess.CalledProcessError as e:
        # Log the error details before raising it
        logging.error(f"Command failed with exit code {e.returncode}")
        logging.error(f"Command: {e.cmd}")
        logging.error(f"Stderr: {e.stderr}")
        logging.error(f"Stdout: {e.stdout}")
        raise # Re-raise the exception to be handled by the calling function

def find_and_rename(start_path: Path, pattern: str, new_name: str, destination_dir: Path):
    """
    Finds a single file matching a pattern, renames it, and moves it.
    Logs a warning if multiple files match the pattern but still proceeds with the first file.
    """
    found_files = list(start_path.rglob(pattern))
    
    if len(found_files) > 1:
        logging.warning(f"Found multiple files ({len(found_files)}) matching pattern '{pattern}' in {start_path}.")
        logging.warning(f"Proceeding with the first file found: {found_files[0]}")

    if found_files:
        source_file = found_files[0]
        destination_file = destination_dir / new_name
        logging.info(f"Found: {source_file}. Moving and renaming to {destination_file}")
        shutil.move(source_file, destination_file)
    else:
        logging.warning(f"Could not find any file matching pattern '{pattern}' in {start_path}")

def find_and_move_all(start_path: Path, pattern: str, destination_dir: Path):
    """
    Finds all files matching a pattern and moves them to a new directory,
    preserving their original names.
    """
    found_files = list(start_path.rglob(pattern))
    if found_files:
        logging.info(f"Found {len(found_files)} files matching '{pattern}'. Moving all to {destination_dir}...")
        for source_file in found_files:
            try:
                # shutil.move can move a file to a directory
                shutil.move(source_file, destination_dir)
            except shutil.Error as e:
                # This can happen if a file with the same name already exists at the destination
                logging.error(f"Could not move file {source_file} to {destination_dir}. A file with the same name might already exist. Error: {e}")
    else:
        logging.warning(f"Could not find any files matching pattern '{pattern}' in {start_path}")

