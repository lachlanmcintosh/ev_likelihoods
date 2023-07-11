"""
Script Description:

This script is designed to scan a directory and its subdirectories for Python files with an '_IO.py' extension. 
The script logs its findings and prompts the user to decide whether to delete these '_IO.py' files.

Usage: python <script_name> --directory <directory_path> [--confirm] [--delete_all]
"""

import os
import sys
import argparse
import logging
from typing import List, Optional

# Set up logging
logging.basicConfig(level=logging.INFO)

def find_io_files(directory_path: str, delete_all: bool) -> Optional[List[str]]:
    """Find _IO.py files in given directory and its subdirectories.

    Args:
        directory_path (str): Path of the directory to search in.
        delete_all (bool): Whether to delete all _IO.py files or just those without corresponding .py file.

    Returns:
        Optional[List[str]]: List of _IO.py files to delete.
    """
    assert os.path.isdir(directory_path), f"Provided path {directory_path} is not a directory"

    files_to_delete = []

    for root, _, files in os.walk(directory_path):
        for file in files:
            if file.endswith('_IO.py'):
                original_file = file[:-6] + ".py"
                original_file_path = os.path.join(os.path.abspath(root), original_file)
                if delete_all or not os.path.isfile(original_file_path):
                    logging.info(f"Marked for deletion: {os.path.join(os.path.abspath(root), file)}")
                    files_to_delete.append(os.path.join(os.path.abspath(root), file))
    return files_to_delete if files_to_delete else None

def delete_files(files: List[str]) -> None:
    """Delete provided list of files.

    Args:
        files (List[str]): List of files to delete.
    """
    for file in files:
        try:
            os.remove(file)
        except Exception as e:
            logging.error(f"Failed to delete {file}. Error: {e}")
        else:
            logging.info(f"Successfully deleted {file}")

def prompt_file_deletion(files_to_delete: List[str]) -> None:
    """Prompt the user for file deletion and delete files if user agrees.

    Args:
        files_to_delete (List[str]): List of files to delete.
    """
    logging.info("\nFound the following _IO.py files:")
    for file in files_to_delete:
        logging.info(file)
    response = input("\nDo you want to delete these files? (yes/no): ")
    if response.lower() == "yes":
        delete_files(files_to_delete)
    else:
        logging.info("\nFiles not deleted.")

def main(args):
    directory_path = args.directory
    delete_all = args.delete_all
    files_to_delete = find_io_files(directory_path, delete_all)

    if files_to_delete:
        prompt_file_deletion(files_to_delete)
    else:
        logging.info("No _IO.py files were found.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--directory', required=True, help='The directory to scan for _IO.py files.')
    parser.add_argument('--confirm', action='store_true', help='Confirm the base directory before proceeding.')
    parser.add_argument('--delete_all', action='store_true', help='Delete all _IO.py files, even if a corresponding .py file exists.')

    args = parser.parse_args()

    if args.confirm:
        print(f"The base directory is set as: {os.path.abspath(args.directory)}")
        proceed = input("Do you wish to proceed? (yes/no): ")
        if proceed.lower() != 'yes':
            print("Aborted.")
            sys.exit(0)
            
    main(args)
