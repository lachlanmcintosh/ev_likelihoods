#!/usr/bin/env python3

"""
Module: Git Operations

This module automates certain git operations within the current directory and its subdirectories. 
The operations include staging deleted files (excluding those ending with '_IO.py') and adding new or modified files with certain extensions (excluding those ending with '_IO.py'). 
"""

import os
import subprocess
from typing import List, Tuple

def stage_deleted_files() -> None:
    """
    Stage deleted files for commit in the current git repository,
    excluding files ending with "_IO.py".
    """
    result = subprocess.run(["git", "ls-files", "--deleted"], capture_output=True, text=True)
    deleted_files = result.stdout.splitlines()

    for file in deleted_files:
        if not file.endswith('_IO.py'):  # Stage only if not ending in '_IO.py'
            subprocess.run(["git", "rm", "--cached", file])

def add_files_with_extensions(dirpath: str, filenames: List[str], extensions: Tuple[str, ...]) -> None:
    """
    Add files with specified extensions in the current git repository,
    excluding files ending with "_IO.py".

    Args:
        dirpath (str): The path of the directory to search files in.
        filenames (List[str]): The names of the files in the directory.
        extensions (Tuple[str, ...]): The file extensions to search for.
    """
    assert os.path.isdir(dirpath), f"Invalid directory: {dirpath}"
    
    for file in filenames:
        if not file.endswith('_IO.py') and os.path.splitext(file)[1] in extensions:
            filepath = os.path.join(dirpath, file)
            assert os.path.isfile(filepath), f"Invalid file: {filepath}"
            subprocess.run(['git', 'add', filepath])

def git_add_files() -> None:
    """
    Perform git add operation on files in the current directory and its subdirectories.
    Stages deleted files (excluding those ending with "_IO.py") and adds new or modified files with specified extensions (excluding those ending with "_IO.py").
    """
    stage_deleted_files()

    extensions = ('.sh', '.py', '.sage', '.R', '.Rnw', '.Rmd')

    # Traverse through all subdirectories and find files with the specified extensions
    for dirpath, _, filenames in os.walk('.'):
        add_files_with_extensions(dirpath, filenames, extensions)

def main() -> None:
    """
    Main function to execute the git add operations.
    """
    git_add_files()

    # The "TODO" file is added here, consider refactoring if used frequently
    subprocess.run(["git", "add", "TODO"])

if __name__ == '__main__':
    main()
