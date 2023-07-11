"""
Split File Script

This Python script splits a given file into smaller parts based on given line numbers. 
The file to be split and the line numbers are passed as command line arguments. 

Additionally, an output directory is specified, where the resulting files are saved. 
If the directory does not exist, the script will create it. 

The resulting files are named as 'run_simulation_part_X.py', where X is the part number. 

Usage: 
python3 split_file.py <file_path> <breakpoints> <output_dir>

- file_path: Path to the file to be split.
- breakpoints: Line numbers where the file should be split. Multiple line numbers can be given separated by space.
- output_dir: Directory to save the output files.
"""

import linecache
import os
import argparse
import logging
from typing import List

# Setup logging
logging.basicConfig(level=logging.INFO)

# Define the argument parser within a function
def setup_arg_parser() -> argparse.ArgumentParser:
    """Creates and returns an argument parser for command-line arguments."""
    parser = argparse.ArgumentParser(description='Break a file into smaller files.')
    parser.add_argument('file_path', type=str, help='Path to the file to be split.')
    parser.add_argument('breakpoints', type=int, nargs='+', help='Line numbers where the file should be split.')
    parser.add_argument('output_dir', type=str, help='Directory to save the output files.')
    return parser

def prepare_output_dir(directory: str) -> None:
    """Checks if the output directory exists, and if not, creates it."""
    if not os.path.exists(directory):
        os.makedirs(directory)

def split_file_into_segments(file_path: str, breaks: List[int], output_dir: str) -> None:
    """Splits a file into segments based on given line numbers and writes them to the output directory."""
    for i in range(len(breaks) - 1):
        start = breaks[i] + 1  # +1 because linecache counts from line 1
        end = breaks[i + 1]
        
        # Define the output file name based on the lines it contains
        output_file = os.path.join(output_dir, f"run_simulation_part_{i+1}.py")

        with open(output_file, "w") as f:
            line = start
            while True:
                if end is not None and line > end:  # If end is defined and line has reached end, break
                    break
                line_content = linecache.getline(file_path, line)
                if line_content == '':  # If line_content is empty, break, this means we've reached end of file
                    break
                f.write(line_content)
                line += 1
        logging.info(f"Segment {i+1} written to {output_file}")

def main() -> None:
    """Main function to execute the script."""
    parser = setup_arg_parser()
    args = parser.parse_args()

    file_path = args.file_path
    breaks = [0] + args.breakpoints
    breaks.append(None)  # None will allow us to read till the end of the file for the last segment
    output_dir = args.output_dir

    prepare_output_dir(output_dir)
    split_file_into_segments(file_path, breaks, output_dir)

if __name__ == "__main__":
    main()
