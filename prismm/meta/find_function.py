# Import necessary modules and functions
import os
import argparse
from typing import Dict, List
from directory_summarizer import summarize_directories

def find_function_locations(summary: Dict[str, Dict[str, Dict[str, List[str]]]], function_name: str) -> Dict[str, List[str]]:
    locations = {}
    for dir_name, files in summary.items():
        for file_name, content in files.items():
            if function_name in content['functions']:
                if dir_name in locations:
                    locations[dir_name].append(file_name)
                else:
                    locations[dir_name] = [file_name]
    return locations

def main() -> None:
    # Initialize argument parser
    parser = argparse.ArgumentParser(description='Find specific function in Python files.')
    parser.add_argument('dirpaths', nargs='*', default=[os.getcwd()], help='The directories to search')
    parser.add_argument('function_name', help='The function name to find')
    args = parser.parse_args()

    # Use function from directory_summarizer to summarize directories
    summary = summarize_directories(args.dirpaths)
    # Find function locations
    locations = find_function_locations(summary, args.function_name)

    if not locations:
        print(f"The function '{args.function_name}' is not found.")
    else:
        print(f"The function '{args.function_name}' is found in the following locations:")
        for dir_name, files in locations.items():
            print(f"Directory: {dir_name}")
            for file_name in files:
                print(f"    File: {file_name}")

if __name__ == "__main__":
    main()
