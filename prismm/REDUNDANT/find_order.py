import os
import sys
import re
from collections import defaultdict


def extract_io_files(file_path):
    with open(file_path, 'r') as file:
        content = file.read()

    input_files = re.findall(r'(?:open|with open)\s*\(\s*[\'\"](.*?)[\'\"]\s*,\s*[\'\"](?:r|rt)[\'\"]', content)
    output_files = re.findall(r'(?:open|with open)\s*\(\s*[\'\"](.*?)[\'\"]\s*,\s*[\'\"](?:w|wt|a|at|x|xt)[\'\"]', content)

    return input_files, output_files


def summarize_io_files(directory):
    io_files = defaultdict(lambda: {"input": set(), "output": set()})

    for root, _, files in os.walk(directory):
        for file in files:
            if not (file.endswith('.py') or file.endswith('.sh') or file.endswith('.sage')):
                continue

            file_path = os.path.join(root, file)
            input_files, output_files = extract_io_files(file_path)

            io_files[file_path]["input"] = set(input_files)
            io_files[file_path]["output"] = set(output_files)

    return io_files


def main(directory):
    io_files_summary = summarize_io_files(directory)

    print("Files read and written by each program:")
    for program, io_files in io_files_summary.items():
        print(f"\nProgram: {os.path.relpath(program, directory)}")
        print("  Input files:")
        for input_file in io_files["input"]:
            print(f"    - {input_file}")
        print("  Output files:")
        for output_file in io_files["output"]:
            print(f"    - {output_file}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python summarize_io_files.py <path_to_directory>")
        sys.exit(1)

    main(sys.argv[1])

