#!/stornext/Home/data/allstaff/l/lmcintosh/mambaforge/envs/sage/bin/python

"""
This module is used to generate test files from existing Python scripts in a specified project directory and its subdirectories.
It provides functionality to remove test functions and their calls, and adds logging for function inputs and outputs.
The module scans through the specified project directory and its subdirectories, identifying Python scripts,
excluding those that already end with '_IO.py'.
For each valid script, it generates a new script with the same name but ending with '_IO.py'.
This new script includes the original code, excluding test functions, and with added input-output logging.

Usage: python script_name.py --directory /path/to/project_directory
"""

import os
import re
import argparse
from contextlib import contextmanager
from typing import Generator, Any
import inspect
import logging

logging.basicConfig(level=logging.INFO)

@contextmanager
def managed_file(file_path: str, mode: str) -> Generator:
    """
    A context manager for opening and closing files.
    :param file_path: Path to the file.
    :param mode: Mode to open the file in.
    :return: Yield the file object.
    """
    file = open(file_path, mode)
    try:
        yield file
    finally:
        file.close()

def strip_test_functions(code: str) -> str:
    """
    Strips test functions from the code.
    :param code: Input code from which to remove test functions.
    :return: Code with test functions removed.
    """
    code_lines = code.split('\n')

    inside_test_function = False
    stripped_lines = []

    for line in code_lines:
        if re.match(r'^def test_', line.strip()):
            inside_test_function = True
        elif re.search(r'\btest_\w+\(', line):  # Remove calls to test functions
            continue
        elif inside_test_function and line.strip() and not line.startswith(" "):
            inside_test_function = False
        if not inside_test_function:
            stripped_lines.append(line)

    return '\n'.join(stripped_lines)

def inject_io_logging(code: str) -> str:
    """
    Injects input-output logging into the code.
    :param code: Input code to which to add input-output logging.
    :return: Code with input-output logging added.
    """
    code_lines = code.split('\n')

    logging_decorator = [
        "from functools import wraps",
        "import os",
        "import inspect",
        "",
        "def log_io(func):",
        "    @wraps(func)",
        "    def wrapper(*args, **kwargs):",
        "        arg_names = inspect.getfullargspec(func).args",
        "        arg_values = {**dict(zip(arg_names, args)), **kwargs}",
        "        result = func(*args, **kwargs)",
        "",
        "        if not os.path.exists('IO_LOGGER'):",
        "            os.mkdir('IO_LOGGER')",
        "        with open(f'IO_LOGGER/{func.__name__}_io_log.txt', 'a') as f:",
        "            f.write(f'Input: {repr(arg_values)}\\n')",
        "            f.write(f'Output: {repr(result)}\\n')",
        "",
        "        return result",
        "    return wrapper",
        "",
    ]

    lines_with_logging = logging_decorator + [""]

    for line in code_lines:
        if line.strip().startswith("def "):
            indentation = line.index("def")
            lines_with_logging.append(" " * indentation + '@log_io')
        lines_with_logging.append(line)

    return '\n'.join(lines_with_logging)

def adjust_imports(code: str, base_dir: str) -> str:
    """
    Adjusts the import statements in the code based on the base directory.
    :param code: Input code in which to adjust imports.
    :param base_dir: Base directory against which to adjust imports.
    :return: Code with adjusted imports.
    """
    code_lines = code.split('\n')
    adjusted_lines = []

    for line in code_lines:
        from_import_match = re.match(r'^from (.*) import (.*)$', line)
        import_match = re.match(r'^import (.*)$', line)
        if from_import_match:
            module, _ = from_import_match.groups()
            if base_dir in module:
                adjusted_lines.append(line.replace(module, module + '_IO'))
            else:
                adjusted_lines.append(line)
        elif import_match:
            module = import_match.groups()[0]
            if base_dir in module:
                adjusted_lines.append(f'import {module}_IO')
            else:
                adjusted_lines.append(line)
        else:
            adjusted_lines.append(line)

    return '\n'.join(adjusted_lines)

def generate_io_scripts(directory_path: str) -> None:
    """
    Generates IO scripts from Python scripts in the specified directory and subdirectories.
    :param directory_path: Path to the base directory to scan for Python scripts.
    """
    for root, _, files in os.walk(directory_path):
        for file in files:
            if file.endswith('.py') and not file.endswith('_IO.py'):
                input_file = os.path.join(root, file)
                output_file = input_file.rsplit(".py", 1)[0] + "_IO.py"

                with managed_file(input_file, 'r') as f:
                    code = f.read()

                code = strip_test_functions(code)
                code = inject_io_logging(code)
                code = adjust_imports(code, directory_path)

                with managed_file(output_file, 'w') as f:
                    f.write(code)
                
                logging.info(f"Generated {output_file} from {input_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate IO scripts from Python scripts in the given directory.')
    parser.add_argument('--directory', required=True, help='Base directory to scan for Python scripts.')
    args = parser.parse_args()

    generate_io_scripts(args.directory)
