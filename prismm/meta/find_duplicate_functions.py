# Built-in library imports
import ast
import sys
import logging
from collections import defaultdict
from typing import List, Dict, Union

# Python logging setup
logging.basicConfig(level=logging.INFO)


def find_duplicate_function_definitions(file_path: str) -> Dict[str, int]:
    """
    Given a Python source file, parse the source code and return any duplicate function definitions.

    Args:
        file_path (str): The path to the Python file to check for duplicate function definitions.

    Returns:
        duplicate_functions (Dict[str, int]): A dictionary with function names as keys and count as values.
    """

    # Validate input types
    assert isinstance(file_path, str), "file_path must be a string"

    # Read the source code from the provided file path
    with open(file_path, "r") as file:
        source_code = file.read()

    # Parse the source code into an AST
    parsed_ast = ast.parse(source_code)

    # Use a defaultdict to track function definition counts
    function_definitions = defaultdict(int)

    # Iterate through the parsed AST and increment the count for each function definition
    for node in ast.walk(parsed_ast):
        if isinstance(node, ast.FunctionDef):
            function_definitions[node.name] += 1

    # Filter function_definitions to only include functions defined more than once
    duplicate_functions = {name: count for name, count in function_definitions.items() if count > 1}

    return duplicate_functions


def report_duplicate_functions(duplicate_functions: Dict[str, int], file_path: str) -> None:
    """
    Report any duplicate function definitions.

    Args:
        duplicate_functions (Dict[str, int]): The dictionary of function names and their counts.
        file_path (str): The path to the Python file checked for duplicate function definitions.
    """
    
    # Validate input types
    assert isinstance(duplicate_functions, dict), "duplicate_functions must be a dictionary"
    assert isinstance(file_path, str), "file_path must be a string"

    if duplicate_functions:
        logging.info(f"Duplicate function definitions found in {file_path}:")
        for function_name, count in duplicate_functions.items():
            logging.info(f"  - {function_name} (defined {count} times)")
    else:
        logging.info(f"No duplicate function definitions found in {file_path}")


def main() -> None:
    """
    Main function to handle script execution.
    """
    # Validate command line arguments
    if len(sys.argv) < 2:
        logging.error("Usage: python find_duplicate_functions.py <path_to_script>")
        sys.exit(1)
    input_script_path = sys.argv[1]
    
    # Find duplicate functions and report them
    duplicate_functions = find_duplicate_function_definitions(input_script_path)
    report_duplicate_functions(duplicate_functions, input_script_path)


if __name__ == "__main__":
    main()
