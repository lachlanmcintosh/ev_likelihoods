import argparse
import os
import glob
import ast
import astor

preamble = """Your primary task is to create a robust testing suite using the pytest library for the existing codebase. The guidelines for this process are as follows:

Pytest Adoption: Adopt pytest as the testing framework for the suite. Familiarize yourself with its features like fixtures and parameterization to write effective tests.

Comprehensive Coverage: With pytest, develop tests for each function in the codebase, covering not only typical usage but also edge and failure cases.

Error & Exception Handling: Tests should validate function inputs and outputs, ensuring correct handling of errors and appropriate exceptions are thrown when necessary.

Scenario Coverage: Create diverse and exhaustive tests, considering a variety of real-world scenarios that the code may encounter.

Assert Statements: Utilize pytest's powerful assert statements for comparing expected and actual outcomes.

Directory Structure and Naming Convention: Maintain a structure where for each module in the codebase, a corresponding 'tests' subdirectory exists. Test files should align with the module they are testing, prefixed with 'test_', e.g., 'test_mymodule.py' for 'mymodule.py'.

The objective is to enhance the reliability and maintainability of the codebase through rigorous testing with pytest. Document each test clearly, detailing its purpose and the particular scenario it's assessing. Your proficiency in crafting test cases is crucial for this task.

Here is the code for the module, followed by any deprecated testing functions if they exist and then followed by typical arguments and their expected output if applicable"""


def print_node_code(node):
    return astor.to_source(node)

def find_corresponding_test_function(function_name, test_functions):
    test_function_name = "test_" + function_name
    if test_function_name in test_functions:
        return print_node_code(test_functions[test_function_name])
    else:
        return f"Test function not found for function: {function_name}\n"

def find_io_logs(function_name, io_log_dir):
    result = ""
    log_file_path = f"{io_log_dir}/{function_name}_io_log.txt"
    log_files = glob.glob(log_file_path)
    if not log_files:
        result += f"No IO log file found for function: {function_name}. Searched path: {log_file_path}\n"

    for log_file in log_files:
        result += f"Log file: {log_file}\n"
        with open(log_file, 'r') as file:
            lines = file.readlines()[:8]
            for line in lines:
                result += line.strip() + "\n"

    return result

def print_code_tests_io(module_filepath, test_module_filepath, io_log_dir):
    result = ""

    with open(module_filepath, 'r') as module_file, open(test_module_filepath, 'r') as test_file:
        module_tree = ast.parse(module_file.read())
        test_module_tree = ast.parse(test_file.read())

    module_functions = {node.name: node for node in ast.walk(module_tree) if isinstance(node, ast.FunctionDef)}
    test_functions = {node.name: node for node in ast.walk(test_module_tree) if isinstance(node, ast.FunctionDef)}

    for function_name, function_node in module_functions.items():
        function_code = print_node_code(function_node)
        test_function_code = find_corresponding_test_function(function_name, test_functions)
        io_logs = find_io_logs(function_name, io_log_dir)

        result += preamble + "\n"  # Print preamble before each function
        result += f"Function {function_name}:\n{function_code}\n"
        result += f"Corresponding test function:\n{test_function_code}\n"
        result += f"IO logs:\n{io_logs}\n"
        result += "=" * 80 + "\n"  # horizontal line to separate functions

    return result

def main():
    parser = argparse.ArgumentParser(description='Analyze a Python module.')
    parser.add_argument('top_dir', help='The top directory where Python modules reside')

    args = parser.parse_args()

    top_dir = args.top_dir  # Get the top directory where the modules reside
    test_module_filepath = "/home/users/allstaff/lmcintosh/clonal_trees/clonal_trees/run_simulation/tests.py"  # Filepath to your test module, replace if different
    io_log_dir = "/home/users/allstaff/lmcintosh/clonal_trees/IO_LOGGER"  # IO logs directory

    for dirpath, dirnames, filenames in os.walk(top_dir):
        if 'tests' not in dirpath:  # Exclude directories named "tests"
            for filename in filenames:
                if filename.endswith(".py") and not filename.endswith("_IO.py") and not filename.startswith("test_"):
                    module_filepath = os.path.join(dirpath, filename)

                    # Prepare the report
                    report = f"Module: {filename}\n"
                    report += "=" * 80 + "\n"  # horizontal line to separate sections

                    # Extract code, tests and IO information for the module
                    report += print_code_tests_io(module_filepath, test_module_filepath, io_log_dir)

                    # Write the report to a file
                    test_dir = os.path.join(dirpath, 'tests')  # Define the tests directory path
                    os.makedirs(test_dir, exist_ok=True)  # Create the tests directory if it does not exist
                    analysis_file_path = os.path.join(test_dir, f"{filename}_analysis.txt")
                    with open(analysis_file_path, 'w') as analysis_file:
                        analysis_file.write(report)

                    print(f"Wrote report to: {analysis_file_path}")  # Added line

if __name__ == "__main__":
    main()
