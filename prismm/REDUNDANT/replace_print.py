import ast
import astunparse
import argparse
import subprocess
import sys
from general_functions import pretty_print

class PrintToPrettyPrint(ast.NodeTransformer):
    def visit_Expr(self, node):
        if isinstance(node.value, ast.Call) and isinstance(node.value.func, ast.Name) and node.value.func.id == "print":
            node.value.func.id = "pretty_print"
        return node

    def visit_FunctionDef(self, node):
        if node.name == "pretty_print":
            return None
        return node

def replace_print_statements(script_content):
    tree = ast.parse(script_content)
    transformer = PrintToPrettyPrint()
    modified_tree = transformer.visit(tree)

    modified_script = astunparse.unparse(modified_tree)
    modified_script = f"from general_functions import *\n{modified_script}"
    return modified_script

def run_all_scripts(test_case, max_number_of_solutions, max_default_path_length):
    python_files = [
        ("run_simulation.py", [str(test_case)]),
        ("run_summed_BP_search.py", [str(test_case), str(max_number_of_solutions), str(max_default_path_length)]),
        ("run_build_trees_and_timings.py", [str(test_case)]),
        ("run_model_evaluation.py", [str(test_case)]),
    ]

    for file_name, args in python_files:
        with open(file_name, 'r') as f:
            script_content = f.read()

        modified_script = replace_print_statements(script_content)
        output_file = file_name.rsplit('.', 1)[0] + '_2.py'

        with open(output_file, 'w') as f:
            f.write(modified_script)

        process = subprocess.run(["python", output_file] + args)
        if process.returncode != 0:
            print(f"Error in {output_file}. Exiting...")
            sys.exit(1)

if __name__ == "__main__":
    test_case = 51
    max_number_of_solutions = 10
    max_default_path_length = 20

    run_all_scripts(test_case, max_number_of_solutions, max_default_path_length)

