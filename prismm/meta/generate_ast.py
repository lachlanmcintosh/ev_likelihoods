"""
Python Function Calls Analyzer

This script performs static analysis on a provided Python script or a directory containing Python scripts to identify all function definitions, line numbers where they start and end, and the functions they call.

The script uses Python's `ast` (Abstract Syntax Trees) module to parse the code and walks the syntax tree to visit function definitions (`FunctionDef`) and function calls (`Call`). The analyzer creates a dictionary where each function is mapped to the line numbers of its definition and a set of functions it calls.

The function information is sorted based on the number of function calls they make in descending order and then printed to the console.

Usage:
    python <script_name.py> <target_file_or_directory>

    where:
    - <script_name.py> is the name of this script
    - <target_file_or_directory> is the name of the Python script or directory to analyze

Output:
    The script prints the name of each function, line numbers where they start and end, and the functions they call.
"""
import ast
import argparse
import os

class FunctionCallVisitor(ast.NodeVisitor):
    def __init__(self):
        self.func_map = dict()
        self.current_function = None

    def visit_FunctionDef(self, node):
        self.current_function = node.name
        self.func_map.setdefault(node.name, {'start_lineno': node.lineno, 'end_lineno': node.body[-1].lineno, 'calls': set()})
        self.generic_visit(node)
        self.current_function = None

    def visit_Call(self, node):
        if isinstance(node.func, ast.Name) and node.func.id in self.func_map and self.current_function is not None:
            # Add to the set of called functions for the current function
            self.func_map[self.current_function]['calls'].add(node.func.id)
        self.generic_visit(node)

def generate_func_map(py_file):
    with open(py_file, 'r') as file:
        root = ast.parse(file.read())

    visitor = FunctionCallVisitor()
    visitor.visit(root)

    sorted_func_map = sorted(visitor.func_map.items(), key=lambda item: len(item[1]['calls']), reverse=True)
    for func, details in sorted_func_map:
        print(f'Function: {func} (start line: {details["start_lineno"]}, end line: {details["end_lineno"]})')
        print(f'    Calls: {", ".join(details["calls"]) if details["calls"] else "None"}')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('path', help='The Python script or directory to analyze')
    args = parser.parse_args()

    if os.path.isfile(args.path):
        generate_func_map(args.path)
    elif os.path.isdir(args.path):
        for dirpath, dirnames, filenames in os.walk(args.path):
            for filename in filenames:
                if filename.endswith('.py'):
                    print(f'Analyzing {filename}:')
                    generate_func_map(os.path.join(dirpath, filename))
    else:
        print('The path provided does not exist')

if __name__ == '__main__':
    main()
