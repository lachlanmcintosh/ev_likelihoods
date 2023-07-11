"""
Python Unused Function Finder

This script performs static analysis on a provided Python script to identify functions that are defined but not called anywhere in the codebase.

The script uses Python's `ast` (Abstract Syntax Trees) module to parse the code and walks the syntax tree to visit function definitions (`FunctionDef`) and function calls (`Call`). It tracks the defined and called functions in two separate sets, and then finds the difference between the two sets to identify unused functions.

Usage:
    python <script_name.py> <target_python_file.py>

    where:
    - <script_name.py> is the name of this script
    - <target_python_file.py> is the name of the Python script to analyze

Output:
    The script prints the names of unused functions, if any, to the console.
"""
import ast
import argparse

class FunctionVisitor(ast.NodeVisitor):
    def __init__(self):
        self.defined = set()
        self.called = set()

    def visit_FunctionDef(self, node):
        self.defined.add(node.name)
        self.generic_visit(node)

    def visit_Call(self, node):
        if isinstance(node.func, ast.Name):
            self.called.add(node.func.id)
        self.generic_visit(node)

def find_unused_functions(py_file):
    with open(py_file, 'r') as file:
        root = ast.parse(file.read())

    visitor = FunctionVisitor()
    visitor.visit(root)

    return visitor.defined - visitor.called

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('script', help='The Python script to analyze')
    args = parser.parse_args()

    unused_funcs = find_unused_functions(args.script)
    for func in unused_funcs:
        print(f'Unused function: {func}')

if __name__ == '__main__':
    main()
