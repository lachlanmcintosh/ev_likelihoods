import os
import ast
import re
import argparse
import logging
from typing import Dict, List, Tuple, Union
from pathlib import Path

# Type aliases for readability
SummaryType = Dict[str, Union[Dict[str, List[str]], List[Tuple[int, str]]]]

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class CodeAnalyzer(ast.NodeVisitor):
    """
    Abstract Syntax Tree (AST) Node Visitor that summarizes Python classes, functions and imports in a file.
    This class overrides the methods of ast.NodeVisitor to visit the functions and classes in a Python file.
    """

    def __init__(self) -> None:
        self.summary: SummaryType = {
            'classes': {},
            'functions': [],
            'imports': []
        }

    def _add_function(self, node: Union[ast.FunctionDef, ast.AsyncFunctionDef]) -> None:
        """
        Add function or async function name to summary and visit child nodes.
        If function starts with 'test_', it will be ignored.
        """
        if not node.name.startswith('test_'):  # Ignore functions starting with 'test_'
            self.summary['functions'].append(node.name)
        self.generic_visit(node)

    def visit_FunctionDef(self, node: ast.FunctionDef) -> None:
        self._add_function(node)

    def visit_AsyncFunctionDef(self, node: ast.AsyncFunctionDef) -> None:
        self._add_function(node)

    def visit_ClassDef(self, node: ast.ClassDef) -> None:
        """
        Add class and its methods to summary and visit child nodes.
        """
        methods = [
            f.name
            for f in node.body
            if isinstance(f, (ast.FunctionDef, ast.AsyncFunctionDef))
            and not f.name.startswith('test_')
        ]
        if methods:
            self.summary['classes'][node.name] = methods
        self.generic_visit(node)


def analyze_python_file(filepath: Path) -> SummaryType:
    """
    Analyzes a Python file and returns a summary of its classes, functions and imports.
    """
    logger.debug(f'Analyzing file: {filepath}')
    data = ""
    import_lines = []
    import_pattern = re.compile(r"^\s*(import|from)\s+[\w\.]+\s*(import)?\s*[\w\.\*\,\(\) ]+$")
    try:
        with open(filepath, 'r') as file:
            for line_number, line in enumerate(file, start=1):
                data += line
                if import_pattern.match(line):
                    import_lines.append((line_number, line.strip()))
    except FileNotFoundError:
        logger.error(f"File {filepath} not found.")
        return {'classes': {}, 'functions': [], 'imports': []}

    try:
        tree = ast.parse(data)
    except SyntaxError as e:
        logger.error(f"Failed to parse {filepath}: {str(e)}")
        return {'classes': {}, 'functions': [], 'imports': []}

    analyzer = CodeAnalyzer()
    analyzer.visit(tree)

    # Add exact lines of import statements
    analyzer.summary['imports'] = import_lines

    return analyzer.summary


def summarize_directories(dirpaths: List[str]) -> Dict[str, Dict[str, SummaryType]]:
    """
    Returns a summary of the classes and functions in all Python files in given directories.
    """
    summary = {}
    for dirpath in dirpaths:
        dirpath_obj = Path(dirpath)
        if not dirpath_obj.exists() or not dirpath_obj.is_dir():
            logger.error(f"{dirpath} is not a valid directory.")
            continue

        root_summary = {
            str(filepath): analyze_python_file(filepath)
            for filepath in dirpath_obj.rglob('*.py')
        }
        if root_summary:  # Ignore empty directories
            summary[dirpath] = root_summary

    return summary


def print_summary(summary: Dict[str, Dict[str, SummaryType]], track_imports: bool) -> None:
    """
    Prints a nicely formatted summary of the classes, functions and imports in Python files.
    """
    print("For context in debugging the below error here is a nicely formatted summary of the classes, functions and all import statements in the module:")
    for dirpath, dir_summary in summary.items():
        print(f"\nDirectory: {dirpath} has the files:")
        for filepath, data in dir_summary.items():
            print(f"   File: {os.path.basename(filepath)}", end="")
            if data['functions']:
                print(f" | Functions: {', '.join(data['functions'])}", end="")
            if data['classes']:
                print(" | Classes: ", end="")
                for classname, methods in data['classes'].items():
                    print(f"{classname} (Methods: {', '.join(methods)}) ", end="")
            if track_imports and data['imports']:
                print(" | Imports: ", end="")
                for line_number, import_line in data['imports']:
                    print(f"Line {line_number}: {import_line}, ", end="")
            print()
        print('----------')


import json

def main() -> None:
    """
    Main function that parses command line arguments and summarizes Python files in given directories.
    """
    parser = argparse.ArgumentParser(description='Summarize Python files in directories.')
    parser.add_argument('dirpaths', type=str, nargs='+', help='The directories to summarize')
    parser.add_argument('--track_imports', action='store_true', help='Flag to track import statements')
    args = parser.parse_args()

    summary = summarize_directories(args.dirpaths)

    # Save summary to a JSON file
    with open('summary.json', 'w') as f:
        json.dump(summary, f, indent=4)

if __name__ == "__main__":
    main()
