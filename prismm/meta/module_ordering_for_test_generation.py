import os
import sys
import ast
from typing import List, Tuple, Union, Generator
import logging

# Setup logger
logging.basicConfig(level=logging.INFO)

class DependencyCounter:
    """Class to manage counting of dependencies in python modules."""
    
    def __init__(self, base_dir_name: str):
        """Initialize with base directory name for the package being developed."""
        self.base_dir_name = base_dir_name
    
    def count_dependencies(self, module_path: str) -> Union[int, float]:
        """
        Count the number of dependencies in a module.

        :param module_path: Path to the python module.
        :return: Number of dependencies or infinity if the module couldn't be loaded.
        """
        try:
            with open(module_path, "r") as source:
                tree = ast.parse(source.read())
            dependencies = 0
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    dependencies += sum(self.base_dir_name in name.name for name in node.names)
                elif isinstance(node, ast.ImportFrom):
                    if node.module and self.base_dir_name in node.module:
                        dependencies += 1
            return dependencies
        except Exception as e:
            logging.error(f"Couldn't load module {module_path} due to: {str(e)}")
            return float('inf')  # Assign maximum dependencies to modules that cannot be loaded

class PythonFileHandler:
    """Class to manage python files in a directory."""

    def __init__(self, dir_to_search: str):
        """Initialize with directory to search."""
        self.dir_to_search = dir_to_search

    def get_py_files(self) -> Generator[str, None, None]:
        """
        Generator to get python files in the directory.

        :return: Full path to a python file.
        """
        # walk through the directory, include subdirectories
        for dirpath, dirnames, filenames in os.walk(self.dir_to_search):
            for filename in filenames:
                if filename.endswith('.py') and not filename.endswith('_IO.py') and not filename.startswith("test"):
                    yield os.path.join(dirpath, filename)

class TestFileHandler:
    """Class to manage test files for python modules."""

    @staticmethod
    def has_test(module_path: str) -> bool:
        """
        Check if a test file exists for the module.

        :param module_path: Path to the python module.
        :return: True if a test file exists, False otherwise.
        """
        base_dir = os.path.dirname(module_path)
        test_dir = os.path.join(base_dir, 'tests')
        test_file_name = 'test_' + os.path.basename(module_path)
        test_file_path = os.path.join(test_dir, test_file_name)
        return os.path.isfile(test_file_path)

    @staticmethod
    def create_test_dir_if_not_exists(module_path: str) -> None:
        """
        Create a test directory if it does not exist.

        :param module_path: Path to the python module.
        """
        base_dir = os.path.dirname(module_path)
        test_dir = os.path.join(base_dir, 'tests')
        if not os.path.exists(test_dir):
            os.makedirs(test_dir)

def main(dir_to_search: str, base_dir_name: str):
    """Main function to execute the script."""

    io_logger_dir = os.path.join(base_dir_name, "IO_LOGGER")
    file_handler = PythonFileHandler(dir_to_search)
    dependency_counter = DependencyCounter(base_dir_name)
    test_file_handler = TestFileHandler()

    py_files = list(file_handler.get_py_files())
    logging.info(f"Found {len(py_files)} python modules")

    module_dependencies = [
        (
            f,
            dependency_counter.count_dependencies(f),
            get_io_function_names_and_logger_files(f)
        ) 
        for f in py_files if not test_file_handler.has_test(f)
    ]

    logging.info(f"Found {len(module_dependencies)} modules without tests")

    module_dependencies.sort(key=lambda x: x[1], reverse=True)

    for module, dependencies, io_logger_files in module_dependencies:
        if dependencies != float('inf'):  # Skip modules that couldn't be loaded
            logging.info(f'{module}: {dependencies} dependencies, IO_LOGGER files: {io_logger_files}')
            test_file_handler.create_test_dir_if_not_exists(module)

if __name__ == "__main__":
    dir_to_search = sys.argv[1]
    base_dir_name = sys.argv[2] if len(sys.argv) > 2 else "clonal_trees"
    main(dir_to_search, base_dir_name)
