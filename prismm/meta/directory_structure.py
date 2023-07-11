import os
import sys

def list_files(startpath):
    for root, dirs, files in os.walk(startpath):
        level = root.replace(startpath, '').count(os.sep)
        indent = ' ' * 4 * (level)
        subindent = ' ' * 4 * (level + 1)
        
        python_files = [f for f in files if f.endswith('.py') and not f.endswith('_IO.py')]
        
        if python_files:
            print('{}{}/'.format(indent, os.path.basename(root)))
            for f in python_files:
                print('{}{}'.format(subindent, f))

if len(sys.argv) != 2:
    print("Please provide the directory path.")
else:
    dir_path = sys.argv[1]
    list_files(dir_path)
