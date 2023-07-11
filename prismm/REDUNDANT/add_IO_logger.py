import sys
import re
from functools import wraps
from contextlib import contextmanager

LOG_INPUT_OUTPUT_DECORATOR = '''
from functools import wraps
import os

def log_input_output(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        result = func(*args, **kwargs)

        if not os.path.exists('IO_LOGGER'):
            os.mkdir('IO_LOGGER')
        with open(f'IO_LOGGER/{func.__name__}_io_log.txt', 'a') as f:
            f.write(f'Input: {repr(args)}\\n')
            f.write(f'Output: {repr(result)}\\n')

        return result
    return wrapper
'''


@contextmanager
def open_file(file, mode):
    f = open(file, mode)
    try:
        yield f
    finally:
        f.close()


def process_code(code):
    code_lines = code.split('\n')

    # Add input-output logging
    modified_lines = [LOG_INPUT_OUTPUT_DECORATOR]
    for line in code_lines:
        match = re.match(r'^def (\w+)\((.*)\):$', line)
        if match:
            name, args = match.groups()
            modified_lines.append(f'@log_input_output\n{line}')
        else:
            modified_lines.append(line)

    return '\n'.join(modified_lines)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python save_input_output.py input_file output_file")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    with open_file(input_file, 'r') as f:
        code = f.read()

    modified_code = process_code(code)

    with open_file(output_file, 'w') as f:
        f.write(modified_code)

