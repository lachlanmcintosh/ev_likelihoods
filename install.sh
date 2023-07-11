#!/bin/bash

# navigate to the directory where setup.py is located
cd ~/PRISMM

# uninstall the current version
pip uninstall -y prismm

# install the updated version
pip install .

# navigate back to the original directory
cd -

