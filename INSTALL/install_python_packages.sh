#!/bin/bash

# Name of the new Conda environment
ENV_NAME="julia-env"

# Create the environment with the latest Python version
conda create -y -n "$ENV_NAME" python

# Activate Conda base so we can use 'conda activate' reliably in scripts
source "$(conda info --base)/etc/profile.d/conda.sh"

# Activate the newly created environment
conda activate "$ENV_NAME"

# Show Python version
echo "Environment '$ENV_NAME' created with Python version:"
python --version

# Install packages
echo ""
echo "Intalling Python packages"
pip install numpy ase pymatgen
echo ""

# Get environment path
ENV_PATH=$(which python)
echo "$ENV_PATH" > python_path.txt

#echo "Numpy, ASE, and Pymatgen packages installed in $ENV_NAME Environment."
# Information for the user
echo "If no errors were reported:"
echo " -> Numpy, ASE, and Pymatgen packages are installed in the $ENV_NAME environment."
echo " -> The $ENV_NAME environment location can be found in the file python_path.txt."
echo ""
echo "If errors were reported:"
echo " -> Verify that Python and Conda are correctly installed and in the User PATH."