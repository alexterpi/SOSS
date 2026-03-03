#!/bin/bash

# 1) Create the Python environment (julia-env) and install the necessary 
#    Python packages (numpy, ase, pymatgen).

echo ""
echo "1) Intalling Python packages"
echo ""

source install_python_packages.sh

# 2) Intall Julia packages and setup Python environment inside Julia.

echo ""
echo "2) Intalling Julia packages"
echo ""

julia --project=. install_julia_packages.jl

# 3) Test installation

echo ""
echo "3) Testing installation"
echo ""

julia --project=. test_installation.jl