# Import Julia packages
try
    using Conda
    using PyCall
    using UnPack
    using Random
    using StatsBase
    using Statistics
    using MutableNamedTuples
    using PyPlot
    using DelimitedFiles
    using Test
    println("Exit to import Julia packages")
catch e
    println("Failed to import Julia packages:")
    println(e)
end

# Check that Python environment is working
importlib_metadata = pyimport("importlib.metadata")

try
    numpy = pyimport("numpy")
    println("Exit to import numpy / numpy version: ", importlib_metadata.version("numpy"))
    global numpy_exit = 1
catch e
    println("Failed to import numpy: ", e)
    global numpy_exit = 0
end

try
    ase = pyimport("ase")
    println("Exit to import ase / ase version: ", importlib_metadata.version("ase"))
    global ase_exit = 1
catch e
    println("Failed to import ase: ", e)
    global ase_exit = 0
end

try
    pymatgen = pyimport("pymatgen")
    println("Exit to import pymatgen / pymatgen version: ", importlib_metadata.version("pymatgen"))
    global pymatgen_exit = 1
catch e
    println("Failed to import pymatgen: ", e)
    global pymatgen_exit = 0
end

if numpy_exit + ase_exit + pymatgen_exit == 3
    println("Python environment properly working")
    println("SOSS was successfully installed\n")
else
    println("An error occurred during the installation process.")
    println("Please review the error messages shown on the screen to verify that the issue is not related to your local installation or configuration of Python, Julia, and/or Conda.")
    println("If you are unable to identify the problem (or if you determine that it is an issue with SOSS itself), you may contact one of the corresponding authors listed in the manuscript (SOSSManuscript.pdf) to verify the correct functioning of SOSS.\n")
    println("Please do not contact the authors before confirming that the issue is not due to a user installation or configuration problem.")
end