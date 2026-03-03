#!/bin/bash

# input files directory
indir=$1

# number of cpus
ncpus=$2

# copy run SOSS script
cp ../SCRIPTS/run_soss.jl .

# start job
echo "\n-> $indir"

# copy input files
cp ../../INPUTS/$indir/* .

# run SOSS
t0=$(date +%s)
julia -t $ncpus --project=../../../INSTALL/. run_soss.jl
t1=$(date +%s)
dt=$((t1 - t0))

# join the structures and optimization files in a single (results) directory
rm -rf RESULTS
mkdir RESULTS
mv OPTFILES STRUCTURES RESULTS

# remove previous results (with the same name) from OUTPUTS directory
rm -rf ../../OUTPUTS/$indir

# move job results to OUTPUTS directory
mv RESULTS ../../OUTPUTS/$indir

# clean RUN directory
rm run_soss.jl
rm -f instruct.vasp settings.jl params.jl seeds.dat

# show success message
#echo "\n$indir was executed successfully\n"
#echo "-> The job took $(echo "scale=2; $dt/60" | bc) minutes, the results are in OUTPUTS/$indir"

# Check if the program ran correctly (this does not mean that the optimization results are correct,
# it simply means that the code did not fail during the optimization)
if [ -d "../../OUTPUTS/$indir/STRUCTURES" ]; then
  # show success message
  echo "\n$indir was executed successfully\n"
  echo "-> The job took $(echo "scale=2; $dt/60" | bc) minutes, the results are in OUTPUTS/$indir\n"
else
  # show failure message
  echo "\nAn error occurred during the execution of SOSS:\n"
  echo "* If the error occurred while running any of the test cases (TEST_1, TEST_2) without modifying the input files,\nSOSS was likely not installed correctly. Please verify that the installation was completed successfully.\n"
  echo "* If the test cases ran successfully and the error appears during a user-defined optimization, there is likely\nan issue with one or more input files (instruct.vasp, settings.jl, params.jl, seeds.dat (optional)). Please\ncheck carefully that all provided input files are correct.\n"
  echo "* If you are unable to identify the problem, you may contact one of the corresponding authors listed in the\nmanuscript (SOSSManuscript.pdf) to verify that the current version of SOSS is functioning correctly. Please ensure that\nthe issue is not due to user input or system requirements/configurations before contacting the authors.\n"
fi