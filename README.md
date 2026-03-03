# SOSS

`SOSS` (Search for Optimal Structure) is an in-house optimization package developed at BCAM (Basque Center for Applied Mathematics) and written in Julia. It is designed to search for ground-state configurations of solid structures using heuristic optimization techniques. The core algorithm exploits the characteristic features of solid-state ionic structures and processes to reduce a high-dimensional optimization problem into a discrete model with precomputed components of the objective function, enabling efficient evaluations.

From a user perspective, the typical workflow of `SOSS` is as follows:

- An atomic structure file is provided as input.
- A settings file is provided as input, specifying the atomic species to be added or removed with respect to the initial structure, together with other optimization options.
- A parameters file associated with the selected optimization solver is provided as input.
- The code is executed.
- The code returns (i) output files that monitor the optimization process and summarize the final results, and (ii) the final structures that minimize the interaction energy of the system. These structures correspond to optimal redistributions of atoms over the lattice sites of the input structure.

For a detailed description of the package, the accompanying user manual (`SOSSManual.pdf`) provides a general introduction to `SOSS` and a complete description of:

- The structure of the package (directories and files)
- The implemented computational modules
- Installation procedures
- Instructions for running test cases
- General usage of the software

Users are strongly encouraged to read the user manual before using `SOSS`.

For a more in-depth description of the theoretical background and practical applications of the code, refer to the manuscript *“Screening Energetically Stable Structures in Solid-State Ionics Applications”* by Alex Teruel et al. (`SOSSManuscript.pdf`).

## Contents

A detailed description of the contents of the three directories that constitute the package (`INSTALL`, `SRC`, and `WORK`) is provided in section 2 of the user manual.

## Computational Modules Documentation

The primary documentation for the `SOSS` computational modules can be found in section 3 of the user manual and in the source code itself (`SRC` directory).

## Installation

Detailed installation instructions and system requirements are provided in sections 4.1, 4.2 and 4.3 of the user manual.

## Test Cases

The package includes two test cases that demonstrate the main functionalities of `SOSS`. These examples are intended to help users verify a correct installation and to serve as starting points for new applications.

Instructions on how to run these examples and interpret the results are provided in section 4.4 of the user manual.

## Usage

The general usage of `SOSS` (input files generation, code execution, and results interpretation) is provided in section 4.5 of the user manual.

## License

`SOSS` is distributed under the **GNU General Public License v3.0**.

You are free to use, modify, and redistribute this software under the terms of the GPLv3.
See the `LICENSE` file for the full license text.

## Reference

This package is distributed as part of a *Computer Physics Communications* (CPC) *Computer Programs in Physics* (CPiP) publication.

If you use `SOSS`, please cite this paper:

-> add reference here <-