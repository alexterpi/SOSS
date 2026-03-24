SOSS
=================================
`SOSS` (Search for Optimally Stable Structures) is an in-house optimization package developed at BCAM (Basque Center for Applied Mathematics) and written in Julia. It is designed to search for ground-state configurations of solid structures using heuristic optimization techniques. The core algorithm exploits the characteristic features of solid-state ionic structures and processes to reduce a high-dimensional optimization problem into a discrete model with precomputed components of the objective function, enabling efficient evaluations.

Key Features: 

- Targeted applications: solid-state ionics (e.g. solid oxide fuel cells, solid-state batteries, supercapacitors, various energy storage devices, and sensors, etc.).
- Type of solid structures: the position of the sites in the lattice remain constant; only some atoms are rearranged in the lattice sites (those that generate vacancies and/or impurities). 
- Input/Output structure format: VASP format (POSCAR).
- Optimization solvers: Simulate Annealing (trajectory-based) **, Discrete Global-based Harmony Search (population-based)?**
- Objective function: Coulomb interaction energy with a cutoff radius **, Ewald summation?**
- Overlap function (optional): an additional repulsion term in the interaction energy to avoid undesired overlaps between ions.
- Ion-site constraint (optional): Classification of lattice sites based on their initial chemical element, enabling selective occupation rules for each ion type during optimization.

From a user perspective, the typical workflow of `SOSS` is as follows:

- Input

  - An initial atomic structure file.
  - A settings file, specifying the atomic species to be added or removed with respect to the initial structure, together with other optimization options.
  - A parameters file associated with the selected optimization solver.
- Optimization
- Output

  - Files that monitor the optimization process and summarize the final results.
  - The final structures that minimize the electrostatic interaction energy of the system. These structures correspond to optimal redistributions of atoms over the lattice sites of the input structure.

.. image:: Scheme.png
   :align: center
   :width: 800px

System Requirements
-------------------

- Julia (latest available version is recommended): `<https://julialang.org/>`_
- Python (latest available version is recommended): `<https://www.python.org/>`_

Documentation
-------------

The User manual (`SOSSManual.pdf`) included in the package provides a general introduction to `SOSS` and a detailed description of:

- The structure of the package (directories and files).
- The implemented computational modules.
- Installation procedures.
- Instructions for running test cases.
- General usage of the software.

Users are strongly encouraged to read the User manual before using `SOSS`.

For a more in-depth description of the theoretical background and practical applications of the code, refer to `“Screening Energetically Stable Structures in Solid-State Ionics Applications” <...>`_ by Alex Teruel et al.

Installation
--------------------

The `SOSS` package can be downloaded from its GitHub repository::

    git clone https://github.com/alexterpi/SOSS


Detailed installation instructions are provided in User manual.

Test Cases
----------

The package includes two test cases that generate optimized undoped (Test 1) and doped (Test 2) `LLZO <https://en.wikipedia.org/wiki/LLZO>`_ structures to demonstrate the main functionalities of `SOSS`.

- Test 1: Optimization without Dopants. In this test, the composition of the material remains unchanged. The objective is to identify the optimal atomic configuration of a specific species already present in the structure. The input structure corresponds to Li\ :sub:`960`\ La\ :sub:`192`\ Zr\ :sub:`128`\ O\ :sub:`768` (2x2x2 supercell of LLZO with all possible Li sites occupied). The optimization aims to obtain the structure Li\ :sub:`448`\ La\ :sub:`192`\ Zr\ :sub:`128`\ O\ :sub:`768` by selecting an energetically favorable distribution of 448 Li atoms across the 960 available Li sites.
- Test 2: Optimization with Aliovalent Doping. This test involves a modification of the composition by introducing dopants. Starting from the same initial structure as in Test 1, a subset of the original atoms is removed, and new dopant species are added. The goal is to obtain the doped structure Li\ :sub:`432`\ La\ :sub:`192`\ Zr\ :sub:`112`\ Ta\ :sub:`16`\ O\ :sub:`768`. During optimization, 432 Li atoms are distributed among the 960 initial Li sites, while 112 Zr and 16 Ta atoms are distributed among the 128 initial Zr sites.

These examples are intended to help users verify a correct installation and to serve as starting points for new applications. Instructions on how to run them and interpret the results are provided in User manual.

License
-------

`SOSS` is distributed under the **GNU General Public License v3.0**.

You are free to use, modify, and redistribute this software under the terms of the GPLv3. See the `LICENSE` file included with the package for the full license text.

Reference
---------

This package is distributed as part of a *Computer Physics Communications* (CPC) *Computer Programs in Physics* (CPiP) publication.

If you use `SOSS`, please cite `this paper <...>`_:


::

    @article{...
    }

Authors
-------

`SOSS` is being developed by:

- Alex Teruel, BCAM - Basque Center for Applied Mathematics, Spain
- Dr. Carlos León-Chinchay, University of Turku, Finland
- Dr. Mauricio Rincon Bonilla, BCAM - Basque Center for Applied Mathematics, Spain
- Dr. Henry Andres Cortes, BCAM - Basque Center for Applied Mathematics, Spain
- Prof. Elena Akhmatskaya, BCAM - Basque Center for Applied Mathematics, Spain

Contact
-------

| If you have questions, please don't hesitate to reach out at: ateruel@bcamath.org
| GitHub: `SOSS GitHub Repository <https://github.com/alexterpi/SOSS>`_
