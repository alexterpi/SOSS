# Optimization Solver
solver = "SANSS"     # Only SANSS solver available.

# Input and Output Atomic Structure Files Format
IN_format  = "VASP"  # Only VASP format available.
OUT_format = "VASP"  # Only VASP format available.

# Note:
# The input atomic structure (instruct.vasp) contains 15 Li atoms per FU
# (Li15 La3 Zr2 O12), but only 7 are actually occupied in the material.
# We use a structure with 15 Li atoms, allowing us to define all the possible
# positions that Li can occupy. However, when determining the initial structure,
# we must indicate the actual number of Li atoms that the structure contains.
# For this reason, we set 7 Li átoms in FU_numbers.

# Initial Structure: Li7 La3 Zr2 O12 (the same as in Test 1)
FU_symbols = ["Li", "La", "Zr", "O"]  # Names (chemical symbols).
FU_numbers = [7, 3, 2, 12]            # Number of atoms and host ions per formula unit.
FU_charges = [1.0, 3.0, 4.0, -2.0]    # Charges (in e units).
FU_protons = [3, 57, 40, 8]           # Atomic numbers.
FU_copies  = 64                       # Number of formula units in the supercell.

# Host ions
Ihost_symbols = ["Li", "Zr"]  # Use Li and Zr sites to perform the optimization.
                              # Li sites define subset 1.
                              # Zr sites define subset 2.
Ihost_numbers = [6.75, 1.75]  # The final structure contains 6.75 Li atoms and 1.75 Zr atoms per FU.

# Added ions / dopants
Iadd_symbols = ["Ta"] # Add Ta ions (Ta-doped LLZO).
Iadd_numbers = [0.25] # The final structure contains 0.25 Ta ions per FU.
Iadd_charges = [5.0]  # Ta charge (in e units).
Iadd_protons = [73]   # Ta atomic number.
Iadd_sites   = [2]    # It is known that Ta ions can only occupy Zr sites -> subset 2.

# Ion-Site constraint
IonSiteConstraint = true   # The constraint must be activated (subsets are needed).

# Cutoff and overlap radius (Å)
r_cutoff  = 7.0   # A cutoff radius of 7 Å proved to work well for this system.
r_overlap = 1.7   # It is known (for this system) that distances below 1.7 Å are forbidden for ions.

# Summary of the optimization:
# The input Atomic Structure File contains a 2x2x2 supercell of Li15 La3 Zr2 O12.
# The Li ions (6.75 per FU) will be distributed among the Li sites (15 per FU).
# The Zr and Ta ions (1.75 and 0.25 per FU) will be distributed among Zr sites (2 per FU).
# The Output Atomic Structure File will contain a 2x2x2 supercell of Li6.25 La3 Zr1.75 Ta0.25 O12.

#===================================================================================================
INFORMATION

- SOLVER

* solver: SANSS
Method used to solve the optimization problem.
Optimization solver.

- IN/OUT DATA

* IN_format: VASP
Format of the input atomic structure file.

* OUT_format: VASP
Format of the output atomic structure file.

- SYSTEM

Note: When we refer to "atoms and host ions in the initial structure",
      we are referring to all the chemical species that are in it.

* FU_symbols: vector of strings
Names (chemical symbols) of atoms and host ions in the initial structure.

* FU_numbers: vector of floats
Number of atoms and host ions (per formula unit) in the initial structure.

* FU_charges: vector of floats
Charges of atoms and host ions in the initial structure.

* FU_protons: vector of integers
Atomic numbers of atoms and host ions in the initial structure.

* FU_copies: integer
Number of formula units contained in the supercell.

* Ihost_symbols: vector of strings
Names (chemical symbols) of host ions (chemical species that are in the
initial structure and can modify their positions during the optimization).

* Ihost_numbers: vector of floats
Number of host ions (per formula unit) in the final structure.

* Iadd_symbols: vector of strings
Names (chemical symbols) of added ions/dopants (chemical species that
are not in the initial structure).

* Iadd_numbers: vector of floats
Number of added ions (per formula unit) in the final structure.

* Iadd_charges: vector of floats
Charges of added ions.

* Iadd_protons: vector of integers
Atomic numbers of added ions.

* Iadd_sites: vector of integers
Subset of host ion sites (following the order defined in Ihost_symbols)
available for added ions.

* IonSiteConstraint: bool
Boolean variable to control the ion-site constraint functionality.
If true, the vector Iadd_sites is used.
If false, the vector Iadd_sites is ignored.

* r_cutoff: float
Cutoff radius of the Coulomb interaction energy.
If it is positive, the code assumes that the entered value is the cutoff in Å.
If it is negative, the code assumes that the entered value is a fraction of the max cutoff.
The max cutoff is defined as half the length of the shortest side of the unit cell.

* r_overlap: float
Overlap radius in Å.
Configurations with ion-ion or ion-atom distances below r_overlap are rejected.
===================================================================================================#