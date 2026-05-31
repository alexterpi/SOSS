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

# Initial Structure: Li7 La3 Zr2 O12
FU_symbols = ["Li", "La", "Zr", "O"]  # Names (chemical symbols).
FU_numbers = [7.0, 3.0, 2.0, 12.0]    # Number of atoms and host ions per formula unit.
FU_charges = [1.0, 3.0, 4.0, -2.0]    # Charges (in e units).
FU_protons = [3, 57, 40, 8]           # Atomic numbers.
FU_copies  = 64                       # Number of formula units in the supercell.

# Host ions
Ihost_symbols = ["Li"]    # Use Li sites to perform the optimization.
Ihost_numbers = [7.0]     # The final structure contains 7.0 Li atoms per FU.

# Added ions / dopants
Iadd_symbols = []   # No dopants.
Iadd_numbers = []
Iadd_charges = []
Iadd_protons = []
Iadd_sites   = []

# Ion-Site constraint
IonSiteConstraint = false   # Ion-Site constraint is not used.

# Energy model
LongRange  = [1, [7.0]] # Coulomb potential with a cutoff radius of 7.0 Å.
ShortRange = [1, [1.7]] # Overlap potential with an overlap radius of 1.7 Å.

# Summary of the optimization:
# The input Atomic Structure File contains a 2x2x2 supercell of Li15 La3 Zr2 O12.
# The Li ions (7 per Fu) will be distributed among the Li sites (15 per FU).
# The output Atomic Structure File will contain a 2x2x2 supercell of Li7 La3 Zr2 O12.

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

* FU_symbols: vector of strings
Names (chemical symbols) of anchored and host ions in the initial structure (chemical species in
the initial structure).

* FU_numbers: vector of floats
Number of anchored and host ions (per formula unit) in the initial structure.

* FU_charges: vector of floats
Charges of anchored and host ions in the initial structure.

* FU_protons: vector of integers
Atomic numbers of anchored and host ions in the initial structure.

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

- CONSTRAINTS

* Iadd_sites: vector of integers
Subset of host ion sites (following the order defined in Ihost_symbols) available for added ions.

* IonSiteConstraint: bool
Boolean variable to control the ion-site constraint functionality.
If true, the vector Iadd_sites is used.
If false, the vector Iadd_sites is ignored.

- ENERGY MODEL

* LongRange  = [LR_type, LR_params], LR_type: integer, LR_params: vector of floats
Vector that controls the long-range interactions model.
Options:
(1) Truncated Coulomb with a cutoff radius r_cutoff (Å)-> LR_type = 1, LR_params = [r_cutoff]
(2) Ewald summation -> LR_type = 2, LR_params = []

* ShortRange  = [SR_type, SR_params], SR_type: integer, SR_params: vector of floats
Vector that controls the short-range interactions model.
Oprions:
(1) No interaction -> SR_type = 0, SR_params = []
(2) Overlap repulsion with an overlap radius r_overlap (Å) -> SR_type = 0, SR_params = [r_overlap]
===================================================================================================#