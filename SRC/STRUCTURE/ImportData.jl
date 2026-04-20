module ImportData

# This module contains the functions used to read the inputs (structure file + settings tuple) and
# generate the data tuple (where the problem instance is defined). 
# The main function is get_data(...), which is defined at the end of the module and is extensively
# commented.
# To understand how the code works, we recommend reading the get_data(...) function and then (for 
# more details) the functions called from there. 
#
# Note: Currently, only VASP format is available.
#
# Note: 1D array -> Vector
#       2D array -> Matrix

# Import Julia modules
using PythonCall
using UnPack
using Random

# Import energy module
include("Energy.jl")
using .Energy

# Define references to Python modules/functions
const ase_io = Ref{Py}(Py(nothing))
const pymatgen_ase = Ref{Py}(Py(nothing))
const np = Ref{Py}(Py(nothing))

# Initialize Python environment
function __init__()
    ase_io[] = pyimport("ase.io.vasp")
    pymatgen_ase[] = pyimport("pymatgen.io.ase")
    np[] = pyimport("numpy")
end

function read_structure(format::String)

    # Function that reads the structure file and return the structure object.
    #
    # Arguments:
    # format -> String containing the format of the structure file.
    #
    # Returns:
    # structure -> ASE structure object.
    
    if format == "VASP"
        structure = ase_io[].read_vasp("instruct.vasp")
        return structure
    else
        println("Format not defined")
    end
end

function whereInt(elements::Array{Int,1}, elementToFind::Int)

    # Function that returns the indexes and number of repetitions of a given integer inside a
    # vector.
    #
    # Arguments:
    # elements ------> Vector of integers.
    # elementToFind -> Integer value we want to find inside the elements vector.
    #
    # Returns:
    # listOfIndexes --------> Vector containing the positions where the integer value was
    #                         found.
    # numberOfCoincidences -> Number of times the element was found.

    listOfIndexes        = findall(x->x==elementToFind, elements)
    numberOfCoincidences = length(listOfIndexes)

    return listOfIndexes, numberOfCoincidences
end

function get_species(symbols)

    # Function that gets a vector containing the chemical symbols of the ions in the initial
    # structure.
    #
    # Arguments:
    # symbols -> Vector containing the chemical symbols of the ions in the initial structure.
    #
    # Returns:
    # species -> Vector containing the chemical symbol of the species that appear in the initial
    #            structure.
    #            Note: It only indicates which species appear in the structure, not which species
    #            are placed at each lattice site.

    Ns = length(symbols)
    species = [symbols[1]]
    for symbol in symbols[2:Ns]
        if symbol ∉ species
            push!(species, symbol)
        end
    end
    
    return species
end

function get_protons(FU_symbols, FU_protons, Ihost_symbols)

    # Function that gets the atomic number (number of protons) of the host ions (interchangeable
    # ions that appear in the initial structure, which are partially removed from the initial
    # structure).
    #
    # Arguments:
    # FU_symbols ----> Vector containing the chemical symbols of the species that appear in the
    #                  initial structure (anchored and host ions).
    #                  Note: Anchored ions are the ions present in the initial structure, whose
    #                  position remains fixed during optimization (immobile species).
    # FU_protons ----> Vector containing the atomic numbers of the species in the initial structure
    #                  (anchored and host ions).
    # Ihost_symbols -> Vector containing the chemical symbols of the host ions.
    #
    # Returns:
    # Ihost_protons -> Vector containing the atomic numbers of the host ions.

    protonsOFspecie = Dict(zip(FU_symbols,FU_protons))

    Ihost_protons = [protonsOFspecie[ion] for ion in Ihost_symbols]

    return Ihost_protons
end

function get_charges(symbols, FU_symbols, FU_charges, Ihost_symbols)

    # Function that gets the charge of the species that appear in the initial structure (anchored
    # and host ions).
    #
    # Arguments:
    # symbols -------> Vector containing the chemical symbols of the ions in the initial structure.
    # FU_symbols ----> Vector containing the chemical symbols of the species that appear in the
    #                  initial structure (anchored and host ions).
    # FU_charges ----> Vector containing the charges of the species that appear in the initial
    #                  structure (anchored and host ions).
    # Ihost_symbols -> Vector containing the chemical symbols of the host ions.
    #
    # Returns:
    # charges -------> Vector containing the charges of the ions in the initial structure. It
    #                  indicates the charge of the ion associated to each lattice site in the
    #                  initial structure.
    # Ihost_charges -> Vector containing the charges of the host ions.

    chargeOFspecie = Dict(zip(FU_symbols,FU_charges))

    charges = zeros(length(symbols))
    for i in eachindex(symbols)
        charges[i] = chargeOFspecie[symbols[i]]
    end

    Ihost_charges = [chargeOFspecie[ion] for ion in Ihost_symbols]

    return charges, Ihost_charges
end

function check_QConservation(FU_symbols, FU_numbers, FU_charges, Ihost_symbols,
                             Ihost_numbers, Ihost_charges, Iadd_numbers, Iadd_charges)

    # Function that checks if the charge is preserved after removing a fraction of the host ions and
    # adding the dopants.
    #
    # Arguments:
    # FU_symbols ----> Vector containing the chemical symbols of the species that appear in the
    #                  initial structure (anchored and host ions).
    # FU_numbers ----> Vector containing the number of ions per formula unit of the species that
    #                  appear in the initial structure (anchored and host ions).
    # FU_charges ----> Vector containing the charges of the species that appear in the initial
    #                  structure (anchored and host ions).
    # Ihost_symbols -> Vector containing the chemical symbols of the host ions.
    # Ihost_numbers -> Vector containing the number of host ions (of each type) per formula unit in
    #                  the final structure.
    # Ihost_charges -> Vector containing the charges of the host ions.
    # Iadd_numbers --> Vector containing the number of dopants (of each type) per formula unit in
    #                  the final structure.
    # Iadd_charges --> Vector containing the charges of the dopants (added ions).
    #
    # Returns:
    # No returns

    # compute difference of energy
    QInitial = 0.0
    for ion in Ihost_symbols
        index = findfirst(x->x==ion, FU_symbols)
        QInitial += FU_numbers[index] * FU_charges[index]
    end
    if isempty(Iadd_numbers)
        QFinal   = sum(Ihost_charges .* Ihost_numbers)
    else
        QFinal   = sum(Ihost_charges .* Ihost_numbers) + sum(Iadd_charges .* Iadd_numbers)
    end
    ΔQ = QFinal - QInitial

    # check conservation
    if isapprox(ΔQ, 0.0, atol = 1e-3)
        QConservation = true
    else
        QConservation = false
    end

    # If charge is not conserved, print a message and continue
    if QConservation == false
        println("WARNING: CHARGE NOT CONSERVED (ΔQ = $ΔQ e)")
    end

end

function get_Sites(symbols::Vector{String}, Ihost_symbols::Vector{String})

    # Function that gets the optimization sites (sites used to place interchangeable ions and/or
    # vacancies during the optimization).
    # The chemical symbols of the anchored/host ions that occupy each site in the initial structure
    # are read, and when the symbols match with the chemical symbols of a host ion type, the sites
    # are selected.
    # The selected sites are called optimization sites.
    # The number of optimization sites is equal to the number of host ions in the initial structure.
    # The optimization site type (subset) is defined according to the host ion type occupying that
    # site in the initial structure.
    #
    # Arguments:
    # symbols -------> Vector containing the chemical symbols of the ions in the initial structure.
    # Ihost_symbols -> Vector containing the chemical symbols of the host ions.
    #
    # Returns:
    # optSites ----> Vector containing the optimization sites.
    # SitesSubset -> Vector containing the optimization sites subset of each optimization site.
    # NoptSites ---> Total number of optimization sites.
    # NoptSites_i -> Vector containing the number of optimization sites of each optimization sites
    #                subset.

    optSites  = Int64[]
    SitesSubset  = Int64[]
    NoptSites_i = Int64[]
    for (i, symbol) in enumerate(Ihost_symbols)
        symbol_sites = findall(x->(x==symbol), symbols)
        type_sites = [i for j in eachindex(symbol_sites)]
        append!(optSites, symbol_sites)
        append!(SitesSubset, type_sites)
        append!(NoptSites_i, length(symbol_sites))
    end
    NoptSites = length(optSites)

    return optSites, SitesSubset , NoptSites, NoptSites_i
end

function get_IonTypes(Iadd_symbols, Iadd_sites , Iadd_numbers, Iadd_charges, Iadd_protons,
                      Ihost_symbols, Ihost_numbers, Ihost_charges, Ihost_protons, FU_copies,
                      IonSiteConstraint::Bool)

    # Function that gets information related to the interchangeable ions.
    # Note:
    # When an interchangeable ion is in the initial structure (is partially removed), we call it
    # host ion.
    # When an interchangeable ion is not in the initial structure (is added), we call it
    # added ion / dopant.
    #
    # Arguments:
    # Iadd_symbols ------> Vector containing the chemical symbols of the dopants (added ions).
    # Iadd_sites --------> Vector containing the optimization sites subset of the dopants (e.i.,the
    #                      optimization sites subset that can be occupied by each added ion type).
    # Iadd_numbers ------> Vector containing the number of dopants (of each type) per formula unit
    #                      in the final structure.
    # Iadd_charges ------> Vector containing the charges of the dopants (added ions).
    # Iadd_protons ------> Vector containing the atomic numbers of the dopants (added ions).
    # Ihost_symbols -----> Vector containing the chemical symbols of the host ions.
    # Ihost_numbers -----> Vector containing the number of host ions (of each type) per formula unit
    #                      in the final structure.
    # Ihost_charges -----> Vector containing the charges of the host ions.
    # Ihost_protons -----> Vector containing the atomic numbers of the host ions.
    # FU_copies ---------> Number of unit formulas contained in the supercell.
    # IonSiteConstraint -> Turn on / off the ion-site constraint functionality.
    #
    # Returns:
    # Ions -------> Vector containing the reference ID of each interchangeable ion type.
    # TIons ------> Vector containing the status (host/added) of each interchangeable ion type.
    # IonsSubset -> Vector containing the optimization sites subset of each interchangeable ion type
    #               (the subset that can be occupied by that ion type).
    # SIons ------> Vector containing the chemical symbol of each interchangeable ion type.
    # LIons ------> Vector containing the number of ions of each interchangeable ion type.
    # QIons ------> Vector containing the charge of each interchangeable ion type.
    # ZIons ------> Vector containing the atomic number of each interchangeable ion type.

    # Initialize vectors
    SIons = String[Iadd_symbols; Ihost_symbols]
    aux_LIons = [Iadd_numbers; Ihost_numbers] * FU_copies
    LIons = Int64[round(Int64, LIon) for LIon in aux_LIons]
    QIons = Float64[Iadd_charges; Ihost_charges]
    ZIons = Int64[Iadd_protons; Ihost_protons]
    Ions  = collect(1:length(SIons))

    # Type of ion -> we use the notation host ion (3) for ions that are present in the initial
    #                structure, and added ion/dopant (4) for ions that are not present in the
    #                initial structure
    TIons = zeros(Int64, length(Ions))
    for i in Ions
        if i <= length(Iadd_symbols)
            TIons[i] = 4 # dopant
        else
            TIons[i] = 3 # host ion
        end
    end

    # Ion-site constraint:
    # The host ions of a given type can only occupy the optimization sites of the subset generated
    # by them (remember that each subset corresponds to a host ion type).
    # The dopants of a given type can only occupy the optimization sites of the subset defined in
    # the "settings.jl" file (Iadd_sites tag).
    # If ion-site constraint is not activated in the "settings.jl" file (IonSiteConstraint = false),
    # then these vectors are ignored and every combination is possible.
    if isempty(Iadd_sites)
        IonsSubset = Int64[]
    else
        IonsSubset = copy(Iadd_sites)
    end
    for i in eachindex(Ihost_symbols)
        append!(IonsSubset, i)
    end

    # if IonSiteConstraint = false, all optimization sites belong to subset 1.
    if IonSiteConstraint == false
        for i in eachindex(IonsSubset)
            IonsSubset[i] = 1
        end
    end

    return Ions, TIons, IonsSubset, SIons, LIons, QIons, ZIons

end

function get_data(settings::NamedTuple)

    # Function that reads the settings tuple and return the data tuple.
    # The settings tuple contains all the information that the user provide, and it is used to
    # read the structure of the material and generate the information that the optimizers need to
    # solve the problem.
    # The data tuple contains the variables that the optimizer needs.
    # 
    # Arguments:
    # settings -> Tuple containing the settings information (explained in the settings module).
    #
    # Returns:
    # data -----> Tuple containing the problem instance (information needed for the solver to
    #             initialize and solve the optimization problem). The description of the constants
    #             stored in the data tuple can be found below.

    println(" ")
    println("-------------------------------------------------------------------------------------------------")
    println("   START IMPORT DATA")
    println("-------------------------------------------------------------------------------------------------")

    @unpack IN_format, Ihost_symbols, Ihost_numbers, Iadd_symbols, Iadd_sites, Iadd_numbers,
            Iadd_charges, Iadd_protons, FU_symbols, FU_numbers, FU_charges, FU_protons, FU_copies,
            LR_type, LR_params, SR_type, SR_params, IonSiteConstraint = settings

    # (1) Read input atomic structure file (ASE object)
    structure = read_structure(IN_format)

    # (2) Get structure information
    # Lattice vectors
    lattice = pyconvert(Matrix{Float64}, structure.cell.array)
    # Positions of lattice sites
    positions = pyconvert(Matrix{Float64}, structure.get_positions())
    # Chemical symbols of the ions in the initial structure (chemical specie associated to each
    # lattice site)
    symbols = pyconvert(Vector{String}, structure.get_chemical_symbols())
    # Chemical symbols of the species that appear in the initial structure
    species = get_species(symbols)
    # Atomic numbers of the ions in the initial structure (atomic number associated to each
    # lattice site)
    numbers = pyconvert(Vector{Int}, structure.get_atomic_numbers())
    # Atomic numbers of host ions
    Ihost_protons = get_protons(FU_symbols, FU_protons, Ihost_symbols)
    # Number of lattice sites
    Ns = length(symbols)
    println("STRUCTURE IMPORTED")

    # Charges of the ions in the initial structure (charge associated to each lattice site) and
    # charge of each host ion type
    charges, Ihost_charges = get_charges(symbols, FU_symbols, FU_charges, Ihost_symbols)
    println("CHARGES DEFINED")
    # Check charge conservation
    check_QConservation(FU_symbols, FU_numbers, FU_charges, Ihost_symbols,
                        Ihost_numbers, Ihost_charges, Iadd_numbers, Iadd_charges)

    # (3) Identify optimization sites and define interchangeable ion types
    # optimization sites, subset of the optimization sites, total number of optimization sites,
    # number of optimization sites of each subset
    optSites, SitesSubset, NoptSites, NoptSites_i = get_Sites(symbols, Ihost_symbols)

    # (interchangeable) ion types: ion type ID, status (host ion / added ion),
    #                   available optimization sites subset, chemical symbol, number of ions,
    #                   charge, atomic number)
    Ions,TIons, IonsSubset,
    SIons, LIons,
    QIons, ZIons = get_IonTypes(Iadd_symbols, Iadd_sites , Iadd_numbers, Iadd_charges, Iadd_protons,
                                Ihost_symbols, Ihost_numbers, Ihost_charges, Ihost_protons,
                                FU_copies, IonSiteConstraint)
    # Total number of interchangeable ions
    L  = sum(Int64, LIons)
    # Number of vacancies (total optimization sites - total interchangeable ions)
    Nv = NoptSites - L 
    # Total number of atoms/ions (anchored + interchangeable atoms/ions = sites - vacancies) in the
    # final structure
    Na = Ns - Nv
    print("OPTIMIZATION SITES DEFINED: ")
    for i in eachindex(Ihost_symbols)
        print("$(Ihost_symbols[i])($(NoptSites_i[i])) ")
    end
    println(" ")
    print("INTERCHANGEABLE IONS DEFINED: ")
    for i in Ions
        print("$(SIons[i])($(LIons[i])) ")
    end
    println(" ")
    if IonSiteConstraint == true 
        print("ION-SITE CONSTRAINTS DEFINED: ")
        for i in eachindex(SIons)
            print("$(SIons[i])(Type $(IonsSubset[i])) ")
        end
        println(" ")
    end

    # (4) Precompute energy components
    E_AA, E_IA, E_II = Energy.get_PrecomputedEnergy(
        structure,
        Ns::Int,
        charges,
        optSites,
        NoptSites,
        Ions,
        TIons,
        QIons,
        LR_type,
        LR_params,
        SR_type,
        SR_params
    )

    # (5) Define data tuple
    data = 
    (
        Ns=Ns,
        NoptSites=NoptSites,
        Nv=Nv,
        Na=Na,
        IonSiteConstraint = IonSiteConstraint,
        Ions=Ions,
        LIons=LIons,
        IonsSubset=IonsSubset,
        SIons=SIons,
        QIons=QIons,
        ZIons=ZIons,
        lattice=lattice,
        positions=positions,
        symbols=symbols,
        species=species,
        optSites=optSites,
        SitesSubset=SitesSubset,
        U=E_IA,
        UionAionBlist=E_II,
        EnergyBase=E_AA
    )

    println("-------------------------------------------------------------------------------------------------")
    println("   END IMPORT DATA")
    println("-------------------------------------------------------------------------------------------------")

    return data
end

end