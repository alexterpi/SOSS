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
using PyCall
using UnPack
using Random

# Import energy module
include("Energy.jl")
using .Energy

# Import python function two swap from ASE format to Pymatgen format
adapt = pyimport("pymatgen.io.ase")


# Define Python functions
function __init__()
    py"""
    import numpy as np
    from ase.io.vasp import read_vasp

    def read_structure_py(format):

        # Function that reads the structure file and return the structure object.
        #
        # Arguments:
        # format -> String containing the format of the structure file.
        #
        # Returns:
        # structure -> ASE structure object.

        if format == "VASP":
            structure = read_vasp("instruct.vasp")
        else:
            print("Format not defined")
        return structure
    """
end

# Transform Python functions to Julia functions
read_structure(format) = py"read_structure_py"(format)

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

function get_cutoff(structure, CutOff)

    # Function that gets the cutoff distance of the interaction energy (neighbors are only
    # considered up to that distance).
    #
    # Arguments:
    # structure -> Pymatgen structure object.
    # CutOff ----> Cutoff distance .
    #              If CutOff < 0 -> The function assumes that the entered value is a fraction of the
    #                               max cutoff.
    #              If CutOff > 0 -> The function assumes that the entered value is the cutoff in Å.
    #
    # Returns:
    # CutOff -> Cutoff distance in Å.
    #
    # Note: If the entered cutoff is larger than the max cutoff, the value is accepted, but a
    #       warning is printed by screen.
    #       Since the code works with PBC, a cutoff value larger than the maximum cutoff value may
    #       cause sites to interact not only with the closest image of their neighbors, but also
    #       with the following images. 

    MaxCutOff = 0.5 * minimum(structure.lattice.abc)

    if CutOff < 0
        CutOff = abs(CutOff) * MaxCutOff
    end
    if CutOff > MaxCutOff
        println("WARNING: CUTOFF LARGER THAN THE EXPECTED MAXIMUM CUTOFF (RMAX = $MaxCutOff Å)")
    end
    return CutOff
end

function getDistNeighbors(structure, CutOff::Float64)

    # Function that gets the neighbors of each site (inside the cutoff distance) and the distances.
    #
    # Arguments:
    # structure -> Pymatgen structure object.
    # CutOff ----> Cutoff distance in Å.
    #
    # Returns
    # neighbors_of -> Vector of vectors, where each subvector contains the neighbors of each site.
    # distSite2Col -> Vector of vectors, where each subvector contains the distance of each site to
    #                 its neighbors.

    # Get the sites, neighbors and distances
    # (https://pymatgen.org/pymatgen.core.html#pymatgen.core.structure.IStructure.get_neighbor_list)
    sites, neighbors, _, distances = structure.get_neighbor_list(CutOff)
    # Shift the indexation (julia = python + 1)
    sites .+= 1
    neighbors .+= 1

    # Generate a vector with the neighbors of each site and a vector with the distances
    neighbors_of = Vector{Int64}[]
    distSite2Col = Vector{Float64}[]
    for site in 1:structure.num_sites
        cols, _ = whereInt(sites, site)
        push!(neighbors_of, neighbors[cols])
        push!(distSite2Col, distances[cols])
    end
    return distSite2Col, neighbors_of
end

function get_AnchAnch_energy(
    optSites::Array{Int,1},
    charges,
    numbers,
    distSite2Col,
    neighbors_of,
    r_overlap)

    # Function that obtains the total interaction energy of the structure without interchangeable
    # ions. That is, this is the total interaction energy of the anchored ions (ions that do not
    # move during the optimization).
    #
    # Arguments:
    # optSites -----> Vector containing the optimization sites. 
    # charges ------> Vector containing the charges of the ions in the initial structure.
    # numbers ------> Vector containing the atomic number of the ions in the initial structure.
    # distSite2Col -> Vector of vectors, where each subvector contains the distance of each site to
    #                 its neighbors.
    # neighbors_of -> Vector of vectors, where each subvector contains the neighbors of each site.
    #
    # Returns:
    # e -> Total anchored-anchored ion pair interaction energy in eV.

    e = 0
    for site in eachindex(charges)
        if site ∉ optSites
            for (col, neighborSite) in enumerate(neighbors_of[site])
                if neighborSite ∉ optSites
                    distance = distSite2Col[site][col]
                    QA = charges[site]
                    QB = charges[neighborSite]
                    ZA = numbers[site]
                    ZB = numbers[neighborSite]
                    e += Energy.energyInteraction(distance, QA, QB, ZA, ZB, 1, 1, r_overlap)
                end
            end
        end
    end
    e *= 0.5
    return e
end

function get_InterAnch_energy(
    optSites,
    neighbors_of,
    distSite2Col,
    charges,
    numbers,
    NoptSites,
    Ions,
    QIons,
    ZIons,
    TIons,
    r_overlap
    )

    # Function that gets the interchangeable-anchored ion pair interaction energy. The matrix
    # U[site s, type ion i] stores the interaction energy between an interchangeable ion of type i
    # in optimization site s with its anchored ion neighbors.
    #
    # Arguments:
    # optSites -----> Vector containing the optimization sites.
    # neighbors_of -> Vector of vectors, where each subvector contains the neighbors of each site.
    # distSite2Col -> Vector of vectors, where each subvector contains the distance of each site to
    #                 its neighbors.
    # charges ------> Vector containing the charges of the ions in the initial structure.
    # numbers ------> Vector containing the atomic number of the ions in the initial structure.
    # NoptSites ----> Total number of optimization sites.
    # Ions ---------> Vector containing the reference ID of each interchangeable ion type.
    # QIons --------> Vector containing the charge of each interchangeable ion type.
    # ZIons --------> Vector containing the atomic number of each interchangeable ion type.
    # TIons --------> Vector containing the status (host/added) of each interchangeable ion type.
    # r_overlap ----> Overlap radius (min possible distance between two interchangeable ions or and
    #                 interchangeable ion and an anchored ion)
    #
    # Returns:
    # U -> Precomputed components of the interchangeable-anchored ion pair interaction energy (in
    #      eV).

    U = zeros(NoptSites, Ions[end])

    for (i, site) in enumerate(optSites) # for all optimization sites
        for (col, neighborSite) in enumerate(neighbors_of[site]) # for all neighbors
            if neighborSite ∉ optSites # only accept if the neighbor contains an anchored ion
                distance = distSite2Col[site][col] # site-neighbor distance
                Q_B = charges[neighborSite] # neighbor charge
                Z_B = numbers[neighborSite] # neighbor atomic number
                for j in Ions # for every possible type of interchangeable ion
                    U[i, j] += Energy.energyInteraction(distance, QIons[j], Q_B, ZIons[j], Z_B,
                                                        TIons[j], 1, r_overlap)
                end
            end
        end
    end

    return U

end

function get_InterInter_energy(
    optSites,
    neighbors_of,
    distSite2Col,
    Ions,
    TIons,
    ZIons,
    QIons,
    r_overlap
    )

    # Function that gets the interchangeable-interchangeable ion pairs interaction energy. The
    # array UionAionBlist[ion i1, ion i2, site s1, site s2] stores the interaction energy between
    # interchangeable ion of type i1 in site s1 and interchangeable ion of type i2 in site s2.
    #
    # Arguments:
    # optSites -----> Vector containing the optimization sites.
    # neighbors_of -> Vector of vectors, where each subvector contains the neighbors of each site.
    # distSite2Col -> Vector of vectors, where each subvector contains the distance of each site to
    #                 its neighbors.
    # Ions ---------> Vector containing the reference ID of each interchangeable ion type.
    # TIons --------> Vector containing the status (host/added) of each interchangeable ion type.
    # ZIons --------> Vector containing the atomic number of each interchangeable ion type.
    # QIons --------> Vector containing the charge of each interchangeable ion type.
    # r_overlap ----> Overlap radius (min possible distance between two interchangeable ions or and
    #                 interchangeable ion and an anchored ion)
    #
    # Returns:
    # UionAionBlist -> Precomputed components of the interchangeable-interchangeable ion pair
    #                  interaction energy (in eV).

    # Get the optimization sites neighbors dictionary
    whereDic = tableWhere(optSites, neighbors_of)
    # Get the energy interaction dictionary
    pairInteracDic = tablePairInterac(optSites, distSite2Col, whereDic, QIons,
                                      ZIons, TIons, Ions, r_overlap)
    # Get the energy interaction array (summing the images if necessary)
    UionAionBlist  = tableUionAionB(optSites, whereDic, pairInteracDic, Ions)

    return UionAionBlist

end

function tableWhere(optSites, neighbors_of)

    # Function that gets the index of each optimization site in the neighbourhood of the rest of
    # optimization sites.
    #
    # Arguments:
    # optSites -----> Vector containing the optimization sites.
    # neighbors_of -> Vector of vectors, where each subvector contains the neighbors of each site.
    #
    # Returns:
    # dic -> Dictionary containing for each optimization site (first element) the index within its
    #        neighbourhood of another optimization site (second element).
    #        Example: If site 3 is the 2nd element (neighbor) in the neighbourhood of site 6, then
    #                 dic[6,3] = 2
    #        Note: If a site and its images fall inside the cutoff radius area of another site, then
    #              the dictionary will return more than one output for the same inputs. This can be
    #              avoided selecting a cutoff radius smaller or equal than a half of the length of
    #              the smaller side of the unit cell.

    dic = Dict()
    for A in optSites
        neighborsOfA = neighbors_of[A]
        for B in optSites
            columns, _   = whereInt(neighborsOfA, B)
            dic[A, B] = columns
        end
    end

    return dic
end

function tablePairInterac(optSites::Array{Int,1}, distSite2Col, whereDic::Dict{Any,Any},
                          QIons::Vector{Float64}, ZIons::Vector{Int64}, TIons::Vector{Int64},
                          Ions::Vector{Int64}, r_overlap::Float64)
    
    # Function that generates a dictionary containing the interaction energy between every possible
    # combination of two optimization sites (which are neighbors) and two interchangeable ion types.
    # 
    # Arguments:
    # optSites -----> Vector containing the optimization sites.
    # distSite2Col -> Vector of vectors, where each subvector contains the distance of each site to
    #                 its neighbors.
    # whereDic -----> Dictionary containing for each optimization site (first element) the index
    #                 within its neighbourhood of another site (second element).
    # QIons --------> Vector containing the charge of each interchangeable ion type.
    # ZIons --------> Vector containing the atomic number of each interchangeable ion type.
    # TIons --------> Vector containing the status (host/added) of each interchangeable ion type.
    # Ions ---------> Vector containing the reference ID of each interchangeable ion type.
    #
    # Returns:
    # dic -> Dictionary containing the interaction energy for each pair of sites and each possible
    #        combination of interchangeable ions in both sites.
    #        Terms: dic[site A, site B, image, ion type i, ion type j].
    #        Note: In principle, two sites should interact only once, but if the cutoff radius is
    #              larger than the maximum cutoff, a site can interact with another site and its
    #              images (due to PBC). For this reason, there is a term "image" in the dicctionary.
    #              If the cutoff is smaller or equal than the maximum cutoff, then there will
    #              be necessarily only one image, and the third term of the dictionary will become
    #              irrelevant.

    dic = Dict()
    for A in optSites
        dist_A_To_col = distSite2Col[A]
        for B in optSites
            for colB in whereDic[A, B]
                for i in eachindex(Ions)
                    for j in eachindex(Ions)
                        e = Energy.energyInteraction(dist_A_To_col[colB], QIons[i], QIons[j],
                                                     ZIons[i],  ZIons[j], TIons[i], TIons[j],
                                                     r_overlap)
                        dic[A, B, colB, Ions[i], Ions[j]] = e
                    end
                end
            end
        end
    end
    return dic
end

function tableUionAionB(optSites::Array{Int,1}, whereDic::Dict{Any,Any},
                        pairInteracDic::Dict{Any,Any}, Ions::Vector{Int64})

    # Function that reads the dictionary generated with the function tablePairInterac(...) and
    # writes it in form of 4D array.
    # Note: If two sites interact more than once (interaction with the site and its images), 
    #       all images are added together (summed) in the same interaction.
    #
    # Arguments:
    # optSites -------> Vector containing the optimization sites.
    # whereDic -------> Dictionary containing for each site (first element) the index within its
    #                   neighbourhood of another site (second element).
    # pairInteracDic -> Dictionary containing the interaction energy for each pair of sites and each
    #                   possible combination of interchangeable ions in both sites. 
    #                   Dic[site A, site B, image, type ion A, tipe ion B].
    # Ions -----------> Vector containing the reference ID of each interchangeable ion type.
    #
    # Returns:
    # array -> Array containing the interaction energy for each pair of optimization sites and each
    #          possible combination of interchangeable ions in both sites.
    #          Terms: array[site A, site B, ion type i, ion type j].

    n = length(optSites)
    m = length(Ions)
    array = zeros( n, n, m, m )
    #
    for l in Ions
        for k in Ions
            for (j, B) in enumerate(optSites)
                for (i, A) in enumerate(optSites)
                    for colB in whereDic[A, B]
                        array[i,j,k,l] += pairInteracDic[ A, B, colB, k, l ]
                    end
                end
            end
        end
    end
    return array
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
            r_cutoff, r_overlap, IonSiteConstraint = settings

    # (1) Read input atomic structure file (ASE object)
    structure = read_structure(IN_format)

    # (2) Get structure information
    # Lattice vectors
    lattice = structure.cell.array
    # Positions of lattice sites
    positions = structure.get_positions()
    # Chemical symbols of the ions in the initial structure (chemical specie associated to each
    # lattice site)
    symbols = structure.get_chemical_symbols()
    # Chemical symbols of the species that appear in the initial structure
    species = get_species(symbols)
    # Atomic numbers of the ions in the initial structure (atomic number associated to each
    # lattice site)
    numbers = structure.get_atomic_numbers()
    # Atomic numbers of host ions
    Ihost_protons = get_protons(FU_symbols, FU_protons, Ihost_symbols)
    # Number of lattice sites
    Ns = length(symbols)
    # show by screen
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
    # show by screen
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

    # (4) Identify neighbors
    # Change to Pymatgen object (we use a pymatgen function to compute neighbors)
    structure = adapt.AseAtomsAdaptor().get_structure(structure)
    # Get the cutoff radius and neighbors (using PBC)
    r_cutoff = get_cutoff(structure, r_cutoff)
    distSite2Col, neighbors_of = getDistNeighbors(structure, r_cutoff)
    println("NEIGHBORS DEFINED")

    # (5) Precompute components of the interaction energy (objective function)
    # anchored-anchored ion pairs energy
    EnergyBase = get_AnchAnch_energy(
        optSites,
        charges,
        numbers,
        distSite2Col,
        neighbors_of,
        r_overlap
    )
    println("ANCHORED-ANCHORED ION PAIRS INTERACTION ENERGY DEFINED")
    # interchangeable-anchored ion pairs energy
    U = get_InterAnch_energy(
        optSites,
        neighbors_of,
        distSite2Col,
        charges,
        numbers,
        NoptSites,
        Ions,
        QIons,
        ZIons,
        TIons,
        r_overlap
    )
    println("INTERCHANGEABLE-ANCHORED ION PAIRS INTERACTION ENERGY DEFINED")
    # interchangeable-interchangeable ion pairs energy
    UionAionBlist = get_InterInter_energy(
        optSites,
        neighbors_of,
        distSite2Col,
        Ions,
        TIons,
        ZIons,
        QIons,
        r_overlap
    )
    println("INTERCHANGEABLE-INTERCHANGEABLE ION PAIRS INTERACTION ENERGY DEFINED")

    # (6) Define data tuple
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
        U=U,
        UionAionBlist=UionAionBlist,
        EnergyBase=EnergyBase
    ) 

    println("-------------------------------------------------------------------------------------------------")
    println("   END IMPORT DATA")
    println("-------------------------------------------------------------------------------------------------")

    return data
end

end