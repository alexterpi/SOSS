module Energy

# import Julia modules

using PythonCall
using LinearAlgebra
using Printf

# define references to Python modules

const pmg_core = Ref{Py}(Py(nothing))
const pmg_ewald = Ref{Py}(Py(nothing))
const pymatgen_ase = Ref{Py}(Py(nothing))

# initialize Python environment

function __init__()
    pmg_core[] = pyimport("pymatgen.core")
    try
        pmg_ewald[] = pyimport("pymatgen.core.ewald")
    catch e
        pmg_ewald[] = pyimport("pymatgen.analysis.ewald")
    end
    pymatgen_ase[] = pyimport("pymatgen.io.ase")
end

# define main function

function get_PrecomputedEnergy(
    structure,
    Ns::Int,
    charges::Array{Float64, 1},
    optSites::Array{Int, 1},
    NoptSites::Int,
    Ions::Array{Int, 1},
    TIons::Array{Int, 1},
    QIons::Array{Float64, 1},
    LR_type::Int,
    LR_params::Array{Float64, 1},
    SR_type::Int,
    SR_params::Array{Float64, 1}
    )

    if LR_type == 1
        E_AA, E_IA, E_II = get_PrecomputedEnergy_Coulomb(structure,
                                                         charges,
                                                         optSites,
                                                         NoptSites,
                                                         Ions,
                                                         TIons,
                                                         QIons,
                                                         LR_params[1],
                                                         SR_type,
                                                         SR_params)
    elseif LR_type == 2
        E_AA, E_IA, E_II = get_PrecomputedEnergy_Ewald(structure,
                                                       Ns,
                                                       charges,
                                                       optSites,
                                                       NoptSites,
                                                       TIons,
                                                       QIons,
                                                       SR_type,
                                                       SR_params)
    else
        println("ERROR: INCORRECT LONG RANGE INTERACTION. CHECK THE SETTINGS FILE.")
        E_AA, E_IA, E_II = 0, 0, 0 
    end

    return E_AA, E_IA, E_II

end

# define long-range interaction functions

function get_PrecomputedEnergy_Coulomb(
    structure,
    charges::Array{Float64, 1},
    optSites::Array{Int, 1},
    NoptSites::Int,
    Ions::Array{Int, 1},
    TIons::Array{Int, 1},
    QIons::Array{Float64, 1},
    r_cutoff::Float64,
    SR_type::Int,
    SR_params::Array{Float64, 1}
    )

    # C = 1/(4πϵ₀) is a conversion constant to enter charges in atomic units (e), distances in Å 
    # and get energy in eV
    C = 14.399645353504035

    # Change to Pymatgen structure object (we use Pymatgen to compute energies)
    structure = pymatgen_ase[].AseAtomsAdaptor().get_structure(structure)

    # check cutoff radius
    r_cutoff = get_cutoff(structure, r_cutoff)

    # get neighbors and distances
    distSite2Col, neighbors_of = getDistNeighbors(structure, r_cutoff)
    println("NEIGHBORS DEFINED")

    # precompute components of the interaction energy

    # E_AA: anchored-anchored ion pairs energy
    E_AA = 0
    for site in eachindex(charges)
        if site ∉ optSites
            for (col, neighborSite) in enumerate(neighbors_of[site])
                if neighborSite ∉ optSites
                    distance = distSite2Col[site][col]
                    QA = charges[site]
                    QB = charges[neighborSite]
                    # Coulomb (long-range) term
                    E_AA += C * QA * QB / distance
                    # short-range term (optional)
                    if SR_type == 1 # overlap function
                        E_AA += overlap(distance, 1, 1, SR_params[1])
                    end
                end
            end
        end
    end
    E_AA *= 0.5
    println("ANCHORED-ANCHORED ION PAIRS INTERACTION ENERGY DEFINED")

    # E_IA: interchangeable-anchored ion pairs energy
    E_IA = zeros(NoptSites, Ions[end])
    for (i, site) in enumerate(optSites) # for all optimization sites
        for (col, neighborSite) in enumerate(neighbors_of[site]) # for all neighbors
            if neighborSite ∉ optSites # only accept if the neighbor contains an anchored ion
                distance = distSite2Col[site][col] # site-neighbor distance
                Q_B = charges[neighborSite] # neighbor charge
                for j in Ions # for every possible type of interchangeable ion
                    # Coulomb (long-range) term
                    E_IA[i, j] += C * QIons[j] * Q_B / distance
                    # short-range term (optional)
                    if SR_type == 1 # overlap function
                        E_IA[i, j] += overlap(distance, TIons[j], 1, SR_params[1])
                    end
                end
            end
        end
    end
    println("INTERCHANGEABLE-ANCHORED ION PAIRS INTERACTION ENERGY DEFINED")

    # E_II: interchangeable-interchangeable ion pairs energy
    # Get the optimization sites neighbors dictionary
    whereDic = Dict()
    for A in optSites
        neighborsOfA = neighbors_of[A]
        for B in optSites
            columns, _   = whereInt(neighborsOfA, B)
            whereDic[A, B] = columns
        end
    end
    # Get the energy interaction dictionary
    pairInteracDic = Dict()
    for A in optSites
        dist_A_To_col = distSite2Col[A]
        for B in optSites
            for colB in whereDic[A, B]
                for i in eachindex(Ions)
                    for j in eachindex(Ions)
                        # Coulomb (long-range) term
                        e = C * QIons[i] * QIons[j] / dist_A_To_col[colB]
                        # short-range term (optional)
                        if SR_type == 1 # overlap function
                            e += overlap(dist_A_To_col[colB], TIons[i], TIons[j], SR_params[1])
                        end
                        pairInteracDic[A, B, colB, Ions[i], Ions[j]] = e
                    end
                end
            end
        end
    end
    # Get the energy interaction array (summing the images if necessary)
    n = length(optSites)
    m = length(Ions)
    E_II = zeros( n, n, m, m )
    for l in Ions
        for k in Ions
            for (j, B) in enumerate(optSites)
                for (i, A) in enumerate(optSites)
                    for colB in whereDic[A, B]
                        E_II[i,j,k,l] += pairInteracDic[ A, B, colB, k, l ]
                    end
                end
            end
        end
    end
    println("INTERCHANGEABLE-INTERCHANGEABLE ION PAIRS INTERACTION ENERGY DEFINED")

    return E_AA, E_IA, E_II

end

function get_PrecomputedEnergy_Ewald(    
    structure,
    Ns::Int,
    charges::Array{Float64, 1},
    optSites::Array{Int, 1},
    NoptSites::Int,
    TIons::Array{Int, 1},
    QIons::Array{Float64, 1},
    SR_type::Int,
    SR_params::Array{Float64, 1}
    )

    # Change to Pymatgen structure object (we use Pymatgen to compute energies)
    structure = pymatgen_ase[].AseAtomsAdaptor().get_structure(structure)
    
    # get number of interchangeable ion types
    NIons = length(QIons)

    # get distances between sites (nearest image distance)
    distances = pyconvert(Matrix{Float64}, structure.distance_matrix)

    # set interchangeable ions charges to 0
    anchored_charges = copy(charges)
    for i in 1:NoptSites
        site = optSites[i]
        anchored_charges[site] = 0.0
    end
    
    println("Computing Ewald Vij matrix...")

    # for each site set charge to 1
    for i in 0:(Ns-1); structure[i].charge = 1.0; end
    # compute Ewald matrix prefactors
    ewald = pmg_ewald[].EwaldSummation(structure, acc_factor=12.0)
    V = pyconvert(Matrix{Float64}, ewald.total_energy_matrix) .* 2.0
    
    println("Computing E_AA...")
    E_AA = 0.0
    for i in 1:Ns
        qi = anchored_charges[i]
        if qi == 0.0 continue end # ignore interchangeable
        for j in 1:Ns
            qj = anchored_charges[j]
            if qj == 0.0 continue end # ignore interchangeable
            # Ewald (long-range) term
            E_AA += qi * qj * V[i, j]
            # short-range term (optional)
            if i ≠ j # no self interaction
                if SR_type == 1 # overlap function
                    E_AA += overlap(distances[i, j], 1, 1, SR_params[1])
                end
            end
        end
    end
    E_AA *= 0.5

    println("Computing E_IA...")
    E_IA = zeros(Float64, NoptSites, NIons)
    for i_idx in 1:NoptSites # optimization sites
        i_glob = optSites[i_idx] 
        for m in 1:NIons # interchangeable ion types
            qm = QIons[m]
            tm = TIons[m] 
            inter_with_anch = 0.0
            for j_glob in 1:Ns # all sites
                q_anch = anchored_charges[j_glob]
                if q_anch == 0.0 continue end # skip optimization sites (only anchored ions)
                # Ewald (long-range) term
                inter_with_anch += qm * q_anch * V[i_glob, j_glob]
                # short-range term (optional)
                if SR_type == 1 # overlap function
                    inter_with_anch += overlap(distances[i_glob, j_glob], tm, 1, SR_params[1])
                end
            end
            # Ewald (self energy)
            E_IA[i_idx, m] = inter_with_anch + 0.5 * qm^2 * V[i_glob, i_glob]
        end
    end

    println("Computing E_II...")        
    E_II = zeros(Float64, NoptSites, NoptSites, NIons, NIons)
    for i in 1:NoptSites # optimization sites
        i_glob = optSites[i]
        for j in 1:NoptSites # optimization sites
            if i ≠ j # no self interaction 
                j_glob = optSites[j]
                for m1 in 1:NIons # interchangeable ion types
                    q1 = QIons[m1]
                    t1 = TIons[m1]
                    for m2 in 1:NIons # interchangeable ion types
                        q2 = QIons[m2]
                        t2 = TIons[m2]
                        # Ewald (long-range) term
                        E_II[i, j, m1, m2] = q1 * q2 * V[i_glob, j_glob]
                        # short-range term (optional)
                        if SR_type == 1 # overlap function
                            E_II[i, j, m1, m2] += overlap(distances[i_glob, j_glob], t1, t2, SR_params[1])
                        end
                    end
                end
            end
        end
    end
    
    return E_AA, E_IA, E_II
end

# define short-range interaction functions

function overlap(r12::Float64, type1::Int64, type2::Int64, r_overlap::Float64)

    # This function computes the overlap interaction energy between two ions.
    # The overlap energy is indeed a constraint added via interaction energy to avoid some 
    # ions from being closer than a defined distance.
    # When two ions (excluding anchored-anchored pairs) are closer than r_overlap, a positive and
    # high energy h_overlap is returned. Otherwise, 0 is returned.
    #
    # Arguments:
    # r12 -------> Distance between ion 1 and ion 2 (Å)
    # type1 -----> Type of ion 1 (1 -> anchored ion, 3 -> host ion, 4 -> added ion)
    # type2 -----> Type of ion 2 (1 -> anchored ion, 3 -> host ion, 4 -> added ion)
    # r_overlap -> Overlap radius (Å)
    #
    # Returns:
    # overlap_energy -> Overlap interaction energy (eV)

    # Barrier height (eV)
    h_overlap = 100.0
    
    # Type of interaction:
    # anchored-anchored(2), anchored-host(4), anchored-added(5), host-host(6), host-added(7),
    # added-added(8)

    TypeInteraction = type1 + type2

    if (r12 <= r_overlap) && (TypeInteraction > 2)
        overlap_energy = h_overlap
    else
        overlap_energy = 0.0
    end

    return overlap_energy

end

# define auxiliar functions

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

function get_cutoff(structure, r_cutoff)

    # Function that gets the cutoff distance of the interaction energy (neighbors are only
    # considered up to that distance).
    #
    # Arguments:
    # structure -> Pymatgen structure object.
    # r_cutoff ----> Cutoff distance .
    #              If r_cutoff < 0 -> The function assumes that the entered value is a fraction of the
    #                               max cutoff.
    #              If r_cutoff > 0 -> The function assumes that the entered value is the cutoff in Å.
    #
    # Returns:
    # r_cutoff -> Cutoff distance in Å.
    #
    # Note: If the entered cutoff is larger than the max cutoff, the value is accepted, but a
    #       warning is printed by screen.
    #       Since the code works with PBC, a cutoff value larger than the maximum cutoff value may
    #       cause sites to interact not only with the closest image of their neighbors, but also
    #       with the following images. 

    r_max = 0.5 * minimum(pyconvert(Vector{Float64}, structure.lattice.abc))

    if r_cutoff < 0
        r_cutoff = abs(r_cutoff) * r_max
    end
    if r_cutoff > r_max
        println("WARNING: CUTOFF LARGER THAN THE EXPECTED MAXIMUM CUTOFF (RMAX = $r_max Å)")
    end
    return r_cutoff
end

function getDistNeighbors(structure, r_cutoff::Float64)

    # Function that gets the neighbors of each site (inside the cutoff distance) and the distances.
    #
    # Arguments:
    # structure -> Pymatgen structure object.
    # r_cutoff ----> Cutoff distance in Å.
    #
    # Returns
    # neighbors_of -> Vector of vectors, where each subvector contains the neighbors of each site.
    # distSite2Col -> Vector of vectors, where each subvector contains the distance of each site to
    #                 its neighbors.

    # Get the sites, neighbors and distances
    # (https://pymatgen.org/pymatgen.core.html#pymatgen.core.structure.IStructure.get_neighbor_list)

    # get neighbors
    res = structure.get_neighbor_list(r_cutoff)

    # convert to julia objects
    sites = pyconvert(Vector{Int64}, res[0]) 
    neighbors = pyconvert(Vector{Int64}, res[1])
    distances = pyconvert(Vector{Float64}, res[3])

    # Shift the indexation (julia = python + 1)
    sites .+= 1
    neighbors .+= 1

    # Generate a vector with the neighbors of each site and a vector with the distances
    neighbors_of = Vector{Int64}[]
    distSite2Col = Vector{Float64}[]
    num_sites = pyconvert(Int, structure.num_sites)
    for site in 1:num_sites
        cols, _ = whereInt(sites, site)
        push!(neighbors_of, neighbors[cols])
        push!(distSite2Col, distances[cols])
    end
    return distSite2Col, neighbors_of
end

end