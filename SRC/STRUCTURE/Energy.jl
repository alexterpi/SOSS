module Energy

# This module contains the functions used to compute the interaction energy between two ions.
# There are 3 types of ions: anchored ion, host ion, added ion (dopant).

function energyInteraction(
    r12::Float64,
    q1::Float64, q2::Float64,
    z1::Int64, z2::Int64,
    type1::Int64, type2::Int64,
    r_overlap::Float64
    )

    # This function computes the interaction energy between two ions:
    # E = E_coulomb + E_overlap
    #
    # Arguments: 
    # r12 -----------> Distance between ions 1 and 2 (Å)
    # q1, q2 --------> Charge of ions 1 and 2 (e)
    # z1, z2 --------> Atomic number of ions 1 and 2
    # r_overlap -----> Overlap radius (Å)
    #
    # Returns:
    # e_total -> Total interaction energy between ions 1 and 2 (eV)

    # Coulomb interaction energy
    e_coulomb = coulomb(r12, q1, q2)

    if r_overlap < 0
        e_overlap = 0.0
    else
        e_overlap = overlap(r12, type1, type2, r_overlap)
    end

    # Total interaction energy
    e_total = e_coulomb + e_overlap

    return e_total
end

function coulomb(r12::Float64, q1::Float64, q2::Float64)

    # This function computes the Coulomb interaction energy between two ions:
    # E = C * (q1 * q2 / r12)
    #
    # Arguments:
    # q1, q2 -> Charge of ions 1 and 2 (e)
    # r12 ----> Distance between ions 1 and 2 (Å)
    #
    # Returns:
    # e_coulomb -> Coulomb interaction energy between ions 1 and 2 (eV)

    # C = 1/(4πϵ₀) is a conversion constant to enter charges in atomic units (e), distances in Å 
    # and get energy in eV
    C = 14.399645353504035

    # Coulomb interaction energy
    e_coulomb = C * q1 * q2 / r12

    return e_coulomb
    
end

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

end 