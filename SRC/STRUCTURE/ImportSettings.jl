module ImportSettings

# This module contains the functions to read the settings file (settings.jl) and generate the
# settings tupple, which stores the settings tags defined in settigs.jl.
# The main function, get_settings(...), calls the functions get_settings_with_dopant(...) (when
# dopants are defined) or get_settings_without_dopant(...) (when dopants are not defined) to
# generate the settings tuple, then displays the settings tags on the screen, and, finally, returns
# the settings tuple.
# Check the settings files provided in the testing directories or the manual for a detailed
# description of the settings tags.

function get_settings_with_dopant()

    # Function that creates the settigs tuple when dopants are defined.
    #
    # Arguments:
    # No arguments.
    #
    # Returns:
    # settigs -> Settings tuple.

    if solver == "SANSS"
        settings =
        (
            solver = solver::String,
            params = "params.jl"::String,
            IN_format = IN_format::String,
            OUT_format = OUT_format::String,
            FU_symbols = FU_symbols::Vector{String},
            FU_numbers = (1.0 * FU_numbers)::Vector{Float64},
            FU_charges = (1.0 * FU_charges)::Vector{Float64},
            FU_protons = FU_protons::Vector{Int64},
            FU_copies = FU_copies::Int64,
            Ihost_symbols = Ihost_symbols::Vector{String},
            Ihost_numbers = (1.0 * Ihost_numbers)::Vector{Float64},
            Iadd_symbols = Iadd_symbols::Vector{String},
            Iadd_sites = Iadd_sites::Vector{Int64},
            IonSiteConstraint = IonSiteConstraint::Bool,
            Iadd_numbers = (1.0 * Iadd_numbers)::Vector{Float64},
            Iadd_charges = (1.0 * Iadd_charges)::Vector{Float64},
            Iadd_protons = Iadd_protons::Vector{Int64},
            r_cutoff = (1.0 * r_cutoff)::Float64,
            r_overlap = r_overlap::Float64
        )
    else
        settings =
        (
            solver = solver::String,
            params = "params.jl"::String,
            IN_format = IN_format::String,
            OUT_format = OUT_format::String,
            FU_symbols = FU_symbols::Vector{String},
            FU_numbers = (1.0 * FU_numbers)::Vector{Float64},
            FU_charges = (1.0 * FU_charges)::Vector{Float64},
            FU_protons = FU_protons::Vector{Int64},
            FU_copies = FU_copies::Int64,
            Ihost_symbols = Ihost_symbols::Vector{String},
            Ihost_numbers = (1.0 * Ihost_numbers)::Vector{Float64},
            Iadd_symbols = Iadd_symbols::Vector{String},
            Iadd_sites = Iadd_sites::Vector{Int64},
            IonSiteConstraint = false::Bool,
            Iadd_numbers = (1.0 * Iadd_numbers)::Vector{Float64},
            Iadd_charges = (1.0 * Iadd_charges)::Vector{Float64},
            Iadd_protons = Iadd_protons::Vector{Int64},
            r_cutoff = (1.0 * r_cutoff)::Float64,
            r_overlap = r_overlap::Float64
        )
    end

    return settings
end

function get_settings_without_dopant()

    # Function that creates the settigs tuple when dopants are not defined.
    #
    # Arguments:
    # No arguments.
    #
    # Returns:
    # settigs -> Settings tuple.

    if solver == "SANSS"
        settings =
        (
            solver = solver::String,
            params = "params.jl"::String,
            IN_format = IN_format::String,
            OUT_format = OUT_format::String,
            FU_symbols = FU_symbols::Vector{String},
            FU_numbers = (1.0 * FU_numbers)::Vector{Float64},
            FU_charges = (1.0 * FU_charges)::Vector{Float64},
            FU_protons = FU_protons::Vector{Int64},
            FU_copies = FU_copies::Int64,
            Ihost_symbols = Ihost_symbols::Vector{String},
            Ihost_numbers = (1.0 * Ihost_numbers)::Vector{Float64},
            Iadd_symbols = Iadd_symbols::Vector{Any},
            Iadd_sites = Iadd_sites::Vector{Any},
            IonSiteConstraint = IonSiteConstraint::Bool,
            Iadd_numbers = Iadd_numbers::Vector{Any},
            Iadd_charges = Iadd_charges::Vector{Any},
            Iadd_protons = Iadd_protons::Vector{Any},
            r_cutoff = (1.0 * r_cutoff)::Float64,
            r_overlap = r_overlap::Float64
        )
    else
        settings =
        (
            solver = solver::String,
            params = "params.jl"::String,
            IN_format = IN_format::String,
            OUT_format = OUT_format::String,
            FU_symbols = FU_symbols::Vector{String},
            FU_numbers = (1.0 * FU_numbers)::Vector{Float64},
            FU_charges = (1.0 * FU_charges)::Vector{Float64},
            FU_protons = FU_protons::Vector{Int64},
            FU_copies = FU_copies::Int64,
            Ihost_symbols = Ihost_symbols::Vector{String},
            Ihost_numbers = (1.0 * Ihost_numbers)::Vector{Float64},
            Iadd_symbols = Iadd_symbols::Vector{Any},
            Iadd_sites = Iadd_sites::Vector{Any},
            IonSiteConstraint = false::Bool,
            Iadd_numbers = Iadd_numbers::Vector{Any},
            Iadd_charges = Iadd_charges::Vector{Any},
            Iadd_protons = Iadd_protons::Vector{Any},
            r_cutoff = (1.0 * r_cutoff)::Float64,
            r_overlap = r_overlap::Float64
        )
    end

    return settings
end

function get_settings()

    # Function that gets the setting tuple by calling get_settings_with_dopant(...) (if dopants are
    # defined) or get_settings_without_dopant(...) (if no dopants are defined), displays the
    # settings tags on the screen, and returns the settigs tuple. 
    #
    # Arguments:
    # No arguments.
    #
    # Returns:
    # settigs -> Setings tuple.

    println(" ")
    println("-------------------------------------------------------------------------------------------------")
    println("   START IMPORT SETTINGS")
    println("-------------------------------------------------------------------------------------------------")
    try
        global settings = get_settings_with_dopant() # dopants
    catch e
        println("OPTIMIZATION WITHOUT DOPANTS")
        global settings = get_settings_without_dopant() # no dopants
    end
    # display settings tags
    for name in eachindex(settings)
        value = settings[name]
        println("$name = $value")
    end
    println("-------------------------------------------------------------------------------------------------")
    println("   END IMPORT SETTINGS")
    println("-------------------------------------------------------------------------------------------------")

    return settings
    
end

end