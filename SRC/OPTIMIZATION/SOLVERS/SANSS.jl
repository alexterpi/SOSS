module SANSS

# SANSS -> Simulated Anealing Structure Search

# This module contains the functions used to find the optimal configutaion (lower energy) of ions
# and vacancies.
# The main function SANSS_set(...), defined at the botton of the file, takes the data and
# params tuples as arguments (check ImportData and ImportSettings modules for more details), 
# initializes the configuration of ions and vacancies at random (appliying constraints if required)
# and optimizes the system using the Simulated Annealing (SA) algorithm, finally, a solution tuple,
# which contains the best solutions (configurations of ions) found and their energies, is returned.

# Reminder: anchored ions -> immobile ions (its position is not optimized)
#           interchangeable ions -> mobile ions (its position is optimized)
# Notation: "ion" (without specifying its status) -> "interchangeable ion".
#           for "anchored ions", the status is always specified.

# load modules
using StatsBase
using Statistics
using Random
using UnPack
using MutableNamedTuples
using PyPlot
using DelimitedFiles
using Printf

####################################################################################################
# functions devoted to writing/plotting data 

function plot_E(
    RunDir::String
    )

    # Function that reads the optimization data from the data files generated during the
    # optimization (acceptance_evol.dat, energy_evol.dat, best_solution.dat) and plots the data.
    # The three plots show:
    #   1. Evolution of the acceptance rate at each main cooling step
    #   2. Evolution of the best and last solution energy at each main cooling step
    #   3. Energy of the best solutions sampled with SOSS
    #
    # Arguments:
    # RunDir -> Directory containing the data files.
    #
    # Returns:
    # No returns, the function generate the plots in the data directory.

    ##########################################
    # get data
    data = readdlm("$RunDir/acceptance_evol.dat")
    acceptance = data[2:end,2]

    data = readdlm("$RunDir/energy_evol.dat")
    bestE = data[2:end,2]
    lastE = data[2:end,3]

    data = readdlm("$RunDir/best_solution.dat")
    sampled_bestE = data[2:end-1,2]
    ##########################################

    # acceptance vs step
    xval = collect(1:length(acceptance))
    xmin = 0.0
    xmax = xval[end] + 1
    ymax = 1.0
    ymin = 0.0
    PyPlot.ioff()
    PyPlot.plot(xval, acceptance, ".k")
    PyPlot.xlim(xmin, xmax)
    PyPlot.ylim(ymin, ymax)
    PyPlot.title("Evolution of Acceptance")
    PyPlot.xlabel("Step")
    PyPlot.ylabel("Acceptance")
    PyPlot.savefig("$RunDir/acceptance_evol.pdf")
    PyPlot.close()

    # energy vs step
    xval = collect(1:length(lastE))
    xmin = 0.0
    xmax = xval[end] + 1
    ymax = maximum([bestE; lastE; bestE]) + 0.02 * abs(maximum([bestE; lastE; bestE]))
    ymin = minimum([bestE; lastE; bestE]) - 0.02 * abs(minimum([bestE; lastE; bestE]))
    PyPlot.ioff()
    PyPlot.plot(xval, lastE, ".k", label = "Last Solution")
    PyPlot.plot(xval, bestE, "-k", label = "Best Solution")
    PyPlot.xlim(xmin, xmax)
    PyPlot.ylim(ymin, ymax)
    PyPlot.title("Evolution of Energy")
    PyPlot.xlabel("Step")
    PyPlot.ylabel("Energy (eV/atom)", color="k")
    PyPlot.legend(loc = "upper right", framealpha = 1, edgecolor = "k")
    PyPlot.savefig("$RunDir/energy_evol.pdf")
    PyPlot.close()

    # best solutions sampled with SOSS
    ΔE = sampled_bestE[end] - sampled_bestE[1]
    xval = collect(1:length(sampled_bestE))
    xmin = 0.0
    xmax = xval[end] + 1
    ymax = maximum(sampled_bestE) + 0.01 * abs(maximum(sampled_bestE))
    ymin = minimum(sampled_bestE) - 0.01 * abs(minimum(sampled_bestE))
    PyPlot.ioff()
    PyPlot.plot(xval, sampled_bestE, ".k", label="ΔE = $ΔE eV/atom")
    PyPlot.xlim(xmin, xmax)
    PyPlot.ylim(ymin, ymax)
    PyPlot.title("Energy of the best solutions sampled with SOSS")
    PyPlot.xlabel("Solution")
    PyPlot.ylabel("Energy (eV/atom)")
    PyPlot.legend(loc="upper right")
    PyPlot.savefig("$RunDir/best_solution.pdf")
    PyPlot.close()
end

function write_E(
    sampled_bestE::Vector{Float64},
    lastE::Vector{Float64},
    bestE::Vector{Float64},
    acceptance::Vector{Float64},
    RunDir::String
    )

    # Function that writes the acceptance evolution, energy evolution and energy of the current best
    # solution in the acceptance_evol.dat, energy_evol.dat, and best_solution.dat files.
    #
    # Arguments:
    # sampled_bestE -> Vector containing the energy of the best solutions sampled.
    # lastE ---------> Vector containing the energy of the last solution at each step. 
    # bestE ---------> Vector containing the energy of the best solution at each step.
    # acceptance ----> Vector containing the the acceptance rate at each step.
    # RunDir --------> Directory containing the data files.
    #
    # Returns:
    # No returns, the function generates the data files in the data directory.

    # acceptance evolution
    file = open("$RunDir/acceptance_evol.dat", "w")
    write(file, "# step, acceptance\n")
    for i in eachindex(acceptance)
        write(file, "$i $(acceptance[i])\n")
    end
    close(file)

    # energy evolution
    file = open("$RunDir/energy_evol.dat", "w")
    write(file, "# step, best E (eV/atom), last E (eV/atom)\n")
    for i in eachindex(lastE)
        write(file, "$i $(bestE[i]) $(lastE[i])\n")
    end
    close(file)

    # energy of the best solutions sampled with SOSS
    file = open("$RunDir/best_solution.dat", "w")
    write(file, "# configuration, E (eV/atom)\n")
    for i in eachindex(sampled_bestE)
        write(file, "$i $(sampled_bestE[i])\n")
    end
    ΔE = sampled_bestE[end] - sampled_bestE[1]
    write(file, "# ΔE = $ΔE eV/atom")
    close(file)
end

function write_opt(
    nRepeats::Int64,
    tempLength::Int64,
    scheme::String,
    IT_scheme::Int64,
    IT_params::Vector{Float64},
    steps::Int64,
    beta::Float64,
    gamma::Int64,
    FirstToImprove::Bool,
    ScanningStyle::Int64,
    MaxStagnation::Int64,
    accepted::Int64,
    EnergyList::Vector{Float64},
    RunList::Vector{Int64},
    time::Float64,
    StepsDir::String,
    DisplayStyle::Int64,
    plots::Bool
    )

    # Function that writes a summary of the optimization.
    # The summary contains the value of the SANSS parameters, the optimal energy for each run,
    # the average optimal energy, the best and worst optimal energy and the optimization time.
    #
    # Arguments:
    # nRepeats -------> Number of runs / optimization cycles.
    # tempLength -----> Number of steps in the secondary cooling scheme.
    # scheme ---------> Type of secondary cooling scheme.
    # IT_scheme ------> Scheme used to compute the initial temperature (Tmax and Tmin).
    # IT_params ------> Vector containing the parameters required by the used initial temperature
    #                   scheme.
    # steps ----------> Number of steps in the main cooling scheme.
    # beta -----------> Intensity of the secondary cooling scheme. T_final = beta * T_initial.
    # gamma ----------> Max number of combinations explored in each move.
    # FirstToImprove -> Skip or continue the exploration of moves when a successfull moves is
    #                   found.
    # ScanningStyle --> Method to decide how the number of combinations explored in each move
    #                   varies during the calculation.
    # MaxStagnation --> Limit of stagnation.
    # accepted -------> Number of accepted optimal configurations.
    # EnergyList -----> Vector containing the energy of the optimal configuration of each run.
    # RunList --------> Vector containing the run ID of each run.
    # time -----------> The time/run in seconds.
    # StepsDir -------> The current steps directory.
    #
    # Returns:
    # No returns, the function generates the data files in the data directory.

    file = open("$StepsDir/summary.dat", "w")

    write(file, "SIMULATED ANNEALING STRUCTURE SEARCH\n\n")
    write(file, "PARAMETERS:\n")
    write(file, "IT_scheme = $IT_scheme\n")
    write(file, "IT_params = $IT_params\n")
    write(file, "steps = $steps\n")
    write(file, "tempLength = $tempLength\n")
    write(file, "scheme = $scheme\n")
    write(file, "beta = $beta\n")
    write(file, "ScanningStyle = $ScanningStyle\n")
    write(file, "gamma = $gamma\n")
    write(file, "FirstToImprove = $FirstToImprove\n")
    write(file, "nRepeats = $nRepeats\n")
    write(file, "MaxStagnation = $MaxStagnation\n")
    write(file, "DisplayStyle = $DisplayStyle\n")
    write(file, "plots = $plots\n\n")

    write(file, "$accepted CONFIGURATIONS ACCEPTED:\n")
    if accepted > 0
        for i in 1:accepted
            e = EnergyList[i]
            n = RunList[i]
            write(file, "RUN #$n: E = $e eV/ATOM\n")
        end
        write(file, "\n")
        average   = mean(EnergyList)
        deviation = sqrt(var(EnergyList, corrected=true, mean=average))
        bestE     = minimum(EnergyList)
        bestRun   = argmin(EnergyList)
        worstE    = maximum(EnergyList)
        worstRun  = argmax(EnergyList)
        write(file, "AVERAGE ENERGY = $average ± $deviation eV/ATOM\n")
        write(file, "BEST ENERGY    = $bestE eV/ATOM [RUN #$bestRun]\n")
        write(file, "WORST ENERGY   = $worstE eV/ATOM [RUN #$worstRun]\n")
    end
    write(file, "\n")
    write(file, "OPTIMIZATION TIME:\n")
    write(file, "$time SECONDS/RUN\n")
    close(file)

    # screen
    #=
    println("\nSUMMARY:")
    println("STEPS, TEMPLENGTH, GAMMA = $steps, $tempLength, $gamma_abs")
    if accepted > 1
        println("ACCEPTED CONFIGURATIONS: $accepted $RunList")
        println("AVERAGE ENERGY         : $average ± $deviation eV/ATOM")
    else
        println("ACCEPTED CONFIGURATIONS: $accepted")
    end
    println("BEST CONFIGURATION     : RUN $bestRun ($bestE eV/ATOM)")
    println("")
    =#

end

####################################################################################################
# Functions devoted to moving particles and computing energies

function get_dE1(
    L::Int,
    config::Array{Int64,2},
    nbrs::Array{Int64,2},
    MoveStyle::Array{Int64,1},
    Nnbrs::Array{Int64,1},
    optSites::Array{Int,1},
    IonsSubset::Vector{Int64},
    U::Array{Float64, 2},
    UionAionBlist::Array{Float64,4},
    rng::MersenneTwister
    )

    # Function that proposes one move and computes its difference of energy.
    # 1st, an ion is selected at random.
    # 2nd, a move is proposed taking into account the move style of the selected ion (ion to
    #      vacancy / swap ions).
    # 3rd, the difference of energy is computed.
    # 4th, the move style, difference of energy, and configuration updates are returned.
    #
    # Arguments:
    # L -------------> Number of total interchangeable ions.
    # config --------> Matrix containing the solution (configuration of ions).
    # nbrs ----------> Matrix containing the neighborhood of the solution.
    # MoveStyle -----> Vector containing the move style of the ions belonging to each optimization
    #                  sites subset.
    # Nnbrs ---------> Vector containing the number of neighbors of the ions
    #                  belonging to each optimization sites subset.
    # optSites ------> Vector containing the optimization sites.
    # IonsSubset ----> Vector containing the subset (of sites) that each ion type
    #                  can occupy.
    # U -------------> Matrix [site s, type ion i] containing the interaction energy between an
    #                  ion of type i in site s with its anchored ion (fixed)
    #                  neighbors.
    # UionAionBlist -> Array [site A, site B, ion type i, ion type j] containing the
    #                  interaction energy between ion of type i in site A and
    #                  ion of type j in site B.
    # rng -----------> Random number generator.
    #
    # Returns:
    # style ---------> Moving style (1. ion to vacancy / 2. swap ions) of the proposed move.
    # dU ------------> Difference of energy of the proposed move.
    # configUpdates -> Tuple containing the configuration updates of the proposed move.

    # select an ion at random
    v = rand(rng, 1:L)
    # optimization site subset (type)
    type = IonsSubset[config[4,v]]
    # move style (1. ion-vacancy / 2. ion-ion)
    style = MoveStyle[type]
    if style == 1
        # select a vacancy at random
        w = rand(rng, 1:Nnbrs[type])
        # ion A goes from site ia to site ib
        ia   = config[2, v]
        ionA = config[4, v]
        ib   = nbrs[type, w]
        B    = optSites[ib]
        # configuration updates
        configUpdates = (v, ib, B, ionA, w, ia, type)
        # compute the difference of energy
        # 1st, energy of ion A with the anchored ions
        dU_InterAnch = -U[ia, ionA] + U[ib, ionA]
        # 2nd, energy of ion A with the rest of ions
        UiA = 0.0
        UiB = 0.0
        for l in 1:L
            k    = config[1,l]
            ilk  = config[2,l]
            ionX = config[4,l]
            if k != v # interaction A-A must be ignored
                UiA += UionAionBlist[ia,ilk,ionA,ionX]
                UiB += UionAionBlist[ib,ilk,ionA,ionX]
            end
        end
        dU_InterInter = UiB - UiA
        # finally, sum both contributions
        dU = dU_InterAnch + dU_InterInter
        # return results
        return style, dU, configUpdates
    elseif style == 2
        # ion 1
        ia   = config[2, v]
        A    = config[3, v]
        ionA = config[4, v]
        # ion 2
        repeat = true
        while repeat
            randnum = rand(rng, 1:Nnbrs[type])
            w = nbrs[type,randnum]
            if (w != v) && (config[4,v] != config[4,w])
                repeat = false
            end
        end
        ib   = config[2, w]
        B    = config[3, w]
        ionB = config[4, w]
        # configuration updates
        configUpdates = (v, ia, A, ionA, w, ib, B, ionB, type)
        # compute the difference of energy
        # 1st, energy of ions A and B with the anchored ions
        dU_InterAnch = (U[ib, ionA] + U[ia, ionB]) - (U[ia, ionA] + U[ib, ionB])
        # 2nd, energy of ions A and B with the rest of ions
        Ui0 = 0.0
        Ui1 = 0.0
        for l in 1:L
            k    = config[1,l]
            ilk  = config[2,l]
            ionX = config[4,l]
            if (k != v) && (k != w) # interactions A-A, A-B, B-B must be ignored
                Ui0 += UionAionBlist[ia,ilk,ionA,ionX] + UionAionBlist[ib,ilk,ionB,ionX]
                Ui1 += UionAionBlist[ib,ilk,ionA,ionX] + UionAionBlist[ia,ilk,ionB,ionX]
            end
        end
        dU_InterInter = Ui1 - Ui0
        # finally, sum both contributions
        dU = dU_InterAnch + dU_InterInter
        # return results
        return style, dU, configUpdates
    end
end

function get_InterAnchEnergy(
    L::Int,
    config::Array{Int64,2},
    U::Array{Float64,2})

    # Function that computes the total interchangeable-anchored ion pairs interaction energy of a
    # given configuration.
    #
    # Arguments:
    # L ------> Number of total interchangeable ions.
    # config -> Matrix containing the solution (configuration of ions).
    # U ------> Matrix [site s, type ion i] containing the interaction energy between an ion of
    #           type i in site s with its anchored ion (fixed) neighbors.
    #
    # Returns:
    # U_InterAnch -> Total interchangeable-anchored ion pairs interaction energy.

    U_InterAnch = 0.0
    for l in 1:L
        optSite = config[2,l]
        ionType = config[4,l]
        U_InterAnch += U[optSite,ionType]
    end
    return U_InterAnch
end

function get_InterInterEnergy(
    L::Int,
    config::Array{Int64,2},
    UionAionBlist::Array{Float64,4}
    )

    # Function that computes the total interchangeable-interchangeable ion pairs interaction energy
    # of a given configuration.
    #
    # Arguments:
    # L -------------> Number of total interchangeable ions.
    # config --------> Matrix containing the solution (configuration of ions).
    # UionAionBlist -> Array [site A, site B, ion type i, ion type j] containing the
    #                  interaction energy between ion of type i in site A and ion of type j
    #                  in site B.
    #
    # Returns:
    # U_InterInter -> Total ion-ion interaction energy.

    U_InterInter = 0.0
    for i in 1:(L-1)
        ia   = config[2,i]
        ionA = config[4,i]
        for j in (i+1):L
            ib   = config[2,j]
            ionB = config[4,j]
            U_InterInter += UionAionBlist[ia, ib, ionA, ionB]
        end
    end
    return U_InterInter
end

function get_TotalEnergy(
    config::Array{Int64,2},
    L::Int64,
    U,
    UionAionBlist::Array{Float64,4},
    EnergyBase::Float64
    )
    
    # Function that computes the total interaction energy of a given configuration.
    #
    # Arguments:
    # config --------> Matrix containing the solution (configuration of ions).
    # L -------------> Number of total interchangeable ions.
    # U -------------> Matrix [site s, type ion i] containing the interaction energy between an
    #                  ion of type i in site s with its anchored ion neighbors.
    # UionAionBlist -> Array [site A, site B, ion type i, ion type j] containing the
    #                  interaction energy between ion of type i in site A and ion of type j
    #                  in site B.
    # EnergyBase ----> Total interaction energy of the structure without interchangeable ions
    #                  (anchored-anchored ion pairs interaction).
    #
    # Returns:
    # e0 -> Total interaction energy.

    U_InterAnch  = get_InterAnchEnergy(L, config, U)                # interchangeable-anchored
    U_InterInter = get_InterInterEnergy(L, config, UionAionBlist)   # interchangeable-interchangeable
    e0           = EnergyBase + U_InterAnch + U_InterInter          # total

    return e0
end

function getNeighborsToScan(
    T::Float64,
    NeighborsToScan::Int64,
    ScanningStyle::Int64
    )

    # Function that returns the number of possible moves that must be scaned.
    # Several styles can be implemented.
    # Current styles:
    # Style 1 -> Constant value.
    #
    # Arguments:
    # T ---------------> Current temperature.
    # NeighborsToScan -> Initial number of moves to be scanned.
    # ScanningStyle ---> Method to decide how the number of combinations explored in each move
    #                    varies during the calculation.
    #
    # Returns:
    # nCfgNeighbors -> Number of moves to be scanned.

    if ScanningStyle == 1
        nCfgNeighbors = NeighborsToScan
    else
        println("ERROR: WRONG 'ScanningStyle' VALUE. CHECK THE 'params.jl' FILE.")
    end

    return nCfgNeighbors
end

function getBestFromScanNeighSimulation(
    L::Int,
    config::Array{Int64,2},
    nbrs::Array{Int64,2},
    MoveStyle::Array{Int64,1},
    Nnbrs::Array{Int64,1},
    optSites::Array{Int,1},
    IonsSubset::Vector{Int64},
    U::Array{Float64, 2},
    UionAionBlist::Array{Float64,4},
    nCfgNeighbors::Int,
    FirstToImprove::Bool,
    rng::MersenneTwister
    )
    
    # Function that tries "nCfgNeighbors" moves and returns the best (lowest difference of energy)
    # one.
    # Note: if FirstToImprove = true, the first move with negative difference of energy is returned.
    #
    # Arguments:
    # config ---------> Matrix containing the solution (configuration of ions).
    # nbrs -----------> Matrix containing the neighborhood of the solution.
    # MoveStyle ------> Vector containing the move style of the ions belonging to each optimization
    #                   sites subset.
    # Nnbrs ----------> Vector containing the number of neighbors of the ions belonging to each
    #                   optimization sites subset.
    # optSites -------> Vector containing the optimization sites.
    # IonsSubset -----> Vector containing the subset (of sites) that each ion type can occupy.
    # U --------------> Matrix [site s, type ion i] containing the interaction energy between an
    #                   ion of type i in site s with its anchored atom (fixed) neighbors.
    # UionAionBlist --> Array [site A, site B, ion type i, ion type j] containing the
    #                   interaction energy between ion of type i in site A and ion of type j
    #                   in site B (integers).
    # nCfgNeighbors --> Number of moves to be scanned.
    # FirstToImprove -> Skip or continue the exploration of moves when a successfull moves is
    #                   found.
    # rng ------------> Random number generator.
    #
    # Returns:
    # style_best ---------> Style of the selected optimization sites subset (1 -> vacancies /
    #                       2 -> no vacancies).
    # dE_best ------------> Difference of energy of the selected move.
    # configUpdates_best -> Tuple containing the configuration updates of the selected move.

    # this is the first neighbor:
    style, dE, configUpdates = get_dE1(L, config, nbrs, MoveStyle, Nnbrs, optSites, IonsSubset, U,
                                       UionAionBlist, rng)
    style_best = style
    dE_best = dE
    configUpdates_best = configUpdates
    # the following neighbors:
    for j in 1:nCfgNeighbors - 1
        if FirstToImprove # if FirstToImprove = true, select the first move with dE < 0
            if dE_best < 0
                return style_best, dE_best, configUpdates_best
            end
        end
        style, dE, configUpdates = get_dE1(L, config, nbrs, MoveStyle, Nnbrs, optSites, IonsSubset, U, UionAionBlist, rng)
        if dE < dE_best
            style_best  = style
            dE_best     = dE
            configUpdates_best = configUpdates
        end
    end
    return style_best, dE_best, configUpdates_best
end

function metropolis(dE::Float64, T::Float64, rng::MersenneTwister)
    
    # Function that accepts/rejects a move with a energy difference dE at temperature T using the
    # metropolis criterion.
    #
    # Arguments:
    # dE --> Difference of energy of the proposed move.
    # T ---> Current temperature.
    # rng -> Random number generator.
    #
    # Returns:
    # true/false -> Accept/Reject

    if dE <= 0.0
        return true
    elseif exp(-dE / T) > rand(rng)
        return true
    else
        return false
    end
end

function move_simulation!(
    L::Int,
    config::Array{Int64,2},
    nbrs::Array{Int64,2},
    MoveStyle::Array{Int64,1},
    Nnbrs::Array{Int64,1},
    p_accepted::Array{Bool,1},
    p_dEsum::Array{Float64,1},
    optSites::Array{Int,1},
    IonsSubset::Vector{Int64},
    U::Array{Float64, 2},
    UionAionBlist::Array{Float64,4},
    T::Float64,
    nCfgNeighbors::Int,
    FirstToImprove::Bool,
    rng::MersenneTwister)
    
    # Function that gets the best move over "nCfgNeighbors" tries, selects/rejects the move using
    # the metropolis criterion and, if the move is accepted, updates the ion configuration.
    #
    # Arguments:
    # L --------------> Number of total interchangeable ions.
    # config ---------> Matrix containing the solution (configuration of ions).
    # nbrs -----------> Matrix containing the neighborhood of the solution.
    # MoveStyle ------> Vector containing the move style of the ions belonging to each optimization
    #                   sites subset.
    # Nnbrs ----------> Vector containing the number of neighbors of the ions belonging to each
    #                   optimization sites subset.
    # p_accepted -----> Single element vector containing the accepted/rejected status (true/false)
    #                   of the proposed move. 
    # p_dEsum --------> Single element vector containing the energy of the current configuration.
    # optSites -------> Vector containing the optimization sites.
    # IonsSubset -----> Vector containing the subset (of sites) that each ion type can occupy.
    # U --------------> Matrix [site s, type ion i] containing the interaction energy between an
    #                   ion of type i in site s with its anchored ion (fixed) neighbors.
    # UionAionBlist --> Array [site A, site B, ion type i, ion type j] containing the
    #                   interaction energy between ion of type i in site A and ion of type j
    #                   in site B.
    # T --------------> Current temperature.
    # nCfgNeighbors --> Number of moves to be scanned.
    # FirstToImprove -> Skip or continue the exploration of moves when a successfull moves is found.
    # rng ------------> Random number generator.
    #
    # Returns:
    # No returns, the function updates the arguments.

    # try nCfgNeighbors new configurations at random and select the best one
    style, dE, configUpdates = getBestFromScanNeighSimulation(L, config, nbrs, MoveStyle,
                                                              Nnbrs, optSites, IonsSubset, U,
                                                              UionAionBlist, nCfgNeighbors,
                                                              FirstToImprove, rng)
    # determine if the new configuration is accepted
    accepted = metropolis(dE, T, rng)
    if accepted
        # update configuration
        if style == 1
            # move ion to vacancy (changes = v, ib, B, ionA, w, ia)
            v, ib, B, ionA, w, ia, type = configUpdates
            config[1,v]  = v
            config[2,v]  = ib
            config[3,v]  = B
            config[4,v]  = ionA
            nbrs[type,w] = ia
        elseif style == 2
            # swap two ions (changes = v, ia, A, ionA, w, ib, B, ionB)
            v, ia, A, ionA, w, ib, B, ionB, type = configUpdates
            config[1,v] = v
            config[2,v] = ib
            config[3,v] = B
            config[4,v] = ionA
            config[1,w] = w
            config[2,w] = ia
            config[3,w] = A
            config[4,w] = ionB
        end
        # update energy
        p_dEsum[1] += dE
    end
    p_accepted[1] = accepted
end

function returnToPreviousOpt!(
    config::Array{Int64,2},
    nbrs::Array{Int64,2},
    config_opt::Array{Int64,2},
    nbrs_opt::Array{Int64,2},
    isOptimal::Array{Bool,1},
    dEsum_opt::Array{Float64,1},
    dEsum::Array{Float64,1},
    e0::Array{Float64,1},
    E::Array{Float64},
    moves::Int
    )

    # Function that checks if the current configuration is optimal, and if it is false, it returns
    # to the optimal solution (best solution found).
    #
    # Arguments:
    # config -----> Matrix containing the current solution (configuration of ions).
    # nbrs -------> Matrix containing the neighborhood of the current solution.
    # config_opt -> Matrix containing the optimal solution (configuration of ions).
    # nbrs_opt ---> Matrix containing the neighborhood of the optimal solution.
    # isOptimal --> True/false if the current configuration is/is not optimal.
    # dEsum_opt --> 1-element vector containing the optimal difference of energy.
    # dEsum ------> 1-element vector containing the current difference of energy.
    # e0 ---------> 1-element vector containing the current energy.
    # E ----------> Vector containing the energy evolution.
    # moves ------> Number of moves.
    #
    # Returns:
    # No returns, the function updates the arguments.

    if !isOptimal[1]
        config .= config_opt
        nbrs .= nbrs_opt
        dEsum[1] = dEsum_opt[1]
        E[moves] = e0[1] + dEsum[1] # impose last E = optimal E
    end
end

####################################################################################################
# Functions devoted to initializing temperature

function RandomWalk(
    L::Int64,
    MoveStyle::Array{Int,1},
    Nnbrs::Array{Int,1},
    optSites::Array{Int,1},
    IonsSubset::Vector{Int64},
    U,
    UionAionBlist::Array{Float64,4},
    EnergyBase::Float64,
    steps::Int64,
    config0::Array{Int64,2},
    nbrs0::Array{Int64,2},
    rng::MersenneTwister)
    
    # Function that performs a Random Walk with a given configuration of ions and vacancies and
    # returns the average, max, and min difference of energy between two consecutive steps.
    #
    # Arguments:
    # L -------------> Number of total interchangeable ions.
    # MoveStyle -----> Vector containing the move style of the ions belonging to each optimization
    #                   sites subset.
    # Nnbrs ---------> Vector containing the number of neighbors of the ions belonging to each
    #                   optimization sites subset.
    # optSites ------> Vector containing the optimization sites.
    # IonsSubset ----> Vector containing the subset (of sites) that each ion type can occupy.
    # U -------------> Matrix [site s, type ion i] containing the interaction energy between an
    #                  ion of type i in site s with its anchore ion (fixed) neighbors.
    # UionAionBlist -> Array [site A, site B, ion type i, ion type j] containing the
    #                  interaction energy between ion of type i in site A and ion of type j
    #                  in site B.
    # EnergyBase ----> Total interaction energy of the structure without interchangeable ions
    #                  (anchore-anchore ion pairs interaction).
    # steps ---------> Random Walk steps.
    # config0 -------> Matrix containing the solution (configuration of ions).
    # nbrs0 ---------> Matrix containing the neighborhood of the solution.
    # rng -----------> Random number generator.
    #
    # Returns:
    # mean_ΔE -> Mean difference of energy between two consecutive steps.
    # max_ΔE --> Max difference of energy between two consecutive steps.
    # min_ΔE --> Min difference of energy between two consecutive steps.

    ########################################################################

    # get initial solution and neighborhood
    config = copy(config0)
    nbrs = copy(nbrs0)
    # get initial energy (e0)
    e0 = get_TotalEnergy(config, L, U, UionAionBlist, EnergyBase)
    # set energy of each step (E) and difference of energy (dEsum) to zero
    E = zeros(steps)
    dEsum = 0.0

    ########################################################################

    # start loop of random moves
    for i in 1:steps

        # propose a move and get de difference of energy and the style (ion->vacancy / ion <-> ion)
        style, dE, configUpdates = get_dE1(L, config, nbrs, MoveStyle, Nnbrs, optSites, IonsSubset,
                                           U, UionAionBlist, rng)

        # update configuration
        if style == 1
            # move ion to vacancy (changes = v, ib, B, ionA, w, ia)
            v, ib, B, ionA, w, ia, type = configUpdates
            config[1,v]  = v
            config[2,v]  = ib
            config[3,v]  = B
            config[4,v]  = ionA
            nbrs[type,w] = ia
        elseif style == 2
            # swap two ions (changes = v, ia, A, ionA, w, ib, B, ionB)
            v, ia, A, ionA, w, ib, B, ionB, type = configUpdates
            config[1,v] = v
            config[2,v] = ib
            config[3,v] = B
            config[4,v] = ionA
            config[1,w] = w
            config[2,w] = ia
            config[3,w] = A
            config[4,w] = ionB
        end

        # update energy
        dEsum += dE
        E[i]   = e0 + dEsum

    end
    # end loop of random moves

    # statistics
    ΔE = [abs(E[i+1] - E[i]) for i = 1:(steps-1)]
    mean_ΔE = mean(ΔE)
    max_ΔE  = maximum(ΔE)
    min_ΔE  = minimum(ΔE)

    # return statistics
    return mean_ΔE, max_ΔE, min_ΔE
end

function InitialTemperature(
    IT_scheme::Int64,
    IT_params::Vector{Float64},
    data::NamedTuple,
    IonSiteConstraint::Bool
    )

    # Function that computes the max (initial) and min (final) temperature of the main cooling
    # scheme.
    # Avalable schemes:
    #   IT_scheme = 1 -> OIT1
    #                    The user decides the values of Tmax and Tmin.
    #                    Tmax = IT_params[1] / Tmin = IT_params[2]
    #   IT_scheme = 2 -> OIT2
    #                    The values of Tmax/Tmin are proportional to the max/min difference of
    #                    energy between two consecutive steps of a Random Walk performed with the
    #                    structure that we want to optimize. The initial configuration of ions and
    #                    vacancies is generated at random.
    #                    Random Walk steps = IT_params[1]
    #                    Tmax = IT_params[2] * ΔE_max ; Tmin = IT_params[2] * ΔE_min
    #   IT_scheme = 3 -> OIT3
    #                    The values of Tmax/Tmin are adjusted to match a certain max/min prob. of
    #                    acceptance.
    #                    The prob. of acceptance is computed as p = exp(-ΔE/T).
    #                    The difference of energy ΔE is the average difference of energy between
    #                    two consecutive steps of a Random Walk performed with the structure that
    #                    we want to optimize. The initial configuration of ions and vacancies is
    #                    generated at random.
    #                    Random Walk steps = IT_params[1]
    #                    Prob_max = IT_params[2]; Prob_min = IT_params[3]
    #                    Tmax = - ΔE / log(Prob_max) ; Tmin = - ΔE / log(Prob_min)
    #
    # Arguments:
    # IT_scheme -> Scheme to be used.
    # IT_params -> Vector containing the parameters required by the used scheme.
    # data ------> Data tuple (defined in ImportData.get_data(...)).
    #
    # Returns:
    # Tmax -> Max (initial) temperature for the main cooling scheme.
    # Tmin -> Min (final) temperature for the main cooling scheme.

    if IT_scheme == 1 # OIT1
        Tmax = IT_params[1]
        Tmin = IT_params[2]
    else # related to the differences of energy between two consecutive steps of a random walk
        # structure data
        @unpack Ions, LIons, IonsSubset, NoptSites, Na, optSites, SitesSubset, U,
                UionAionBlist, EnergyBase = data
        # total number of interchangeable ions
        L = sum(LIons)
        # initialize random number generator
        rng_rw = MersenneTwister(1)
        # fill optimization sites at random
        if IonSiteConstraint == false # without constraints
            config = get_init_sol(LIons, NoptSites, optSites, rng_rw)
            nbrs, MoveStyle, Nnbrs = get_nbrs(config, L, NoptSites)
        else # with constraints
            config = get_init_sol_constrained(LIons, NoptSites, optSites, 
                                              SitesSubset, IonsSubset, rng_rw)
            nbrs, MoveStyle, Nnbrs = get_nbrs_constrained(config, SitesSubset,
                                                          IonsSubset, NoptSites)
        end

        # random walk
        steps = round(Int64, IT_params[1])
        mean_ΔE, max_ΔE, min_ΔE = RandomWalk(L, MoveStyle, Nnbrs, optSites, IonsSubset, U,
                                             UionAionBlist, EnergyBase, steps, config, nbrs, rng_rw)
        if IT_scheme == 2 # OIT2
            Tmax = IT_params[2] * max_ΔE
            Tmin = IT_params[2] * min_ΔE
        elseif IT_scheme == 3 # OIT3
            Tmax = - mean_ΔE / log(IT_params[2])
            Tmin = - mean_ΔE / log(IT_params[3])
        end

        # display random walk summary
        println("RANDOM WALK:")
        println("ΔE_MAX = $max_ΔE")
        println("ΔE_MIN = $min_ΔE")
        println("ΔE_AVE = $mean_ΔE")
        println("")
    end

    # display initial temperature
    println("INITIAL TEMPERATURE:")
    println("TMAX = $Tmax")
    println("TMIN = $Tmin")
    println("")

    return Tmax, Tmin
end

function getListT!(T0::Float64,
    beta::Float64,
    tempLength::Int64,
    scheme::String,
    listT::Array{Float64,1})
    
    # Function that generates a vector containing the temperature values of the secondary loop.
    # The function reads the reference temperature of the main cooling scheme (T0), and computes the
    # secondary cooling scheme.
    # The currently available schemes are:
    #   scheme = "constant" -> No cooling, the reference temperature is used at each step of the
    #                          secondary cooling
    #   scheme = "linear" ---> Linear cooling, the temperature decreases linearly from T0 to
    #                          T0 x beta.
    #
    # Arguments:
    # T0 ---------> Reference temperature of the main cooling scheme.
    # beta -------> Intensity of the secondary cooling scheme: Tf = T0 x beta.
    # tempLength -> Number of steps (temperatures) of the secondary cooling scheme.
    # scheme -----> Type of secondary cooling.
    # listT ------> Vector containing the temperatures of the secondary cooling scheme.
    #
    # Returns:
    # No returns, the function updates the value of the arguments (listT). 

    if scheme == "constant"
        for i in 1:tempLength
            listT[i] = T0
        end
    elseif scheme == "linear"
        Tf = T0 * beta
        dT = (Tf - T0) / tempLength
        for i in 1:tempLength
            listT[i] =  T0 + ( (i - 1) * dT)
        end
    end
end

####################################################################################################
# Functions devoted to initializing the random configuration of ions/vacancies

function get_init_sol_constrained(
    LIons::Vector{Int64},
    NoptSites::Int,
    optSites::Array{Int64,1},
    SitesSubset::Array{Int64,1},
    IonsSubset::Array{Int64,1},
    rng::MersenneTwister)

    # Function that generates a (constrained) random configuration of ions.
    # That is, generates a matrix, where each column corresponds to an ion (ith ion -> ith column).
    # The matrix (config) has 4 rows, each row determine the following information:
    # config[1,:] ion ID
    # config[2,:] optimization site ID (where the ion is located)
    # config[3,:] site ID (where the ion is located) <- it is equivalent to the optimization site ID
    #                                                   (but it is convenient to have both stored)
    # config[4,:] type of ion
    # 
    # Arguments:
    # LIons -------> Vector containing the number of ions of each type of ion.
    # NoptSites ---> Number of optimization sites.
    # optSites ----> Vector containing the optimization sites.
    # SitesSubset -> Vector containing the subset (of sites) to which each site belongs.
    # IonsSubset --> Vector containing the subset (of sites) that each ion type can occupy.
    # rng ---------> Random number generator.
    #
    # Returns:
    # config -> Matrix containing a random solution (configuration of ions).

    # total number of interchangeable ions
    L = sum(Int64, LIons)

    # init arrays
    empty_sites = [i for i in 1:NoptSites]
    config = zeros(Int, 4, L)
    ion = 1; count = 0

    for l in 1:L

        # determine the ion type
        count += 1
        if count > LIons[ion]
            ion  += 1
            count = 1
        end

        # get the optimization site subset
        required_subset = IonsSubset[ion]

        # find the optimization sites that belong to this subset
        possible_sites = findall(x->x==required_subset, SitesSubset)

        # find vacancies (available sites) on these optimization sites
        available_sites = intersect(possible_sites, empty_sites)
        
        # select an available site at random
        r = rand(rng, available_sites)

        # remove the selected site from the empty sites list
        filter!(x->x!=r, empty_sites)

        # config
        config[1,l] = l             # ion id
        config[2,l] = r             # optimization site id
        config[3,l] = optSites[r]   # site id
        config[4,l] = ion           # type of ion

    end
    return config
end

function get_init_sol(
    LIons::Vector{Int64},
    NoptSites::Int,
    optSites::Array{Int64,1},
    rng::MersenneTwister)

    # Function that generates a (constrained) random configuration of ions.
    # That is, generates a matrix, where each column corresponds to an ion (ith ion -> ith column).
    # The matrix (config) has 4 rows, each row determine the following information:
    # config[1,:] ion ID
    # config[2,:] optimization site ID (where the ion is located)
    # config[3,:] site ID (where the ion is located) <- it is equivalent to the optimization site ID
    #                                                   (but it is convenient to have both stored)
    # config[4,:] type of ion
    # 
    # Arguments:
    # LIons -------> Vector containing the number of ions of each type of ion.
    # NoptSites ---> Number of optimization sites.
    # optSites ----> Vector containing the optimization sites.
    # SitesSubset -> Vector containing the subset (of sites) to which each site belongs.
    # IonsSubset --> Vector containing the subset (of sites) that each ion type can occupy.
    # rng ---------> Random number generator.
    #
    # Returns:
    # config -> Matrix containing a random solution (configuration of ions).

    L = sum(Int64,LIons)

    # random configuration
    il = sample(rng,1:NoptSites,L,replace=false)

    config = zeros(Int,4,L)
    ion = 1; count = 0
    for (l,r) in enumerate(il)

        # determine the ion type
        count += 1
        if count > LIons[ion]
            ion  += 1
            count = 1
        end

        # config
        config[1,l] = l             # ion id
        config[2,l] = r             # optimization site id
        config[3,l] = optSites[r]   # site id
        config[4,l] = ion           # type of ion

    end
    return config
end

function get_nbrs_constrained(config::Array{Int64,2},
                                  SitesSubset::Array{Int64,1},
                                  IonsSubset::Array{Int64,1},
                                  NoptSites::Int)
    
    # Function that generates the neighborhood of a given solution (configuration of ions), when
    # the ion-site constraint is activated (there are different optimization sites subsets).
    # That is, generates a matrix containing the neighbors of the ions belonging to each
    # optimization sites subset.
    # The 1st dimension (row) determines the the optimization sites subset.
    # Style 1: when the subset contains vacancies (more sites than ions -> NoptSites[i] > L[i]),
    #          the row indicates the vacancies belonging to the subset, so each element of that row
    #          indicates the site where a vacancy is located.
    #          The moves of ions belonging to style 1 subsets consist of moving an ion to a vacancy.
    # Style 2: when the subset contains no vacancies (NoptSites[i] = L[i]), the row indicates the
    #          ions belonging to the subset, so each element of that row indicates the site where
    #          an ion is located.
    #          The moves of ions belonging to style 2 subsets consists of swaping two ions.
    # The number of rows of the matrix is equal to the number of subsets.
    # The number of columns of the matrix is equal to the number of neighbors of the subset with
    # the larger neighborhood. The extra columns of the subsets with smaller neighborhoods are
    # filled with zeros, which are ignored during the optimization.
    # The function also generates two auxiliar vectors: one to identify the "style" of each subset
    # (row), and one to identify the number of neighbors (non-zero elements) of each subset.
    # 
    # Arguments:
    # config ------> Matrix containing the random configuration of ions.
    # SitesSubset -> Vector containing the subset (of sites) to which each site belongs.
    # IonsSubset --> Vector containing the subset (of sites) that each ion type can occupy.
    # NoptSites ---> Number of optimization sites.
    #
    # Returns:
    # nbrs ------> Matrix containing the neighborhood of the given solution (configutaion of ions).
    # MoveStyle -> Vector containing the move style of the ions belonging to each optimization sites
    #              subset.
    # Nnbrs -----> Vector containing the number of neighbors of the ions belonging to each
    #              optimization sites subset.

    # total number of interchangeable ions
    L = length(config[1,:])

    # number of optimization sites subsets
    Nsubsets = maximum(IonsSubset)

    # subset to which each ion belongs
    IonsSubset_list = IonsSubset[config[4,:]]

    nbrs_aux = []
    for i in 1:Nsubsets # loop over subsets

        # number of sites belonging to the current subset
        Nsites = count(x->x==i, SitesSubset)

        # number of ions belonging to the current subset
        Nions = count(x->x==i, IonsSubset_list)

        subset_nbrs = Int64[]
        if Nions < Nsites # there are vacancies
            style = 1 # store vacancies
            for r in 1:NoptSites # loop over optimization sites
                if r ∉ config[2,:] # only select vacancies belonging to the current subset
                    if i == SitesSubset[r]
                        append!(subset_nbrs, r)
                    end
                end
            end
        else # there are no vacancies
            style = 2 # store ions
            for ion in 1:L # loop over ions
                if IonsSubset[config[4,ion]] == i # only select ions belonging to the current subset
                    append!(subset_nbrs, ion)
                end
            end
        end
        push!(nbrs_aux, [style, subset_nbrs])
    end

    # write in matrix style
    nrows = length(nbrs_aux) # number of rows = number of subsets
    Nnbrs = zeros(Int64, nrows)
    for i in eachindex(nbrs_aux) # loop over subsets
        Nnbrs[i] = length(nbrs_aux[i][2]) # number of neighbors (vacancies or ion)
    end
    ncols = maximum(Nnbrs) # number of columns = max(number of neighbors)
    nbrs = zeros(Int64, nrows, ncols) # matrix of neighbors
    MoveStyle = zeros(Int64, nrows) # vector of move styles (1 or 2)
    for i in eachindex(nbrs_aux) # loop over subsets
        MoveStyle[i] = nbrs_aux[i][1] # set the style of each subset
        for j in 1:Nnbrs[i] # loop over neighbors
            nbrs[i,j] = nbrs_aux[i][2][j] # fill the matrix elements with the neighbors
        end
    end

    return nbrs, MoveStyle, Nnbrs
end

function get_nbrs(config::Array{Int64,2},
                       L::Int,
                       NoptSites::Int)
    

    # Function that generates the neighborhood of a given solution (configuration of ions), when
    # the ion-site constraint is deactivated (there are no subsets).
    # Same description as get_nbrs_constrained(...) but without splitting the sites/ions in
    # different subsets. In this case, the matrix is a row matrix (1 row), and any ion can visit any
    # optimization site. 
    # 
    # Arguments:
    # config ------> Matrix containing the random configuration of ions.
    # SitesSubset -> Vector containing the subset (of sites) to which each site belongs.
    # L -----------> Number of total interchangeable ions.
    # NoptSites ---> Number of optimization sites.
    #
    # Returns:
    # nbrs ------> Matrix containing the neighborhood of the given solution (configutaion of ions).
    # MoveStyle -> Vector containing the move style of the ions belonging to each optimization sites
    #              subset.
    # Nnbrs -----> Vector containing the number of neighbors of the ions belonging to each
    #              optimization sites subset.

    subset_nbrs = Int64[] # only 1 subset
    if L < NoptSites # there are vacancies
        style = 1
        for r in 1:NoptSites # loop over optimization sites
            if r ∉ config[2,:] # only select vacancies
                append!(subset_nbrs, r)
            end
        end
    else # there are no vacancies
        style = 2
        for ion in 1:L # loop over ions
            append!(subset_nbrs, ion)
        end
    end
    nbrs_aux = [[style, subset_nbrs]]

    # write in matrix style
    Nnbrs = zeros(Int64, 1) # unique subset
    Nnbrs[1] = length(nbrs_aux[1][2])
    ncols = Nnbrs[1] # number of columns = number of neighbors
    nbrs = zeros(Int64, 1, ncols) # unique subset -> 1 row
    MoveStyle = zeros(Int64, 1)
    MoveStyle[1] = nbrs_aux[1][1] # set the style
    for j in 1:Nnbrs[1]
        nbrs[1,j] = nbrs_aux[1][2][j] # fill the matrix elements with the neighbors
    end

    return nbrs, MoveStyle, Nnbrs
end

####################################################################################################
# SANSS inner loop

function SANSS_opt(
    L::Int64,
    Na::Int64,
    gamma::Int64,
    ScanningStyle::Int64,
    FirstToImprove::Bool,
    optSites::Array{Int,1},
    IonsSubset::Vector{Int64},
    U,
    UionAionBlist::Array{Float64,4},
    EnergyBase::Float64,
    Tmax::Float64,
    Tmin::Float64,
    beta::Float64,
    tempLength::Int64,
    scheme::String,
    steps::Int64,
    MaxStagnation::Int64,
    config::Array{Int64,2},
    nbrs::Array{Int64,2},
    MoveStyle::Array{Int64,1},
    Nnbrs::Array{Int64,1},
    run::Int64,
    RunDir::String,
    DisplayStyle::Int64,
    rng::MersenneTwister
    )
    
    # Function that gets the initial configuration and its neighborhood, the system information
    # and the SA parameters, and computed the optimal configuration of ions using the SA algorithm.
    # Scheme:
    # 1st, the necessary arrays and scalars are initializated.
    # 2nd, the main cooling scheme loop starts.
    # 3rd, the secondary cooling scheme temperatures are computed (getListT!(...)).
    # 4th, the secondary cooling scheme loop starts.
    # 5th, a move is performed (move_simulation!(...)).
    # 6th, if the move is accepted, update some arrays and scalars.
    # 7th, the secondary cooling scheme loop ends.
    # 8th, the main cooling scheme temperature is updated (geometric cooling).
    # 9th, if the last configuration is not optimal, it is replaced by the optimal one.
    # 10th, the main cooling loop ends.
    # 11th, the optimal configuration (best configuration found) is returned.
    #
    # Arguments:
    # L --------------> Number of total interchangeable ions.
    # Na -------------> Number of total atoms/ions (interchangeable + anchored ions).
    #                   Na = number of sites - number of vacancies
    # gamma ----------> Max number of combinations explored in each move.
    # ScanningStyle --> Method to decide how the number of combinations explored in each move
    #                   varies during the calculation.
    # FirstToImprove -> Skip or continue the exploration of moves when a successfull moves is
    #                   found.
    # optSites -------> Vector containing the optimization sites.
    # SitesSubset ----> Vector containing the subset (of sites) to which each site belongs.
    # IonsSubset -----> Vector containing the subset (of sites) that each ion type can occupy.
    # U --------------> Matrix [site s, type ion i] containing the interaction energy between an
    #                   ion of type i in site s with its anchored ion (fixed) neighbors.
    # UionAionBlist --> Array [site A, site B, ion type i, ion type j] containing the
    #                   interaction energy between ion of type i in site A and ion of type j
    #                   in site B.
    # EnergyBase -----> Total interaction energy of the structure without interchangeable ions.
    # Tmax -----------> Initial temperature of the first step of the main cooling scheme.
    # Tmin -----------> Initial temperature of the last step of the main cooling scheme.
    # beta -----------> Intensity of the secondary cooling scheme. T_final = beta * T_initial.
    # tempLength -----> Number of steps in the secondary cooling scheme.
    # scheme ---------> Type of secondary cooling scheme.
    # steps ----------> Number of steps in the main cooling scheme.
    # config ---------> Matrix containing the solution (configuration of ions).
    #                   Note: To understand how it is constructed, go to the function
    #                         get_init_sol_constrained(...).
    # nbrs -----------> Matrix containing the neighborhood of the solution.
    #                   Note: To understand how it is constructed, go to the function
    #                         get_nbrs_constrained(...).
    # MoveStyle ------> Vector containing the move style of the ions belonging to each optimization
    #                   sites subset.
    #                   Note: To understand how it is constructed, go to the function
    #                         get_nbrs_constrained(...).
    # Nnbrs ----------> Vector containing the number of neighbors of the ions belonging to each
    #                   optimization sites subset.
    #                   Note: To understand how it is constructed, go to the function
    #                         get_nbrs_constrained(...). 
    # run ------------> Run (optimization cycle) ID.
    # RunDir ---------> Run directory.
    # DisplayStyle ---> Display style.
    # rng ------------> Random number generator.
    #
    # Returns:
    # config ---> Matrix containing the best solution (configuration of ions) found.
    # E[moves] -> Energy of the best solution found.

    # initialize scalars
    maximumMoves::Int64    = steps * tempLength 
    OptConfig::Bool        = true
    keepLooping::Bool      = true
    stagnation::Int64      = 0
    moves::Int64           = 0
    iter::Int64            = 0
    naccepted::Int64       = 0
    lastAccepted::Int64    = 0
    alpha::Float64         = (Tmin / Tmax) ^ (1 / (steps - 1))
    NeighborsToScan::Int64 = gamma

    # initialize arrays
    accepted      = zeros(Bool, 1)
    isOptimal     = zeros(Bool, 1)
    acceptance    = zeros(Float64, 0)
    lastE         = zeros(Float64, 0)
    bestE         = zeros(Float64, 0)
    sampled_bestE = zeros(Float64, 0)
    E_evol        = zeros(Float64, tempLength)
    listT         = zeros(Float64, tempLength)
    record_T      = zeros(Float64, maximumMoves)
    E             = zeros(Float64, maximumMoves)
    T_opt         = zeros(Float64, 1)

    # set best configuration to initial configuration
    nbrs_opt = copy(nbrs)
    config_opt = copy(config)

    # compute energy of initial configuration
    e0 = [get_TotalEnergy(config, L, U, UionAionBlist, EnergyBase)]

    # set initial energy difference and best energy difference to 0
    dEsum     = zeros(Float64, 1)
    dEsum_opt = zeros(Float64, 1)

    # initialize reference temperature
    T0::Float64 = Tmax

    if DisplayStyle > 0
        println(@sprintf("\nRUN %03d: RANDOM CONFIGURATION => E = %12.6e EV/ATOM", run, e0[1]/Na))
    end

    for step in 1:steps # start main cooling loop
        if keepLooping

            # generate temperatures for the secondary cooling loop
            getListT!(T0, beta, tempLength, scheme, listT)

            # reset variables
            OptConfig = false
            iter      = 0
            naccepted = 0
            E_evol   .= 0

            for T in listT # start secondary cooling loop
                if keepLooping

                    # update variables
                    moves += 1
                    iter  += 1
                    accepted[1] = false
                    if moves > 1
                        Eold = E[moves - 1]
                    else
                        Eold = e0[1]
                    end

                    # compute the number of possible moves that must be explored
                    nCfgNeighbors = getNeighborsToScan(T, NeighborsToScan, ScanningStyle)

                    # perform one move
                    move_simulation!(L, config, nbrs, MoveStyle, Nnbrs, accepted, dEsum,
                                     optSites, IonsSubset, U, UionAionBlist, T, nCfgNeighbors,
                                     FirstToImprove, rng)

                    if accepted[1]
                        # update system
                        lastAccepted = moves
                        naccepted += 1
                        if dEsum[1] < dEsum_opt[1]
                            # update system
                            isOptimal[1] = true
                            OptConfig = true
                            stagnation = 0
                            dEsum_opt[1] = dEsum[1]
                            T_opt[1] = T
                            config_opt[:,:] = config[:,:]
                            nbrs_opt[:,:] = nbrs[:,:]
                            # save optimal configuration energy
                            append!(sampled_bestE, e0[1] + dEsum[1])
                            # show results by screen
                            if DisplayStyle == 1
                                println(@sprintf("RUN %03d: STEP = %010d, T = %12.6e, E = %12.6e EV/ATOM", run, step, T, sampled_bestE[end]/Na))
                            end
                        else
                            isOptimal[1] = false
                        end
                    else
                        isOptimal[1] = false
                    end
                    
                    # update system 
                    E[moves] = e0[1] + dEsum[1]
                
                    # save temperature evolution
                    record_T[moves] = T

                    # check minimum energy
                    if T <= 1.0
                        keepLooping = false
                    end

                    # save energy evolution
                    E_evol[iter] = e0[1] + dEsum[1]
                end
            end # end secondary cooling loop

            # update stagnation
            if !OptConfig
                stagnation += 1
            end

            # in the 1st step, the best energy is the optimal energy by definition
            if step == 1
                append!(sampled_bestE , minimum(E_evol))
            end

            # save results after secondary cooling loop (acceptance, last energy, best energy)
            append!(acceptance, naccepted/tempLength)
            append!(lastE, E_evol[iter])
            append!(bestE, sampled_bestE[end])

            # return to the previous optimal configuration
            returnToPreviousOpt!(config, nbrs, config_opt, nbrs_opt, isOptimal, dEsum_opt, dEsum, e0, E, moves)

            # show results by screen
            if DisplayStyle == 2
                println(@sprintf("RUN %03d: STEP = %010d, T0 = %12.6e, E = %12.6e EV/ATOM, STAGNATION = %010d",
                                 run, step, T0, bestE[step]/Na, stagnation))
            end

            # if keepLooping = false (T < 1) display message
            if keepLooping == false && DisplayStyle > 0
                println(@sprintf("RUN %03d: TEMPERATURE < 1 => STOP OPTIMIZATION", run))
            end

            # check stagnation
            if MaxStagnation > 0
                if (stagnation == MaxStagnation) && (keepLooping == true)
                    keepLooping = false
                    if DisplayStyle > 0
                        println("LIMIT OF STAGNATION REACHED => STOP OPTIMIZATION")
                    end
                end
            end

            # update main cooling temperature (geometric cooling scheme)
            T0 *= alpha

        end

    end # end main cooling loop

    # show results by screen
    if (DisplayStyle == 0)
        step = length(bestE)
        println(@sprintf("\nRUN %03d: FINAL ENERGY = %12.6e EV/ATOM", run, bestE[step] / Na))
    end

    # generate optimization files
    write_E(sampled_bestE/Na, lastE/Na, bestE/Na, acceptance, RunDir)

    # return optimal configuration and its energy
    return config, E[moves]
end

####################################################################################################
# SANSS outer loop (main function)

function SANSS_set(data::NamedTuple, params::NamedTuple)

    # Function that gets the data tuple (check ImportData module for more details) and params tuple 
    # (check SelectSolver module for more details), generates the initial random configuration of
    # ions and vacancies and optimizes the configuration using the SA algorithm.
    # Scheme:
    # 1st, the initial and final temperatures to be used in the SA algorithm are computed.
    # 2nd, the runs loop (which can be executed in parallel) starts.
    # 3rd, an initial random configuration is generated for each run.
    # 4th, the optimal configuration is computed (check SANSS_opt) for each run (if steps = 0, the
    #      random configuration is not optimized).
    # 5th, the runs loop ends. 
    # 6th, some statistics are performed.
    # 7th, the optimal configurations are returned.
    #
    # Arguments:
    # data ---> Named tuple (the data tuple) defined in ImportData.get_data(...).
    # params -> Named tuple (the params tuple) defined in SelectSolver.get_params(...).
    #           This tuple contains the parameters defined in the parameters file params.jl.
    #           For more information about the SA parameters, consult the SOSS manual or the
    #           the params.jl file provided in the WORK/INPUTS/TEST_1 directory.
    #
    # Returns:
    # solution -> Tuple containing the solution information.
    #             The solution information consists on the vector containing the ID number of each
    #             run, the vector containing the final energy of each run, and the vector of vectors
    #             containing the final configuration of ions for each run.

    # load solver parameters
    @unpack IT_scheme, IT_params, tempLength, beta, scheme, gamma,
            ScanningStyle, FirstToImprove, FirstToImprove, steps,
            MaxStagnation, nRepeats, DisplayStyle, plots = params

    # load system data
    @unpack IonSiteConstraint, Ions, LIons, IonsSubset, NoptSites, Na, optSites,
            SitesSubset , U, UionAionBlist, EnergyBase   = data

    # check if there is a seeds file or not
    if isfile("seeds.dat")
        seeds = readdlm("seeds.dat", Int64)
        seeds_configs = seeds[1:nRepeats,1] # read initial configuration rng seeds from file
        seeds_optruns = seeds[1:nRepeats,2] # read solver rng seeds from file
    else
        seeds_configs = rand(1:10000000, nRepeats) # generate initial configuration rng at random
        seeds_optruns = rand(1:10000000, nRepeats) # generate solver rng seeds at random

        # write seeds.jl file
        seedsdata = hcat(seeds_configs, seeds_optruns)
        writedlm("seeds.dat", seedsdata)
    end

    # remove (if necessary) and make an optimization files directory
    OptDIR = "OPTFILES"
    if isdir(OptDIR)
        rm(OptDIR, recursive=true)
    end
    mkdir(OptDIR)

    # get the total number of ions
    L = sum(LIons)

    # initialize general variables
    Tmax       = 0.0
    Tmin       = 0.0
    accepted   = 0
    nonzero    = zeros(Bool, nRepeats)
    runtime    = 0.0
    RunList    = zeros(Int64, nRepeats)
    TimeList   = zeros(Float64, nRepeats)
    EnergyList = zeros(Float64, nRepeats)
    ConfigList = [zeros(Int64, L) for i in 1:nRepeats]

    # get the initial/final (max/min) temperature
    Tmax, Tmin = InitialTemperature(IT_scheme, IT_params, data, IonSiteConstraint)

    # loop of runs (in parallel, but not distributed)
    println("START RUNS")
    Threads.@threads for run in 1:nRepeats

        # make run directory
        RunDir = "$OptDIR/RUN_$run"
        if steps > 0
            mkdir(RunDir)
        end

        # initialize rng
        rng_configs = MersenneTwister(seeds_configs[run])
        rng_optruns = MersenneTwister(seeds_optruns[run])

        # generate initial solution (configuration/distribution of ions) and its neighborhood
        if IonSiteConstraint == false # without ion-site constraint
            config = get_init_sol(LIons, NoptSites, optSites, rng_configs) # initial configuration
            nbrs, MoveStyle, Nnbrs = get_nbrs(config, L, NoptSites) # neighborhood
        else # with ion-site constraint
            config = get_init_sol_constrained(LIons, NoptSites, optSites, SitesSubset, IonsSubset,
                                              rng_configs) # initial configuration
            nbrs, MoveStyle, Nnbrs = get_nbrs_constrained(config, SitesSubset,
                                                          IonsSubset, NoptSites) # neighborhood
        end

        # run optimization
        if steps > 0
            # search optimal optimization
            stats = @timed config, energy = SANSS_opt(L, Na, gamma, ScanningStyle,
                                                     FirstToImprove, optSites,
                                                     IonsSubset, U, UionAionBlist, EnergyBase, 
                                                     Tmax, Tmin, beta, tempLength,
                                                     scheme, steps, MaxStagnation, config,
                                                     nbrs, MoveStyle, Nnbrs, run, RunDir,
                                                     DisplayStyle, rng_optruns)

            # store final configuration
            accepted += 1
            RunList[run]    = run
            TimeList[run]   = stats.time
            EnergyList[run] = energy
            ConfigList[run] = config[3,:]
        else # steps = 0 -> random (initial) configuration
            energy = get_TotalEnergy(config, L, U, UionAionBlist, EnergyBase)
            println(@sprintf("RUN %03d: RANDOM CONFIGURATION GENERATED", run))
            accepted += 1
            RunList[run]    = run
            EnergyList[run] = energy
            ConfigList[run] = config[3,:]
        end
    end

    # nonzero elements (accepted final configurations)
    nonzero = RunList .!= 0
    # check nonzero elements
    if length(nonzero) == 0
        println("None of the final configurations were accepted.
        Check the optimization settings and solver parameters used.)")
    end
    # execution time
    if steps > 0
        runtime = mean(TimeList[nonzero])
    else
        runtime = 0.0
    end

    # plots (optional)
    if plots && steps > 0
        for run in 1:nRepeats
            RunDir = "$OptDIR/RUN_$run"
            plot_E(RunDir)
        end
    end

    # write optimization summary
    write_opt(
        nRepeats, tempLength, scheme, IT_scheme, IT_params, steps, beta,
        gamma, FirstToImprove, ScanningStyle, MaxStagnation, accepted,
        EnergyList[nonzero]/Na, RunList[nonzero], runtime, OptDIR, DisplayStyle, plots
    )

    # solution tuple
    soltution = (
        runs = RunList[nonzero],               # accepted runs
        energies = EnergyList[nonzero],        # energies of the accepted runs
        configs = ConfigList[nonzero]          # ion distributions of the accepted runs
    )

    return soltution
end

####################################################################################################

end