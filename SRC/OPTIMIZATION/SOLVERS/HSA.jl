module HSA

#using PyPlot
using PythonCall
using PythonPlot
using StatsBase
using Statistics
using Random
using DelimitedFiles
using Test
using UnPack
using MutableNamedTuples

# Define references to Python modules/functions
const plt = Ref{Py}(Py(nothing))

# Initialize Python environment
function __init__()
    plt[] = pyimport("matplotlib.pyplot")
end

function get_ilEnergy_ion_atoms(
    LIons::Vector{Int64},
    il::Array{Int64,1},
    U::Array{Float64,2}
    )

    # total number of ions
    L = sum(Int64, LIons)

    # sum all ion-atoms contributions
    ion = 1; count = 0
    s = 0.0
    for i in 1:L
        count += 1
        if count > LIons[ion]
            ion += 1
            count = 1 
        end
        s += U[il[i], ion]
    end

    return s
end

function get_ilEnergy_ion_ion(
    LIons::Vector{Int64},
    il::Array{Int64,1},
    UionAionBlist::Array{Float64,4},
    )

    L = sum(Int64, LIons)
    ion_i = 1; count_i = 0
    Uion = 0.0
    for i in 1:(L-1)
        count_i += 1
        if count_i > LIons[ion_i]
            ion_i += 1
            count_i = 1
        end
        if (count_i + 1) > LIons[ion_i]
            ion_j = ion_i + 1; count_j = 0
        else
            ion_j = ion_i; count_j = count_i
        end
        for j in (i+1):L
            count_j += 1
            if count_j > LIons[ion_j]
                ion_j += 1
                count_j = 1
            end
            Uion += UionAionBlist[il[i], il[j], ion_i, ion_j]
        end
    end

    return Uion
end

function getInitMemory(
    LIons::Vector{Int64},
    Ne::Int,
    Nmem::Int,
    Nv_list::Array{Int64,1},
    il::Array{Int64,1},
    U::Array{Float64,2},
    UionAionBlist::Array{Float64,4},
    EnergyBase::Float64
    )
    # total ions
    L  = sum(Int64, LIons)
    # harmony memory
    hM_il = zeros(Int64, Nmem, L) # 4×L×nWalkers ::Array{Float64,3}: 
    hM_ev = zeros(Int64, Nmem, Ne) # 4×L×nWalkers ::Array{Float64,3}: 
    hM_U  = zeros(Float64, Nmem)
    worstMemIndex = -1 # int that can take a value from 1 to Nmem
    worstEnergy = 0.0

    for memRow in 1:Nmem
        get_il_proposal!(L, Nv_list, il)
        copy_il_to_hM!(L, il, hM_il, memRow)

        # determine energy of il
        u = get_il_Ut_plus_Uion(LIons, il, U, UionAionBlist, EnergyBase)
        
        hM_U[memRow] = u
        if isWorstEnergy(memRow, u, worstEnergy)
            worstMemIndex = memRow
            worstEnergy   = hM_U[memRow]
        end
    end
    #
    return hM_il, hM_U, worstMemIndex, worstEnergy
end

function get_il_proposal!(
    L::Int,
    Nv_list::Array{Int64,1},
    il::Array{Int64,1}
    )
    shuffle!(Nv_list)
    for i in 1:L
        il[i] = Nv_list[i]
    end
end

function vacancy_to_site(il::Vector{Int64}, removedSites::Vector{Int64})

    sites = zeros(Int64, length(il))
    for (i, vacancy) in enumerate(il)
        sites[i] = removedSites[vacancy]
    end

    return sites
end

function copy_il_to_hM!(
    L::Int,
    il::Array{Int64,1},
    hM_il::Array{Int64,2},
    memRow::Int
    )
    for i in 1:L
        hM_il[memRow, i] = il[i]
    end
end

function isWorstEnergy(
    memRow::Int,
    e::Float64,
    worstEnergy::Float64
    )
    return (memRow > 1) ? (e > worstEnergy) : true
end

function wrap( x::Int, Nv::Int )
    if x > Nv
        return x - Nv
    elseif x < 1
        return x + Nv
    else
        return x
    end
end

function getPitchedValues(
    ia::Int,
    delta::Int,
    Nv::Int
    )

    if rand() < 0.5
        ia1 = wrap( ia - delta, Nv)
        ia2 = wrap( ia + delta, Nv)
    else
        ia1 = wrap( ia + delta, Nv)
        ia2 = wrap( ia - delta, Nv)
    end
    return ia1, ia2
end

function getShuffled_MemColumn!(
    Nmem::Int,
    col::Int,
    memList::Array{Int,1},
    hM_il::Array{Int,2},
    )

    for row in 1:Nmem
        memList[row] = hM_il[row, col]
    end
    shuffle!(memList)
end

function get_il_Ut_plus_Uion(
    LIons::Vector{Int64},
    il::Array{Int64,1},
    U::Array{Float64,2},
    UionAionBlist::Array{Float64,4},
    energyBase::Float64
    )

    # determine energy of il (ion-atom and ion-ion energy)
    Ut   = get_ilEnergy_ion_atoms(LIons, il, U)
    Uion = get_ilEnergy_ion_ion(LIons, il, UionAionBlist)

    # total energy = atom-atom energy + ion-atom energy + ion-ion energy
    Etot = energyBase + Ut + Uion

    return Etot
end

function compute_statistics(
    ndata::Int64,
    bestE::Vector{Float64},
    worstE::Vector{Float64},
    acceptance::Vector{Bool}
    )
    # best energy
    data = bestE[end - ndata + 1 : end]
    bestE_mean = mean(data)
    bestE_std  = std(data, corrected=true, mean=bestE_mean)
    # worst energy
    data = worstE[end - ndata + 1 : end]
    worstE_mean = mean(data)
    worstE_std  = std(data, corrected=true, mean=worstE_mean)
    # acceptance
    data = acceptance[end - ndata + 1 : end]
    acceptance_mean = mean(data)
    return bestE_mean, bestE_std, worstE_mean, worstE_std, acceptance_mean
end

function plot1(
    namePlot1::String,
    nIterations::Int,
    bestE_hist::Array{Float64,1},
    worstE_hist::Array{Float64,1},
    par_hist::Array{Float64,1},
    hcmr_hist::Array{Float64,1},
    hms::Int64
    )

    x        = collect(1:nIterations)
    x1_ymax  = maximum(worstE_hist)
    x1_ymin  = minimum(bestE_hist) - 0.01*abs(minimum(bestE_hist))

    plt[].close()
    plt[].ioff()
    fig, ax1 = plt[].subplots()
    ax1.set_xlabel("ITERATION")
    ax1.set_ylabel("ENERGY (eV/ATOM)", color="k")
    ax1.set_ylim([x1_ymin, x1_ymax])
    ax1.plot(x, bestE_hist, "-r", label = "BEST E")
    ax1.plot(x, worstE_hist,"-.k", label = "WORST E")
    plt[].legend(loc = "center left", framealpha = 1, edgecolor = "k")
    ax2 = ax1.twinx()
    ax2.set_ylabel("PROBABILITY", color="b")
    ax2.set_ylim([0.0, 1.0])
    ax2.tick_params(axis="y", colors="b")
    ax2.plot(x, par_hist, "--b", linewidth=1.4, label = "PAR")
    ax2.plot(x, hcmr_hist, ":g", linewidth=1.4, label = "HCMR")
    #fig.legend(loc = "upper center")
    plt[].legend(loc = "center right", framealpha = 1, edgecolor = "b")
    plt[].title("HMS = $hms")
    fig.savefig(namePlot1)
    plt[].close(fig)
end

function plot2(
    x_vals::Array{Int64,1},
    y1mean::Array{Float64,1},
    y1dev::Array{Float64,1},
    y2mean::Array{Float64,1},
    y2dev::Array{Float64,1},
    y3mean::Array{Float64,1},
    y3dev::Array{Float64,1}
    )

    ax1_ymin = minimum(y1mean) - maximum(y1dev) - 0.01*(minimum(y1mean) - maximum(y1dev))
    ax1_ymax = maximum(y2mean) + maximum(y2dev) + 0.01*(maximum(y2mean) + maximum(y2dev))
    ax2_ymax = maximum(y3mean) + maximum(y3dev) + 0.01*(maximum(y3mean) + maximum(y3dev))

    plt[].close()
    plt[].ioff()
    fig, ax1 = plt[].subplots()
    ax1.set_xlabel("HMS")
    ax1.set_ylabel("mean energy (eV/ATOM)", color="k")
    ax1.set_ylim([ax1_ymin, ax1_ymax])
    plt[].errorbar(x_vals, y2mean, yerr=y2dev, color="k", label = "worst E")
    plt[].errorbar(x_vals, y1mean, yerr=y1dev, color="r", label = "best E")
    plt[].legend(loc = "upper left", framealpha = 1)
    ax2 = ax1.twinx()
    ax2.set_ylabel("time (seconds)", color="b")
    ax2.set_ylim([0.0, ax2_ymax])
    ax2.tick_params(axis="y", colors="blue")
    plt[].errorbar(x_vals, y3mean, yerr=y3dev, color="b", label = "time")
    plt[].legend(loc = "upper right", framealpha = 1)
    fig.savefig("OPTFILES/hsa_evolution.png")
    plt[].close(fig)
end

function writedata1(
    file::String,
    hms::Int64,
    E_opt::Float64,
    ndata::Int64,
    mean_bestE::Float64,
    std_bestE::Float64,
    mean_worstE::Float64,
    std_worstE::Float64,
    acceptance::Float64
    )
    f = open(file, "a")

    write(f, "HMS = $hms\n")
    write(f, "E_OPT = $E_opt eV/ATOM\n")
    write(f, "LAST $ndata CONFIGURATIONS\n")
    write(f, "E_BEST  = $mean_bestE ± $std_bestE eV/ATOM\n")
    write(f, "E_WORST = $mean_worstE ± $std_worstE eV/ATOM\n")
    write(f, "ACCEPTANCE = $acceptance\n")
    write(f, "\n")
    close(f)

end

function writedata2(
    nruns::Int64,
    niter::Int64,
    hcmr::Float64,
    par::Float64,
    x_vals::Array{Int64,1},
    y1mean::Array{Float64,1},
    y1dev::Array{Float64,1},
    y2mean::Array{Float64,1},
    y2dev::Array{Float64,1},
    y3mean::Array{Float64,1},
    y3dev::Array{Float64,1},
    RunList::Vector{Int64},
    HmsList::Vector{Int64},
    EnergyList::Vector{Float64}
    )

    file = open("OPTFILES/OPTFILE.dat", "w")

    # parameters
    write(file, "HSA\n")
    write(file, "\n")
    write(file, "PARAMETERS:\n")
    write(file, "NUMBER OF RUNS = $nruns\n")
    write(file, "NUMBER OF ITERATIONS = $niter\n")
    write(file, "HMS = ")
    for hms in x_vals
        write(file, "$hms ")
    end
    write(file, "\n")
    write(file, "HCMR = $hcmr\n")
    write(file, "PAR = $par\n")
    write(file, "\n")

    # best configurations
    opt_run = argmin(EnergyList)
    write(file, "BEST CONFIGURATIONS\n")
    for i = 1:nruns
        write(file, "RUN #$(RunList[i]): HMS = $(HmsList[i]), E = $(EnergyList[i]) eV/ATOM")
        if i == opt_run
            write(file, " [OPTIMAL CONFIGURATION]")
        end
        write(file, "\n")
    end
    write(file, "\n")

    # average over N last configurations
    write(file, "STATISTICS\n")
    for i in 1:length(x_vals)
        write(file, "STEP $i:\n")
        write(file, "HMS = $(x_vals[i])\n")
        write(file, "E_BEST  = $(y1mean[i]) ± $(y1dev[i]) eV/ATOM\n")
        write(file, "E_WORST = $(y2mean[i]) ± $(y2dev[i]) eV/ATOM\n")
        write(file, "TIME    = $(y3mean[i]) ± $(y3dev[i])SECONDS/RUN\n")
    end
    close(file)
end

function displaysummary1(
    run::Int64,
    hms::Int64,
    E_opt::Float64,
    ndata::Int64,
    mean_bestE::Float64,
    std_bestE::Float64,
    mean_worstE::Float64,
    std_worstE::Float64,
    acceptance::Float64
    )
    println("RUN #$run:")
    println("HMS = $hms")
    println("E_OPT = $E_opt eV/ATOM")
    println("LAST $ndata CONFIGURATIONS")
    println("E_BEST  = $mean_bestE ± $std_bestE eV/ATOM")
    println("E_WORST = $mean_worstE ± $std_worstE eV/ATOM")
    println("ACCEPTANCE = $acceptance")
end

function displaysummary2(
    nruns::Int64,
    RunList::Vector{Int64},
    HmsList::Vector{Int64},
    EnergyList::Vector{Float64},
    x_vals::Vector{Int64},
    bestE_mean::Vector{Float64},
    bestE_std::Vector{Float64},
    worstE_mean::Vector{Float64},
    worstE_std::Vector{Float64},
    time_mean::Vector{Float64},
    time_std::Vector{Float64}
    )
    println("----------------------------------------------------------------------------------")
    println("SUMMARY")
    println("----------------------------------------------------------------------------------")

    println("AFTER $nruns RUNS:")
    println("BEST CONFIGURATIONS")
    for i in 1:nruns
        println("RUN #$(RunList[i]): HMS = $(HmsList[i]), E = $(EnergyList[i]) eV/ATOM")
    end
    println("STATISTICS")
    for (k, hms) in enumerate(x_vals)
        println("STEP $k:")
        println("HMS     = $hms")
        println("E_BEST  = $(bestE_mean[k]) ± $(bestE_std[k]) eV/ATOM")
        println("E_WORST = $(worstE_mean[k]) ± $(worstE_std[k]) eV/ATOM")
        println("TIME    = $(time_mean[k]) ± $(time_std[k]) SECONDS/RUN")
    end
end

function newHarmony!(    
    L::Int,
    Nv::Int,
    Nmem::Int,
    hcmr::Float64,
    par::Float64,
    il::Array{Int,1},
    Nv_list::Array{Int,1},
    memoryList::Array{Int,1},
    hM_il::Array{Int,2},
    )

    boolNotRefilled = ones(Bool, Nv) # = [true, true, true, ...]
    shuffle!(Nv_list)

    for col in 1:L
        a_final = -1
        if rand() <= hcmr
            # memList will be a column of hM_il, shuffled            
            getShuffled_MemColumn!(Nmem, col, memoryList, hM_il) # mutates memList

            pitch = (rand() <= par)

            if !pitch
                repeat = true
                count = 0
                while repeat && (count < Nmem)
                    count += 1
                    a = memoryList[count]
                    if boolNotRefilled[a]
                        a_final = a
                        repeat = false
                    end
                end
                pitch = repeat
            end
            #
            if pitch
                a = memoryList[1]
                repeat = true
                delta = 0
                while repeat && (delta < Nv)
                    delta += 1
                    a1, a2 = getPitchedValues(a, delta, Nv)
                    if boolNotRefilled[a1]
                        a_final = a1
                        repeat = false
                    elseif boolNotRefilled[a2]
                        a_final = a2
                        repeat = false
                    end
                end
            end
        else
            count = 0
            repeat = true
            while repeat && (count < Nv)
                count += 1
                a = Nv_list[count]
                if boolNotRefilled[a]
                    a_final = a
                    repeat = false
                end
            end
        end
        il[col] = a_final
        boolNotRefilled[a_final] = false
    end
end

function loop_par!(
    LIons::Vector{Int64},
    Nv::Int,
    Ne::Int,
    Nmem::Int,
    Nv_list::Array{Int,1},
    il::Array{Int,1},
    U::Array{Float64,2},
    UionAionBlist::Array{Float64,4},
    hcmr::Float64,
    par::Float64,
    tempListNmem::Array{Int,1},
    nIterations::Int,
    EnergyBase::Float64,
    bestEnergyMem::Array{Float64,1},
    worstEnergyMem::Array{Float64,1},
    lpar::Array{Float64,1},
    lhcmr::Array{Float64,1},
    acceptance_hist::Array{Bool,1},
    optimal::MutableNamedTuple,
    removedSites::Vector{Int64}
    )

    # total number of ions
    L  = sum(Int64, LIons)

    # initial random configurations and worst configurartion
    hM_il, hM_U, worstRow, worstEnergy = getInitMemory(LIons, Ne, Nmem, Nv_list, il, U, UionAionBlist, EnergyBase)
    
    # reset optimal energy
    optimal.energy = EnergyBase

    for i in 1:nIterations

        lpar[i]  = par
        lhcmr[i] = hcmr

        # new configurartion (il)
        newHarmony!(L, Nv, Nmem, hcmr, par, il, Nv_list, tempListNmem, hM_il)
    
        # determine energy of il
        u = get_il_Ut_plus_Uion(LIons, il, U, UionAionBlist, EnergyBase)

        # if il's energy is better than the worst row in memory, replace it
        if u < worstEnergy
            copy_il_to_hM!(L, il, hM_il, worstRow)
            hM_U[worstRow] = u
            acceptance_hist[i] = true
        else
            acceptance_hist[i] = false
        end
        
        # determine again the new worstRow and worstEnergy
        worstRow          = argmax(hM_U)
        worstEnergy       = hM_U[worstRow]
        worstEnergyMem[i] = worstEnergy
        bestEnergyMem[i]  = hM_U[argmin(hM_U)]

        # check if last configuration is the best configuration
        if u < optimal.energy
            il_sites = vacancy_to_site(il, removedSites)
            optimal.energy = u
            optimal.config = il_sites
        end

    end
end

function optimize_system(data::NamedTuple, params::NamedTuple)

    # solver parameters
    @unpack nruns = nruns, lastN, nIterations, hcmr, par, hms_0, hms_diff, hms_steps, display_mode = params
    # structure data
    @unpack Ions, LIons, optSites, NoptSites, Nv, Na, U, UionAionBlist, EnergyBase = data

    # New terminology, update: removedSites -> optSites, Nv -> NoptSites, Ne -> Nv (and associated variables)
    removedSites = optSites
    Ne = Nv
    Nv = NoptSites
    #------------------------#

    L = sum(Int64, LIons)

    # optimization main directory
    if isdir("OPTFILES")
        rm("OPTFILES", recursive=true)
    end
    mkdir("OPTFILES")

    # optimization files and subdirectories
    for i in 1:nruns
        if !isdir("OPTFILES/RUN$i")
            mkdir("OPTFILES/RUN$i")
        end
    end

    # hms loop steps
    x_vals = [hms_0 + hms_diff * (i-1) for i = 1:hms_steps]

    # vacancies list
    Nv_list::Array{Int64,1} = collect(1:Nv)

    # configuration
    il::Array{Int64,1} = zeros(Int64, L)

    # historial
    bestE_hist::Array{Float64,1}   = zeros(Float64, nIterations)
    worstE_hist::Array{Float64,1}  = zeros(Float64, nIterations)
    par_hist::Array{Float64,1}     = zeros(Float64, nIterations)
    hcmr_hist::Array{Float64,1}    = zeros(Float64, nIterations)
    acceptance_hist::Array{Bool,1} = zeros(Bool   , nIterations)

    # statistics over last configurations
    bestE_mean      = zeros(Float64, lastN)
    bestE_std       = zeros(Float64, lastN)
    worstE_mean     = zeros(Float64, lastN)
    worstE_std      = zeros(Float64, lastN)
    acceptance_mean = zeros(Float64, lastN)

    # statistics over runs
    yBestMean::Array{Float64,1}  = zeros(Float64, hms_steps)
    yBestDev::Array{Float64,1}   = zeros(Float64, hms_steps)
    yWorstMean::Array{Float64,1} = zeros(Float64, hms_steps)
    yWorstDev::Array{Float64,1}  = zeros(Float64, hms_steps)
    yTimeMean::Array{Float64,1}  = zeros(Float64, hms_steps)
    yTimeDev::Array{Float64,1}   = zeros(Float64, hms_steps)

    # optimal configuration last run
    optimal = MutableNamedTuple(
        energy  = EnergyBase::Float64,
        config  = zeros(Int64, L)
    )

    # solutions for each run
    RunList    = collect(1:nruns)
    HmsList    = zeros(Int64, nruns)
    EnergyList = EnergyBase * ones(Float64, nruns)
    ConfigList = [Vector{Int64}(undef, L) for _ in 1:nruns]

    # start the optimization
    if nruns > 0
        for (k, Nmem) in enumerate(x_vals) # start HMS loop
            println("----------------------------------------------------------------------------------")
            println("STEP $k: HMS = $Nmem")
            println("----------------------------------------------------------------------------------")
            tempListNmem = zeros(Int, Nmem)
            results = zeros(3, nruns)
            for i in 1:nruns # start runs loop
                
                # compute new configuration
                time = @elapsed loop_par!(LIons, Nv, Ne, Nmem, Nv_list, il, U,
                                          UionAionBlist, hcmr, par, tempListNmem, nIterations,
                                          EnergyBase, bestE_hist, worstE_hist, par_hist, hcmr_hist,
                                          acceptance_hist, optimal, removedSites)

                # update best solution
                if optimal.energy < EnergyList[i]
                    EnergyList[i] = optimal.energy
                    ConfigList[i] = optimal.config
                    HmsList[i]    = Nmem
                end

                # save historial for statistics over runs
                results[1, i] = bestE_hist[end]
                results[2, i] = worstE_hist[end]
                results[3, i] = time
                
                # statistics over last configurations
                bestE_mean , bestE_std ,
                worstE_mean, worstE_std,
                acceptance_mean          = compute_statistics(lastN, bestE_hist, worstE_hist, acceptance_hist)

                # write/plot results
                file = "OPTFILES/RUN$i/SUMMARY.dat"
                fig  = "OPTFILES/RUN$i/HSM$Nmem.png"
                writedata1(file, Nmem, optimal.energy/Na, lastN,
                           bestE_mean/Na, bestE_std/Na,
                           worstE_mean/Na, worstE_std/Na,
                           acceptance_mean)
                plot1(fig, nIterations, bestE_hist/Na, worstE_hist/Na, par_hist, hcmr_hist, Nmem)

                # show results on screen
                if display_mode == 1
                    displaysummary1(i, Nmem, optimal.energy/Na, lastN, bestE_mean/Na, bestE_std/Na,
                                    worstE_mean/Na, worstE_std/Na, acceptance_mean)
                    println("")
                end
                
            end # end runs loop

            # statistics over runs
            yBestMean[k]  = mean( results[1, 1:nruns] )
            yBestDev[k]   = std(  results[1, 1:nruns] )
            yWorstMean[k] = mean( results[2, 1:nruns] )
            yWorstDev[k]  = std(  results[2, 1:nruns] )
            yTimeMean[k]  = mean( results[3, 1:nruns] )
            yTimeDev[k]   = std(  results[3, 1:nruns] )

        end # end HMS loop

        # write/plot results
        writedata2(nruns, nIterations,hcmr, par,
                   x_vals, yBestMean/Na, yBestDev/Na, yWorstMean/Na, yWorstDev/Na,
                   yTimeMean, yTimeDev, RunList, HmsList, EnergyList/Na)
        if hms_steps > 1 & nruns > 1
            plot2(x_vals, yBestMean/Na, yBestDev/Na, yWorstMean/Na, yWorstDev/Na,
                  yTimeMean, yTimeDev)
        end

        # show results on screen
        if display_mode == 1
            displaysummary2(nruns, RunList, HmsList, EnergyList/Na, x_vals,
                            yBestMean/Na, yBestDev/Na, yWorstMean/Na, yWorstDev/Na,
                            yTimeMean, yTimeDev)
        end

    else
        Nmem = 10
        tempListNmem = zeros(Int, Nmem)
        time = @elapsed loop_par!(LIons, Nv, Ne, Nmem, Nv_list, il, U, UionAionBlist, hcmr, par, tempListNmem, nIterations, EnergyBase, bestE_hist, worstE_hist, par_hist, hcmr_hist, acceptance_hist, optimal, removedSites)
        plot1("hsa_", nIterations, bestE_hist, worstE_hist, par_hist, hcmr_hist, Nmem)
    end

    # final solutions
    solution = (
        runs = RunList,
        energies = EnergyList,
        configs = ConfigList
    )

    return solution
end

end