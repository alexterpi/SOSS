#===================================================
| SIMULATED ANNEALING STRUCTURE SEARCH PARAMETERS: |
| 1. DO NOT MODDIFY THE VARIABLE TYPES             |
| 2. FIND INFORMATION ABOUT THE PARAMETERS BELOW   |
===================================================#

# Initial Temperature
IT_scheme::Int64 = 1
IT_params::Vector{Float64} = [1000, 10]

# Temperature Length
tempLength::Int64 = 10_000

# Main Cooling Scheme
# N/A

# Secondary Cooling Scheme
beta::Float64 = 0.1
scheme::String = "linear"

# Exploration Criterion
gamma::Int64 = 100
ScanningStyle::Int64 = 1
FirstToImprove::Bool = true

# Stopping Criterion
steps::Int64 = 10
MaxStagnation::Int64 = 0

# Recurrence
nRepeats::Int64 = 1

# Display
DisplayStyle::Int64 = 2

# Plots
plots::Bool = true

#===================================================
| Note:                                            |
|  To reproduce some specific results, a RNG seeds |
| file (seeds.dat) can be provided.                |
|  The file must contain (at least) "nRepeats"     |
| double-column rows. Each column must be an       |
| integer, and a space must separate the columns.  |
|  The two numbers in each row are used as seeds   |
| to initialize the initial solution RNG (first    |
| number) and the optimization solver RNG (second  |
| number).                                         |
|  If the file contains more than "nRepeats" rows, |
| only the first "nRepeats" rows are used, and the |
| rest are ignored. If the file contains fewer     |
| than "nRepeats" rows, an error occurs.           |
===================================================#
#===================================================
|--------------------------------------------------|
| nRepeats                                         |
| def: Number of runs (optimization cycles).       |
|..................................................|
| steps                                            |
| def: Number of steps in the Main Cooling Scheme. |
| typical values: 10 - 1000                        |
|..................................................|
| IT_scheme                                        |
| def: Initial (and final) Temperature scheme.     |
| typical values: 1 -> The user selects the        |
|                      max and min temperature     |
|                      (OIT1)                      |
|                 2 -> The max/min temperature is  |
|                      proportional to the max/min |
|                      difference of energy        |
|                      obtained during two         |
|                      consecutive steps of a      |
|                      random walk performed with  |
|                      the target system (OIT2)    |
|                 3 -> The max/min temperature is  |
|                      selected to match some      |
|                      max/min probability of      |
|                      acceptance, defined as      |
|                      p = exp(-T/ΔE), where ΔE is |
|                      the average difference of   |
|                      energy during two           |
|                      consecutive steps of a      |
|                      random walk performed with  |
|                      the target system (OIT3)    |
| recommended: IT_scheme = 1                       |
| info: The max/min temperature is the initial     |
|       temperature in the initial/final step of   |
|       the Main Cooling Scheme.                   |
|..................................................|
| IT_params                                        |
| def: Parameters asociated to each IT scheme.     |
|                                                  |
| if IT_scheme = 1 -> IT_params[1] = Tmax          |
|                     IT_params[2] = Tmin          |
|                                                  |
| typical values: IT_params[1] = 1000              |
|                 IT_params[2] = 10                |
|                                                  |
| if IT_scheme = 2 -> IT_params[1] = random walk   |
|                                    steps         |
|                     IT_params[2] = scale         |
|                                    constant      |
|                                                  |
| if IT_scheme = 3 -> IT_params[1] = random walk   |
|                                    steps         |
|                     IT_params[2] = max prob.     |
|                     IT_params[3] = min prob.     |
|..................................................|
| tempLength                                       |
| def: Number of steps in the Secondary Cooling    |
|      Scheme.                                     |
| typical values: 1000 - 10000                     |
|..................................................|
| beta                                             |
| def: Intensity of the secondary cooling.         |
|      T_final = beta · T_initial                  |
| default value: 0.1                               |
| info: Ignored when using OSCS1.                  |
|..................................................|
| scheme                                           |
| def: Secondary Cooling Scheme.                   |
| possible values: constant (OSCS1), linear (OSCS2)|
| recommended: linear                              |
|..................................................|
| gamma                                            |
| def: Maximum number of explored solutions before |
|      selecting the candidate solution (ONE1).    |
| typical values: 100 - 1000                       |
|..................................................|
| ScanningStyle                                    |
| def: Update of the number of explored solutions  |
|      (ONE1) during the optimization              |
| possible values: 1 (no update)                   |
|..................................................|
| FirstToImprove                                   |
| def: Parameter to activate ONE2 (choose the first|
|      improving solution during the Neighbourhood |
|      Exploration).                               |
| possible values: true  -> activate ONE2          |
|                  false -> continue with ONE1     |
|                           (explore gamma         |
|                           solutions and choose   |
|                           the best one)          |
|..................................................|
| MaxStagnation                                    |
| def: Limit of stagnation.                        |
| if MaxStagnation > 0 (max stagnation value)      |
| if MaxStagnation = 0 (stagnation is ignored)     |
| recommended value: 0                             |
| note: At every main cooling step the improvement |
|       of the energy is checked, if there is not  |
|       improvement the stagnation is increased by |
|       1.                                         |
|       The optimization stops when MaxStagnation  |
|       is reached.                                |
|..................................................|
| DisplayStyle                                     |
| def: Display style.                              |
| typical values: 0 -> minimalist                  |
|                 1 -> optimal configurations      |
|                 2 -> steps                       |
| note: When running multiple optimizations in     |
|       parallel, we recommend to set it to 0.     |
|..................................................|
| plots                                            |
| def: Whether or not to plot the data in the      |
|      files generated during the optimization.    |
| typical values: true  -> plot                    |
|                 false -> do not plot             |
| note: When running multiple optimizations in     |
|       in parallel, plotting the data can produce |
|       errors in some environments, so we         |
|       recommend to set it to false (and do the   |
|       plots after the optimization).             |
|--------------------------------------------------|
===================================================#