#=================================================
| HARMONY SEARCH ALGORITHM PARAMETERS:           |
| 1. DO NOT MODDIFY THE VARIABLE TYPES           |
| 2. FIND INFORMATION ABOUT THE PARAMETERS BELOW |
=================================================#
nruns::Int64        = 1 # 10
lastN::Int64        = 20_000 # 20_000
nIterations::Int64  = 100_000 # 100_000
hcmr::Float64       = 0.9
par::Float64        = 0.7
hms_0::Int64        = 10
hms_diff::Int64     = 5
hms_steps::Int64    = 5 # 9
display_mode::Int64 = 1
#================================================= 
|------------------------------------------------|
| nruns                                          |
| def: number of runs for each value of hms      |
|................................................|
| lastN                                          |
| def: number of last configurations considered  |
|      in single-run statistics                  |
|................................................|
| nIterations                                    |
| def: number of iterations for each run         |
| typical values: 30k - 100k                     |
|................................................|
| hcmr                                           |
| def: Harmony Memory Considering Rate           |
| typical values: 0.70 - 0.95                    |
|................................................|
| par                                            |
| def: Pitching Adjust Rate                      |
| typical values: 0.20 - 0.50                    |
|................................................|
| hms (hms_0, hms_diff, hms_steps)               |
| hms[i] = hms_0 + (i-1) * hms_diff              |
| for i = 1,...,hms_steps                        |
| def: Harmony Matrix Size                       |
| typical values: 10 - 50                        |
|................................................|
| display                                        |
| def: show information by screen                |
| 0 -> off                                       |
| 1 -> on                                        |
|------------------------------------------------|
=================================================#