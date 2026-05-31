# load modules
include("../../../SRC/STRUCTURE/ImportSettings.jl")
include("../../../SRC/STRUCTURE/ImportData.jl")
include("../../../SRC/STRUCTURE/ExportData.jl")
include("../../../SRC/OPTIMIZATION/SelectSolver.jl")
using .ImportSettings
using .ImportData
using .ExportData
using .SelectSolver

# load optimization settings
Base.include(ImportSettings, "settings.jl")
# load optimization solver parameters
Base.include(SelectSolver, "params.jl")

# import settings
settings = ImportSettings.get_settings()
# import data
t1 = time()
data = ImportData.get_data(settings)
t2 = time()
println("Import data: $(t2 - t1) seconds")
# solve the problem
t3 = time()
solution = SelectSolver.solve(data, settings)
t4 = time()
println("Optimization: $(t4 - t3) seconds")
# export results
ExportData.save_structure(data, solution, settings)