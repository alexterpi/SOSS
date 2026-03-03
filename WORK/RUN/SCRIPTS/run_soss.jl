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
data = ImportData.get_data(settings)
# solve the problem
solution = SelectSolver.solve(data, settings)
# export results
ExportData.save_structure(data, solution, settings)