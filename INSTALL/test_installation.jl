# ==========================================================
# Test Installation Script
# ==========================================================
using TOML
using Pkg
using Printf

println("---------------------------------------------------------")
println("Starting SOSS Full Installation Test")
println("---------------------------------------------------------")

# --- PART 1: Validate Julia Packages (Project.toml) ---
println("Step 1: Checking Julia dependencies...")

try
    # Parse the Project.toml file to identify required dependencies
    project_data = TOML.parsefile(joinpath(pwd(), "Project.toml"))
    julia_deps = collect(keys(project_data["deps"]))
    
    println("   Detected $(length(julia_deps)) Julia packages.")
    
    for pkg in julia_deps
        # Dynamically evaluate 'using' for each package to verify availability
        eval(Meta.parse("using $pkg"))
        println("   Julia package: $pkg is ready.")
    end
    println("All Julia packages loaded successfully.\n")
catch e
    # Updated error message to focus on system requirements and prerequisites
    println("\nERROR during installation:")
    println(e)
    println("\nPlease verify that the system meets the installation requirements.")
    exit(1)
end

# --- PART 2: Validate Python Environment (CondaPkg.toml) ---
println("Step 2: Checking Python environment...")

# Import interoperability tools already validated in the previous step
using PythonCall

try
    toml_path = joinpath(pwd(), "CondaPkg.toml")
    if !isfile(toml_path)
        println("ERROR: 'CondaPkg.toml' not found.")
        exit(1)
    end

    # Parse Python dependencies defined in CondaPkg.toml
    config = TOML.parsefile(toml_path)
    pkg_dict = config["deps"]
    pkg_names = collect(keys(pkg_dict))

    # Use importlib.metadata to query installed package versions in Python
    importlib = pyimport("importlib.metadata")
    
    all_python_ok = true
    println("-"^57)
    @printf("%-20s | %-15s | %-10s\n", "Python Package", "Status", "Version")
    println("-"^57)

    for pkg_name in pkg_names
        status = "FAILED"
        version = "N/A"
        try
            # Attempt to import and retrieve the version string
            p = pyimport(pkg_name)
            version = string(importlib.version(pkg_name))
            status = "OK"
        catch
            all_python_ok = false
        end
        @printf("%-20s | %-15s | %-10s\n", pkg_name, status, version)
    end
    println("-"^57)

    if all_python_ok
        println("\nSUCCESS: SOSS is fully installed and functional.")
    else
        println("\nWARNING: Some Python dependencies failed.")
        println("Please verify that the system meets the installation requirements.")
    end

catch e
    println("\nCRITICAL ERROR checking Python:")
    println(e)
    println("\nPlease verify that the system meets the installation requirements.")
end

println("---------------------------------------------------------")
