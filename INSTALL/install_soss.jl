# ==========================================================
# SOSS Installation Script
# ==========================================================
using Pkg

println("---------------------------------------------------------")
println("⚙️  Starting SOSS Automatic Installation...")
println("---------------------------------------------------------")

try
    # 1. Activate the project environment in the current directory
    Pkg.activate(".")
    
    # 2. Downlad and install Julia and Python dependencies (from Project.toml and CondaPkg.toml)
    println("Downloading and configuring dependencies (Julia & Python)...")
    Pkg.instantiate()
    
    # 3. Precompile modules
    println("⚡ Precompiling modules...")
    Pkg.precompile()

    println("\nSUCCESS: SOSS has been installed correctly.")
catch e
    println("\nERROR during installation:")
    println(e)
    println("\nPlease verify that the system meets the installation requirements.")
end
