# Import Julia package manager module
using Pkg

# Install Julia packages from Project.toml
Pkg.instantiate()

# Get Python path from python_path.txt
python_path = strip(read("python_path.txt", String))

# Setup Python path
ENV["PYTHON"] = python_path
Pkg.build("PyCall")

# Exit Julia
exit()