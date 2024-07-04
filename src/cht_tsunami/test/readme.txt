
# git clone: 
# https://github.com/Deltares-research/cht_tsunami
# https://github.com/Deltares-research/CoastalHazardsToolkit

#Open miniconda powershell (or anaconda powershell)

# Create new environment
conda create -n tsunami python=3.10

# Activate environment
conda activate tsunami

# Make editable package of cht_tsunami
pip install -e c:\work\checkouts\git\cht_tsunami

# Make editable package of cht (used for SFINCS fileio etc)
pip install -e c:\work\checkouts\git\coastalhazardstoolkit

# dependencies needed
conda install dask 
conda install openpyxl