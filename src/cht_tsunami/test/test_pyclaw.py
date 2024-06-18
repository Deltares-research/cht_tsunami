import matplotlib.pyplot as plt

from cht_tsunami.tsunami import Tsunami
from cht.sfincs2.sfincs import SFINCS

#file_name = "c:\\work\\checkouts\\git\\cht_tsunami\\src\\clawpack\\geoclaw\\data\\info_sz.dat.txt"
source_file = "C:/Users/Lauren/Desktop/Tools/cht_tsunami/src/cht_tsunami/test/alaska1964.csv"

# Generate tsunami
ts = Tsunami(source_file=source_file, plot=False, smoothing=True)

ts.write("tsunami.nc")

# Move this next plotting bit somehow to Tsunami object? 
# E.g. with: ts.plot() ?

fig, ax = plt.subplots()
x = ts.data["x"].values 
y = ts.data["y"].values
z = ts.data["dZ"].values
c = ax.pcolor(x, y, z, cmap="seismic", vmin=-20.0, vmax=20.0)
ax.set_title('dZ')
ax.axis('equal')
fig.colorbar(c, ax=ax)
fig.tight_layout()
plt.show()

# Read SFINCS model
sf = SFINCS(root=r"C:/Users/Lauren/Desktop/Tools/SFINCS_Tsunami/models/global", mode="r")

# Interpolate the data to the mesh
sf.initial_conditions.interpolate(ts.data, var_name="dZ")

# Write the initial conditions to the SFINCS input file
sf.input.variables.ncinifile = "sfincs_ini.nc"
sf.initial_conditions.write()

# Write new sfincs input file
sf.input.write()
