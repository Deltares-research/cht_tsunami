import matplotlib.pyplot as plt

from cht_tsunami.tsunami import Tsunami
from cht.sfincs2.sfincs import SFINCS

#file_name = "c:\\work\\checkouts\\git\\cht_tsunami\\src\\clawpack\\geoclaw\\data\\info_sz.dat.txt"
source_file = "alaska1964.csv"

# Generate tsunami
ts = Tsunami(source_file=source_file, plot=True)

x = ts.data["x"].values
y = ts.data["y"].values
z = ts.data["dZ"].values

fig, ax = plt.subplots()
c = ax.pcolor(x, y, z, cmap="seismic", vmin=-20.0, vmax=20.0)
ax.set_title('dZ')
ax.axis('equal')
fig.colorbar(c, ax=ax)
fig.tight_layout()

plt.show()

# # Read SFINCS model
# sf = SFINCS(file_name="sfincs.inp")

# # Interpolate the data to the mesh
# sf.initial_conditions.interpolate(ts.data)

# # Write the initial conditions to the SFINCS input file
# sf.input.inifile = "sfincs_ini.nc"
# sf.initial_conditions.write()
