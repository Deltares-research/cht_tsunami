import matplotlib.pyplot as plt

from cht_tsunami.tsunami import Tsunami, get_source_file_details
from cht_sfincs.sfincs import SFINCS

source_details = get_source_file_details()

# Depending on the source type, initialize the Tsunami object accordingly
if source_details["source_type"] == "GEOTIFF":
    ts = Tsunami(geo_file=source_details["file_name"])
elif source_details["source_type"] == "SIFTCSV":
    ts = Tsunami(csv_file=source_details["file_name"], plot=False, smoothing=True)
elif source_details["source_type"] == "AUS_PTHA":
    ts = Tsunami(event_excel_file=source_details["event_excel_file"], 
                 event_row_number=source_details["event_row_number"], 
                 statistics_excel_file=source_details["statistics_excel_file"], 
                 plot=False, smoothing=False)

ts.write("tsunami.nc")

ts.plot_fullfault()

# Read SFINCS model
sf = SFINCS(root=r"C:/Users/Lauren/Desktop/Tools/SFINCS_Tsunami/models/global", mode="r")

# Interpolate the data to the mesh
sf.initial_conditions.interpolate(ts.data, var_name="dZ")

# Write the initial conditions to the SFINCS input file
sf.input.variables.ncinifile = "sfincs_ini.nc"
sf.initial_conditions.write()

# Write new sfincs input file
sf.input.write()
