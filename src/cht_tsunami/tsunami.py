import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import pandas as pd
import rasterio
from scipy.ndimage import gaussian_filter

from clawpack.geoclaw import dtopotools

# Lauren - user defined input, geotiff, siftcsv, or aus_ptha
def get_source_file_details():
    source_type = input("Enter the source file type (geotiff, siftcsv, or aus_ptha): ").strip().upper()

    if source_type == "GEOTIFF":
        file_name = input("Enter the name of the geotiff file: ").strip()
        return {"source_type":source_type, "file_name": file_name}
    
    elif source_type == "SIFTCSV":
        file_name = input("Enter the name of the csv file: ").strip()
        return {"source_type": source_type, "file_name": file_name}   
    
    elif source_type == "AUS_PTHA":
        event_excel_file = input("Enter the name of the Australian PTHA Event excel file: ").strip()
        event_row_number = int(input("Enter the event row number: ").strip())
        statistics_excel_file = input("Enter the name of the Australian PTHA Statistics excel file: ").strip()
        return {
            "source_type": source_type,
            "event_excel_file": event_excel_file,
            "event_row_number": event_row_number,
            "statistics_excel_file": statistics_excel_file
        } 
    
    else:
        print("Invalid source file type. Please enter either geotiff, sift, csv or aus_ptha.")
        return get_source_file_details()

class Tsunami:
    def __init__(self, geo_file=None, csv_file=None, event_excel_file=None, ucsb_file=None,
                 event_row_number=None, statistics_excel_file=None,
                 plot=False,
                 smoothing=False,
                 compute=True,
                 dx=1.0/60.0,
                 buffer_size=5.0):
        if geo_file is not None:
            self.read_geotiff(geo_file)
            compute = False
        elif csv_file is not None:
            self.read_csvfault(csv_file, plot=plot)
        elif ucsb_file is not None:
            self.read_ucsb(ucsb_file, plot=plot)
        elif event_excel_file is not None and event_row_number is not None and statistics_excel_file is not None:
            self.read_ptha(event_excel_file, event_row_number, statistics_excel_file, plot=plot)

        if compute:
            self.compute(smoothing=smoothing, dx=dx, buffer_size=buffer_size)

    def read_fault_file(self, file_name):
        """Read the fault file and store the data in the class"""
        # Get extension of the file
        extension = file_name.split(".")[-1].lower()
        if extension == "csv":
            self.read_csvfault(file_name)
        elif extension == "geo" or extension == "tif" or extension == "tiff":
            self.read_geotiff(file_name)
        else:
            # 
            print("File type not supported")

    def read_geotiff(self, geo_file):
        with rasterio.open(geo_file) as src:
            data = src.read(1)  # Reading the first band
            transform = src.transform
            x = transform * (np.arange(src.width) + 0.5, np.zeros(src.width))
            y = transform * (np.zeros(src.height), np.arange(src.height) + 0.5)
            ds = xr.Dataset()
            ds['x'] = xr.DataArray(x[0], dims=['x'])
            ds['y'] = xr.DataArray(y[1], dims=['y'])
            ds['dZ'] = xr.DataArray(data, dims=['y', 'x'])
            ds.attrs['crs'] = src.crs.to_string()
        self.data = ds

    def read_csvfault(self, csv_file, plot=False):
        """Use clawpack scripts to read in the input file and store the data in the class"""
        self.fault = dtopotools.Fault()
        column_map = {"longitude":0, "latitude":1, "depth":2, "length":3, "width":4, "strike":5,  "dip":6, "rake":7, "slip":8}
        input_units = {'length': 'km', 'width': 'km', 'depth': 'km', 'slip': 'm', 'mu': 'Pa'}
        self.fault.read(csv_file, column_map, skiprows=1, delimiter=",", input_units=input_units,coordinate_specification="noaa sift")
        if plot:
            self.fault.plot_subfaults()

    def read_ucsb(self, ucsb_file, plot=False):
        """Use clawpack scripts to read in the input file and store the data in the class"""
        self.fault = dtopotools.UCSBFault()
        column_map = {"longitude":0, "latitude":1, "depth":2, "length":3, "width":4, "strike":5,  "dip":6, "rake":7, "slip":8}
        input_units = {'length': 'km', 'width': 'km', 'depth': 'km', 'slip': 'm', 'mu': 'Pa'}
        self.fault.read(ucsb_file, rupture_type="static")
        if plot:
            self.fault.plot_subfaults()

    def read_ptha(self, event_excel_file, event_row_number, statistics_excel_file, plot=False):
        """Read the event and statistics Excel files for PTHA and store the data in the class"""
        event_df = pd.read_excel(event_excel_file)
        stats_df = pd.read_excel(statistics_excel_file)

        event_row = event_df.iloc[event_row_number - 1] # Adjust for zero-based index
        event_indices = list(map(int, event_row['event_index_string'].strip().split('-')[:-1]))  # Remove the trailing empty string
        slips = list(map(float, event_row['event_slip_string'].strip().split('_')[:-1]))  # Remove the trailing empty string

        input_units = {'length': 'km', 'width': 'km', 'depth': 'km', 'slip': 'm', 'mu': 'Pa'}

        fault_segments = []
        for index, slip in zip(event_indices, slips):
            stats_row = stats_df.iloc[index - 1]  # Adjust for zero-based index
            subfault = dtopotools.SubFault()
            subfault.longitude=stats_row['lon_c']
            subfault.latitude=stats_row['lat_c']
            subfault.depth=stats_row['depth']
            subfault.strike=stats_row['strike']
            subfault.dip=stats_row['dip']
            subfault.rake=stats_row['rake']
            subfault.slip=slip
            subfault.length=stats_row['length']
            subfault.width=stats_row['width']
            subfault.coordinate_specification="centroid" # THIS IS A GUESS
            fault_segments.append(subfault)

        self.fault = dtopotools.Fault(subfaults=fault_segments,input_units=input_units)
        if plot:
            self.fault.plot_subfaults()


    def compute(self, dx=1.0/60.0, smoothing=True, buffer_size=5.0):
        """Compute tsunami with Okada (1985). Store initial water surface displacement in XArray"""
        x, y = self.fault.create_dtopo_xy(buffer_size=buffer_size, dx=dx)
        dtopo = self.fault.create_dtopography(x, y)

        # Take the last displacement
        dZ = np.squeeze(dtopo.dZ[-1,:,:])

        if smoothing:
            dZ = gaussian_filter(dZ, sigma=5)

        # x, y = np.meshgrid(dtopo.x, dtopo.y)
        x = dtopo.x
        y = dtopo.y

        # Create an xarray dataset with the surface displacement dtopo.dZ. X and Y are the coordinates dtopo.x and dtopo.y         
        ds = xr.Dataset()
        ds['x'] = xr.DataArray(x, dims=['x'])
        ds['y'] = xr.DataArray(y, dims=['y'])
        ds['dZ'] = xr.DataArray(dZ, dims=['y', 'x'])

        ds.rio.write_crs(4326, inplace=True)

        self.data = ds

    def interpolate(self, mesh):
        # mesh is an xugrid grid object, return xugrid data array

        # Interpolate the data to the mesh


        pass
    
    def plot_fullfault(self):
        fig, ax = plt.subplots()
        x = self.data["x"].values 
        y = self.data["y"].values
        z = self.data["dZ"].values
        max_abs_dz = max(abs(z.min()), abs(z.max()))
        c = ax.pcolor(x, y, z, cmap="seismic", vmin=-max_abs_dz, vmax=+max_abs_dz)
        ax.set_title('dZ')
        ax.axis('equal')
        fig.colorbar(c, ax=ax)
        fig.tight_layout()
        plt.show()


    def write(self, file_name, format="netcdf"):
        if format == "netcdf":
            self.data.to_netcdf(file_name)
        elif format == "geotiff":
            # Write geotiff file named file_name using PIL
            pass
        elif format == "xugrid":
            pass
        else:
            raise ValueError("Format not supported")
        
        