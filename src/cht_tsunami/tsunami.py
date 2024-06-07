import xarray as xr
import numpy as np

from clawpack.geoclaw import dtopotools

class Tsunami:
    def __init__(self, source_file=None, plot=False):
        if source_file is not None:
            # Read input file (assuming it is a SIFT fault for now)
            self.read_fault(source_file, plot=plot)
            # Compute the initial water surface displacement
            self.compute()

    def read_geotiff(self, geotiff_file):
        # read geotiff

        # Create an xarray dataset with the surface displacement dtopo.dZ. X and Y are the coordinates dtopo.x and dtopo.y         
        ds = xr.Dataset()
        ds['dZ'] = xr.DataArray(dZ, dims=['x', 'y'])
        ds['x'] = xr.DataArray(x, dims=['x', 'y'])
        ds['y'] = xr.DataArray(y, dims=['x', 'y'])

        self.data = ds


    def read_fault(self, source_file, plot=False):
        """Use clawpack scripts to read in the input file and store the data in the class"""
        self.fault = dtopotools.Fault()
        # Assume for now that this is a CSV fault
        column_map = {"longitude":0, "latitude":1, "depth":2, "length":3, "width":4, "strike":5,  "dip":6, "rake":7, "slip":8}
        input_units = {'length': 'km', 'width': 'km', 'depth': 'km', 'slip': 'm', 'mu': 'Pa'}
        self.fault.read(source_file, column_map, skiprows=1, delimiter=",", input_units=input_units)
        if plot:
            self.fault.plot_subfaults()

#        # Assume for now that this is a SIFT fault
#        self.fault = dtopotools.SiftFault(source_file=source_file)

    def read_oz_ptha_files(self):
        pass


    def compute(self, dx=100.0, dy=100):
        """Compute tsunami with Okada (1985). Store initial water surface displacement in XArray"""
        x, y = self.fault.create_dtopo_xy()
        dtopo = self.fault.create_dtopography(x, y)

        # Take the last displacement
        dZ = np.squeeze(dtopo.dZ[-1,:,:])
        x, y = np.meshgrid(dtopo.x, dtopo.y)

        # Create an xarray dataset with the surface displacement dtopo.dZ. X and Y are the coordinates dtopo.x and dtopo.y         
        ds = xr.Dataset()
        ds['dZ'] = xr.DataArray(dZ, dims=['x', 'y'])
        ds['x'] = xr.DataArray(x, dims=['x', 'y'])
        ds['y'] = xr.DataArray(y, dims=['x', 'y'])
        # plt = dtopo.plot_dZ_colors(0.0)
        # plt.draw()

        self.data = ds

    def interpolate(self, mesh):
        # mesh is an xugrid grid object, return xugrid data array

        # Interpolate the data to the mesh


        pass


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
        
        