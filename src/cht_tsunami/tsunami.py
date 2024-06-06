import xarray as xr

from clawpack.geoclaw import dtopotools

class Tsunami:
    def __init__(self, file_name=None):
        if file_name is not None:
            # Read input file
            self.read_input(file_name)
            # Compute the initial water surface displacement
            self.compute(file_name)

    def read_input(self, file_name):
        """Use clawpack scripts to read in the input file and store the data in the class"""
        pass

    def compute(self, dx=100.0, dy=100):
        """Compute tsunami with Okada (1985). Store initial water surface displacement in XArray"""

        ds = xr.Dataset()
        ds['z'] = xr.DataArray([0.0], dims=['x', 'y'])

        self.dataset = ds
        pass

    def interpolate(self, mesh):
        # mesh is an xugrid grid object, return xugrid data array

        # Interpolate the data to the mesh


        pass

