"""Tsunami source generation using Okada (1985) fault models.

Provides the :class:`Tsunami` class for computing initial water surface
displacement from earthquake fault parameters. Supports reading fault
definitions from CSV (NOAA SIFT format), GeoTIFF, UCSB, and Australian
PTHA Excel files.  Uses Clawpack's ``dtopotools`` for the Okada computation.
"""

from pathlib import Path
from typing import Optional, Union

import numpy as np
import rasterio
import xarray as xr
from scipy.ndimage import gaussian_filter

from clawpack.geoclaw import dtopotools


class Tsunami:
    """Tsunami initial condition generator based on Okada (1985).

    The workflow is:

    1. Create a :class:`Tsunami` instance.
    2. Set fault parameters via :meth:`set_subfault` or read from file
       via :meth:`read_fault_file`.
    3. Call :meth:`compute` to generate the sea-floor displacement.
    4. Access the result via :attr:`data` (an ``xarray.Dataset`` with
       variables ``x``, ``y``, ``dZ``).

    Attributes
    ----------
    fault : dtopotools.Fault or None
        Clawpack fault object containing the subfaults.
    data : xr.Dataset or None
        Computed displacement field with ``dZ(y, x)`` in metres.
    """

    def __init__(self) -> None:
        self.fault: Optional[dtopotools.Fault] = None
        self.data: Optional[xr.Dataset] = None

    # ------------------------------------------------------------------
    # Fault definition — single subfault (for GUI use)
    # ------------------------------------------------------------------

    def set_subfault(
        self,
        longitude: float,
        latitude: float,
        depth: float = 10.0,
        length: float = 100.0,
        width: float = 50.0,
        strike: float = 0.0,
        dip: float = 15.0,
        rake: float = 90.0,
        slip: float = 5.0,
    ) -> None:
        """Define a single rectangular subfault with Okada parameters.

        Parameters
        ----------
        longitude, latitude : float
            Epicentre location in decimal degrees.
        depth : float
            Depth to top of fault plane in km.
        length : float
            Along-strike length in km.
        width : float
            Down-dip width in km.
        strike : float
            Strike angle in degrees (0–360, clockwise from north).
        dip : float
            Dip angle in degrees (0–90).
        rake : float
            Rake angle in degrees (−180 to 180).
        slip : float
            Slip magnitude in metres.
        """
        input_units = {
            "length": "km",
            "width": "km",
            "depth": "km",
            "slip": "m",
            "mu": "Pa",
        }
        subfault = dtopotools.SubFault()
        subfault.longitude = longitude
        subfault.latitude = latitude
        subfault.depth = depth
        subfault.length = length
        subfault.width = width
        subfault.strike = strike
        subfault.dip = dip
        subfault.rake = rake
        subfault.slip = slip
        subfault.coordinate_specification = "top center"
        self.fault = dtopotools.Fault(subfaults=[subfault], input_units=input_units)

    # ------------------------------------------------------------------
    # Read fault from file
    # ------------------------------------------------------------------

    def read_fault_file(self, file_name: Union[str, Path]) -> None:
        """Read fault parameters from a file.

        Dispatches to the appropriate reader based on file extension:
        ``.csv`` → :meth:`read_csvfault`,
        ``.tif`` / ``.tiff`` / ``.geo`` → :meth:`read_geotiff`,
        ``.ucsb`` → :meth:`read_ucsb`.

        Parameters
        ----------
        file_name : str or Path
            Path to the fault definition file.
        """
        ext = str(file_name).rsplit(".", 1)[-1].lower()
        if ext == "csv":
            self.read_csvfault(file_name)
        elif ext in ("geo", "tif", "tiff"):
            self.read_geotiff(file_name)
        elif ext in ("ucsb", "usb"):
            self.read_ucsb(file_name)
        else:
            raise ValueError(f"Unsupported fault file extension: .{ext}")

    def read_geotiff(self, geo_file: Union[str, Path]) -> None:
        """Read a pre-computed displacement field from a GeoTIFF.

        Parameters
        ----------
        geo_file : str or Path
            Path to the GeoTIFF file.
        """
        with rasterio.open(geo_file) as src:
            data = src.read(1)
            transform = src.transform
            x = transform * (np.arange(src.width) + 0.5, np.zeros(src.width))
            y = transform * (np.zeros(src.height), np.arange(src.height) + 0.5)
            ds = xr.Dataset()
            ds["x"] = xr.DataArray(x[0], dims=["x"])
            ds["y"] = xr.DataArray(y[1], dims=["y"])
            ds["dZ"] = xr.DataArray(data, dims=["y", "x"])
            ds.attrs["crs"] = src.crs.to_string()
        self.data = ds

    def read_csvfault(self, csv_file: Union[str, Path]) -> None:
        """Read fault parameters from a NOAA SIFT CSV file.

        Parameters
        ----------
        csv_file : str or Path
            Path to the CSV file with columns: longitude, latitude,
            depth, length, width, strike, dip, rake, slip.
        """
        self.fault = dtopotools.Fault()
        column_map = {
            "longitude": 0,
            "latitude": 1,
            "depth": 2,
            "length": 3,
            "width": 4,
            "strike": 5,
            "dip": 6,
            "rake": 7,
            "slip": 8,
        }
        input_units = {
            "length": "km",
            "width": "km",
            "depth": "km",
            "slip": "m",
            "mu": "Pa",
        }
        self.fault.read(
            csv_file,
            column_map,
            skiprows=1,
            delimiter=",",
            input_units=input_units,
            coordinate_specification="noaa sift",
        )

    def read_ucsb(self, ucsb_file: Union[str, Path]) -> None:
        """Read fault parameters from a UCSB file.

        Parameters
        ----------
        ucsb_file : str or Path
            Path to the UCSB fault file.
        """
        self.fault = dtopotools.UCSBFault()
        self.fault.read(ucsb_file, rupture_type="static")

    def read_ptha(
        self,
        event_excel_file: Union[str, Path],
        event_row_number: int,
        statistics_excel_file: Union[str, Path],
    ) -> None:
        """Read fault parameters from Australian PTHA Excel files.

        Parameters
        ----------
        event_excel_file : str or Path
            Path to the PTHA event Excel file.
        event_row_number : int
            1-based row number of the event.
        statistics_excel_file : str or Path
            Path to the PTHA statistics Excel file.
        """
        import pandas as pd

        event_df = pd.read_excel(event_excel_file)
        stats_df = pd.read_excel(statistics_excel_file)

        event_row = event_df.iloc[event_row_number - 1]
        event_indices = list(
            map(int, event_row["event_index_string"].strip().split("-")[:-1])
        )
        slips = list(map(float, event_row["event_slip_string"].strip().split("_")[:-1]))

        input_units = {
            "length": "km",
            "width": "km",
            "depth": "km",
            "slip": "m",
            "mu": "Pa",
        }

        fault_segments = []
        for index, slip in zip(event_indices, slips):
            stats_row = stats_df.iloc[index - 1]
            subfault = dtopotools.SubFault()
            subfault.longitude = stats_row["lon_c"]
            subfault.latitude = stats_row["lat_c"]
            subfault.depth = stats_row["depth"]
            subfault.strike = stats_row["strike"]
            subfault.dip = stats_row["dip"]
            subfault.rake = stats_row["rake"]
            subfault.slip = slip
            subfault.length = stats_row["length"]
            subfault.width = stats_row["width"]
            subfault.coordinate_specification = "centroid"
            fault_segments.append(subfault)

        self.fault = dtopotools.Fault(subfaults=fault_segments, input_units=input_units)

    # ------------------------------------------------------------------
    # Compute displacement
    # ------------------------------------------------------------------

    def compute(
        self,
        dx: float = 1.0 / 60.0,
        smoothing: bool = True,
        buffer_size: float = 5.0,
        sigma: int = 5,
    ) -> None:
        """Compute sea-floor displacement using the Okada (1985) model.

        The result is stored in :attr:`data` as an ``xarray.Dataset``
        with variable ``dZ(y, x)`` in metres, georeferenced to EPSG:4326.

        Parameters
        ----------
        dx : float
            Grid spacing in degrees.
        smoothing : bool
            Apply Gaussian smoothing to the displacement field.
        buffer_size : float
            Buffer around the fault in degrees.
        sigma : int
            Standard deviation for the Gaussian filter (in grid cells).
        """
        if self.fault is None:
            raise RuntimeError(
                "No fault defined. Call set_subfault() or read_fault_file() first."
            )

        x, y = self.fault.create_dtopo_xy(buffer_size=buffer_size, dx=dx)
        dtopo = self.fault.create_dtopography(x, y)

        dZ = np.squeeze(dtopo.dZ[-1, :, :])

        if smoothing:
            dZ = gaussian_filter(dZ, sigma=sigma)

        ds = xr.Dataset()
        ds["x"] = xr.DataArray(dtopo.x, dims=["x"])
        ds["y"] = xr.DataArray(dtopo.y, dims=["y"])
        ds["dZ"] = xr.DataArray(dZ, dims=["y", "x"])
        ds.rio.write_crs(4326, inplace=True)

        self.data = ds

    # ------------------------------------------------------------------
    # Output
    # ------------------------------------------------------------------

    def write(self, file_name: Union[str, Path], format: str = "netcdf") -> None:
        """Write the displacement field to a file.

        Parameters
        ----------
        file_name : str or Path
            Output file path.
        format : str
            Output format: ``"netcdf"`` or ``"geotiff"``.
        """
        if self.data is None:
            raise RuntimeError("No data to write. Call compute() first.")
        if format == "netcdf":
            self.data.to_netcdf(file_name)
        else:
            raise ValueError(f"Format not supported: {format}")
