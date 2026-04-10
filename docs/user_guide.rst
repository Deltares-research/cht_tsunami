User guide
==========

Overview
--------

A tsunami simulation workflow with ``cht_tsunami`` typically involves:

1. **Define a fault** -- either manually (Okada parameters) or from the GEM
   fault database.
2. **Compute displacement** -- the Okada (1985) model translates fault slip
   into a seafloor vertical displacement field.
3. **Use the displacement** -- the resulting 2-D array can be used as an
   initial condition for a shallow-water propagation model (e.g. SFINCS,
   Delft3D, or GeoClaw).

GEM fault database
------------------

The `GEM Global Active Faults database
<https://github.com/GEMScienceTools/gem-global-active-faults>`_ contains
fault traces with geometry, dip, rake, and slip-rate attributes for thousands
of active faults worldwide.

.. code-block:: python

   from cht_tsunami.faults import get_faults

   # Download and cache the GEM database (GeoDataFrame)
   faults = get_faults()

   print(f"Total faults: {len(faults)}")
   print(faults.columns.tolist())

The result is a :class:`geopandas.GeoDataFrame` with LineString geometries
and attributes including ``fault_name``, ``dip_dir``, ``lower_seis_depth``,
``upper_seis_depth``, ``net_slip_rate``, and ``rake``.

Converting faults to Okada parameters
--------------------------------------

:func:`~cht_tsunami.faults.get_okada_params_from_fault` converts a GEM fault
record into the parameters needed by the Okada displacement model:

.. code-block:: python

   from cht_tsunami.faults import get_faults, get_okada_params_from_fault

   faults = get_faults()
   fault = faults[faults["fault_name"] == "Manila Trench"].iloc[0]

   params = get_okada_params_from_fault(fault, slip=8.0)

The returned dict contains:

- ``longitude``, ``latitude`` -- fault centroid
- ``depth`` -- depth to the top of the rupture (km)
- ``length``, ``width`` -- fault dimensions (km)
- ``strike`` -- strike angle (degrees, clockwise from north)
- ``dip`` -- dip angle (degrees, 0 = horizontal, 90 = vertical)
- ``rake`` -- slip direction (degrees, 90 = pure thrust, 0 = left-lateral)
- ``slip`` -- total slip (metres, as provided by the user)

Computing seafloor displacement
---------------------------------

The :class:`~cht_tsunami.tsunami.Tsunami` class wraps the Okada computation:

.. code-block:: python

   import numpy as np
   from cht_tsunami.tsunami import Tsunami

   ts = Tsunami()

   # Set fault parameters (can call multiple times for multiple subfaults)
   ts.set_subfault(
       longitude=120.5,
       latitude=15.0,
       depth=15.0,
       length=300.0,
       width=80.0,
       strike=0.0,
       dip=20.0,
       rake=90.0,
       slip=10.0,
   )

   # Define output grid
   lon = np.linspace(118.0, 123.0, 500)
   lat = np.linspace(12.0, 18.0, 600)

   # Compute vertical displacement (metres)
   dz = ts.compute_displacement(lon, lat)

The displacement is computed using the Okada (1985) analytical solution
as implemented in Clawpack/GeoClaw. The result is a 2-D NumPy array
suitable for use as a tsunami initial condition.

Multiple subfaults
------------------

For complex rupture scenarios, call ``set_subfault()`` multiple times before
computing displacement. Each subfault contributes linearly to the total
displacement field:

.. code-block:: python

   ts = Tsunami()

   # Northern segment
   ts.set_subfault(longitude=142.0, latitude=39.0, depth=10.0,
                   length=100.0, width=50.0, strike=200.0,
                   dip=12.0, rake=90.0, slip=12.0)

   # Southern segment
   ts.set_subfault(longitude=142.5, latitude=37.5, depth=15.0,
                   length=120.0, width=60.0, strike=195.0,
                   dip=15.0, rake=85.0, slip=8.0)

   dz = ts.compute_displacement(lon, lat)
