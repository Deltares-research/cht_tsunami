Getting started
===============

Installation
------------

Install from GitHub:

.. code-block:: bash

   pip install git+https://github.com/Deltares-research/cht_tsunami.git

This also installs `Clawpack <https://www.clawpack.org/>`_ as a dependency,
which is used internally for the Okada displacement computation.

For development:

.. code-block:: bash

   git clone https://github.com/Deltares-research/cht_tsunami.git
   cd cht_tsunami
   pip install -e .

Quick example -- compute displacement from a fault
----------------------------------------------------

Define a simple fault and compute the seafloor displacement on a regular grid:

.. code-block:: python

   import numpy as np
   from cht_tsunami.tsunami import Tsunami

   # Create a Tsunami object with a regular grid
   ts = Tsunami()

   # Define grid (longitude, latitude)
   lon = np.linspace(140.0, 145.0, 500)
   lat = np.linspace(35.0, 42.0, 700)

   # Set a single subfault (Okada parameters)
   ts.set_subfault(
       longitude=142.37,    # fault centroid longitude
       latitude=38.32,      # fault centroid latitude
       depth=10.0,          # depth to top of fault (km)
       length=200.0,        # fault length (km)
       width=100.0,         # fault width (km)
       strike=198.0,        # strike angle (degrees)
       dip=10.0,            # dip angle (degrees)
       rake=90.0,           # rake angle (degrees, 90 = pure thrust)
       slip=15.0,           # slip magnitude (m)
   )

   # Compute displacement field
   dz = ts.compute_displacement(lon, lat)

   # dz is a 2-D array (len(lat) x len(lon)) of vertical displacement in metres
   print(f"Max uplift:     {dz.max():.2f} m")
   print(f"Max subsidence: {dz.min():.2f} m")

Quick example -- use GEM fault database
-----------------------------------------

Look up a fault from the GEM Global Active Faults database and compute
its Okada parameters:

.. code-block:: python

   from cht_tsunami.faults import get_faults, get_okada_params_from_fault

   # Download the GEM fault database (cached after first call)
   faults_gdf = get_faults()

   # Find a specific fault by name
   fault = faults_gdf[faults_gdf["fault_name"] == "Nankai"].iloc[0]

   # Convert to Okada parameters
   params = get_okada_params_from_fault(fault, slip=5.0)

   print(params)
   # {'longitude': ..., 'latitude': ..., 'depth': ..., 'length': ...,
   #  'width': ..., 'strike': ..., 'dip': ..., 'rake': ..., 'slip': 5.0}
