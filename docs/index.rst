cht_tsunami
###########

``cht_tsunami`` computes tsunami initial conditions from earthquake fault
parameters using the Okada (1985) elastic dislocation model. It provides
tools for accessing the GEM Global Active Faults database and translating
fault geometry into seafloor displacement fields suitable for tsunami
propagation models.

Under the hood, the actual displacement computation uses
`Clawpack/GeoClaw <https://www.clawpack.org/geoclaw/>`_, but users interact
only with the :class:`~cht_tsunami.tsunami.Tsunami` class and the
:mod:`~cht_tsunami.faults` utilities.

.. toctree::
   :maxdepth: 2
   :caption: Contents

   getting_started
   user_guide
   api
   changelog
