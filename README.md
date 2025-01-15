# Gravitational waves signals for neutron stars binaries with a cosmological constant in their cores
 Pipeline to generate simulated gravitational waves signals for neutron star binaries, described by a novel prescription for the equation of state (EOS), where a vacuum energy phase transition is triggered in the stellar core, presented in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.111.024001.

The work uses the PyCBC software package (https://pycbc.org/). The masses and tidal deformabilities of the two components of the binary are drawned from "TOV" data (containgin masses and tidal deformabilities for each EOS configuration). These data were generated on a separate analysis. The necessary code for obtaining them (written in Julia) is available in a dedicated repository (https://github.com/GiuliaVentagli/NSCC).

# Content
Here is a list of the files:
- GW_dataGen.ipynb: This is the main notebook (Python) which generate the gravitational waves signals. It uses threaded pool to efficiently generate signals from a large number of EOSs.
- model_GW_dataGen.py: The Python file containing the required functionality needed for the notebook "GW_dataGen.ipynb".
- TOVs_Thread.zip: a zipped folder contanining the following files:
  - matrixrho10000.dat and matrixcs10000.dat, datasets of random values for mass density and speed of sound used to construct the phenomenological EOS at high-density, used to generate the TOV files, and used here to extract the EOS nuclear properties.
  - TOV data, the nomenclature is TOV_\<lowdensityEOS\>_\<Lambda\>_\<z\>.csv, where lowdensityEOS is the low density EOS prescription used to create the data, Lambda is the value of the cosmological constant in the core and z is a bookkeeping parameter, corresponding to the zth line in matrixrho10000.dat and matrixcs10000.dat. For example, TOV_sly_-194.0_36.csv corresponds to a configuration with the SLy EOS for low density regime, a value of the cosmological constant in the core of -194.0 and a high density EOS described by the 36th line in matrixrho10000.dat and matrixcs10000.dat.
