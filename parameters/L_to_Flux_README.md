# Luminosity to Flux tabular table description

The Luminosity to Flux tabular tables are located at ``./parameters/Models/L_to_Flux/`` and contains the mapping
 between Luminosity to Flux. The tabular files are in fits format, and 4 example files are provided.
 The Luminosity band should be consistent with the scaling relation used to assign the luminosities, and the count rate
 should be in the energy band of interest, and the unit for the flux should be in [ergs/s/cm^2] unit.

The code expects 3D array, where the first dimension defines the plasma temperature,
 the second column defines the redshift of object, and finally the third column defines the hydrogen column density.

Caveat
------

The current version does not take the hydrogen column density into account so it takes the 1st entry and downgrade the 
matrix in to 2D matrix. This is there for the future release of the code. 
