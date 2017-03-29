# ``Map_Parameters.json`` parameter description

- ``Unit`` is the pixel unit. This is a dummy variable and he code does not
 look it up. It is just used for the reference and plotting purpose. 
 
- ``Flux_energy_band`` is the pixel unit. This is a dummy variable and he code does not
 look it up. It is just used for the reference and
  plotting purpose and labeling convention.

- ``DEC_map_size`` is the total height of the map. It is in degree unit. 
The code make a surface brightness map of with in ``[-[DEC_map_size]/2, +[[DEC_map_size]/2]`` range.

- ``RA_map_size`` is the total width of the map. It is in degree unit. 
The code make a surface brightness map of with in ``[-[RA_map_size]/2, +[[RA_map_size]/2]`` range.

- ``Num_pixels`` is the number of pixel on each side of the map. 
So the pixel resolution would be
 ``[DEC_map_size]/[Num_pixels] x [RA_map_size]/[Num_pixels]`` square degree

