# ``Input_Parameters.json`` parameter description




- ``Energy_Unit`` is the energy unit, the only viable option is ergs in this version.

- ``Mass_limit`` is the mass limit, halos above this mass get assigned X-ray observables. 
The mass is in unit of Msun/h.

- ``F_limit`` flux limit, cluster below this flux limit does not contribute to
 the surface brightness map. The unit is in ergs/s/cm^2.
 
- ``Lxm_model`` specifies the Lx - Mass scaling model. There should be a 
``[Lxm_model]_parameters.json`` file under ``./parameters/Lxm/`` directory which 
specifies the scaling relatoins. The user may use the ``default`` setup.

- ``Txm_model`` specifies the Tx - Mass scaling model. There should be a 
``[Txm_model]_parameters.json`` file under ``./parameters/Txm/`` directory which 
specifies the scaling relatoins. The user may use the ``default`` setup.

- ``rTL_correlation`` specifies the correlation coefficient between Lx and Tx at 
fixed halo mass.

- ``L_to_Flux_model_file_name`` specifies the file that contains the mapping
 between Luminosity to the count rate. 
 The details of the tabular table are explained [...]. This should provided by
 the user and depends on the input luminosity and the count rate in the energy 
 band of interest.
 
- ``Flux_energy_band`` and string which specifies the energy band of the output flux.
This is a dummy variable and the code does not look it up. It is just used for the
 reference and plotting purpose. 

- ``Surface_Brightness`` specifies the surface brightness profile model. 
In the current version, only the beta profile is implemented. 

- ``SB_beta_bar`` is the expected beta, which is used in the beta profile.

- ``SB_beta_sig`` is the lognormal scatter for the value of beta. The beta value 
changes from a cluster to another.

- ``xc_bar`` is the expected value of the core radius over R500. 

- ``xc_sig`` is the lognormal scatter for the x_c. The x_c value 
changes from a cluster to another.

- ``xc_Beta_r`` is the correlation between beta and x_c.
  
 
