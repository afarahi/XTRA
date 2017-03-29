# X-TRA: A Template-based Approach for Modeling X-ray Survey Yields of Groups and Clusters

X-TRA - X-ray Template Realization Algorithm - a template approach to realizing large-scale,
synthetic X-ray emission maps from groups and clusters of galaxies. The method applies scaling
relations and template emission shapes to an input set of halo masses and redshifts 
(derived from N-body simulation lightcone outputs), producing raw surface brightness maps in the
X-ray wavelength. The approach offers a fast, flexible means to characterize transfer functions
and systematic error sources in X-ray surveys.


Installation
------------ 

First, clone the repo with 

    $ git clone "https://github.com/afarahi/XXX.git"

The code is developed in Python 2.7, so should work out of the box, 
if you have standard astrophysical libraries installed already.

The dependencies are ``numpy``, ``matplotlib``, and ``astropy`` libraries. 

Running
-------

Running the code is simple. 

1. First, the user need to specify the input/output directories.

    ``Directories.json`` specifies where the code need to look for input catalog and where it
    saves the outputs. ``Halo-dir`` in the input directory for the halo catalog. 
    ``Cluster_dir`` is where it saves the cluster catalog after assigning the Luminosity, Temperature, and flux.
    ``SB_Map_dir`` is where it saves the surface brightness maps, and finally 
    ``Event_Map_dir`` is where it saves the mocked event files.

2. ``python main.py 1 [input-filename]`` : make the cluster catalog and save the output at ``Cluster_dir``

3. ``python main.py 2 [input-filename]`` : read the cluster catalog from
  ``Cluster_dir`` and makes the surface brightness map and save it as a fit file at ``SB_Map_dir``
  
4. ``python main.py 3 [input-filename]`` : read the surface brightness map and produces
 tiled synthetic images of noised XMM exposures 


Specifying Model Parameters
---------------------------

1. ``./parameters/Cosmological_Parameters.json`` : specifies the cosmological parameters. It will be used as reference 
 and to calculate the evolution factor and the angular distance.
  
2. ``./parameters/Input_Parameters.json`` : specifies the free parameters of the model, e.g. scaling relation, 
the surface brightness profile shape. For more detail look at the documentation provided
 at [./parameters/Input_README.md](./parameters/Input_README.md).

3. ``./parameters/Map_Parameters.json`` : specifies the scales of surface brightness maps. 
For more detail look at the documentation provided 
at [./parameters/Map_README.md](./parameters/Map_README.md).

4. ``./parameters/Event_Map_Parameters.json`` : specifies the tiling scheme and the exposure time. 
For more detail look at the documentation provided at 
[./parameters/Event_Map_README.md](./parameters/Event_Map_README.md).

Caveats
-------

- The current version of the code expects the surface brightness maps be in unit of [ergs/s/cm^2/pixle] in 
soft band [0.5-2.0] keV in order to make the event maps. If the user does not want to generate the event maps, then 
any unit or X-ray band would works fine.

- The current version of the code assumes a 256 x 256 pixel paintings with 4.3" resolution. The event map properties
are fixed. If you are interested in different setting the best way would be reaching out to Arya Farahi, 
<a href="mailto:aryaf@umich.edu">aryaf@umich.edu</a>.

Using code
-----------------

We are glad, you're interested in using code to X-TRA. If there is anything we can help you with feel free to reach us
at <a href="mailto:aryaf@umich.edu">aryaf@umich.edu</a> for any question or request. 

If you want to modify the code for your own purpose, we would much rather you email us asking to collaborate than modify
 the code as a black box.

**License**

This project is [licensed](./LICENSE.md) under the terms of the MIT license.

**Citation**

 Please cite XXX, if you make use of these codes in your work.
 It is appreciated if you email <a href="mailto:aryaf@umich.edu">aryaf@umich.edu</a> to tell us that you are using it. 