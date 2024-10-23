---
author: Daniele Guffanti - University & INFN Milano-Bicocca (daniele.guffanti@mib.infn.it)
date: 2024-10-22 08:25
---
# VisibilityMap Tools

This folder contains a few example macros for creating 3D/2D visibility maps
using the `PhotonVisibilityExport_module`. 

Running the macro `ttree_to_th3.C` will produce a file containing the detector visibility 
map stored in a three-dimensional THnF and the individual visibility maps for each OpDet
as `THnSparseF` objects.

The file `DUNEVisUtils.hpp` containes a few useful functions to handle these histograms, 
in particular to "rebin" the maps for a better display and to profile them along a given 
dimension. 

Example fhicl configurations can be found in the `fcl` folder. Note that in case of 
Xe doping, the visivility must be exported for Ar and Xe separately and combined
in a later stage with the appropriate weights. 
