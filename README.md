# Dynamic bathymetry and land-sea mask tool for the coupled model POEM

This repository contains the code of for the sea level change module of the POEM climate model.
The code is currently designed for use with the ocean model MOM6, ice sheet model PISM and solid earth model
VILMA. 

## Prerequisites

This tool is written is Python 3 and requires the following packages:

- Numpy
- NetCDF4
- Re (regular expressions)
- Sys
- Copy
- Time
- Math

In addition, it requires the tool regrid runoff created by Alistair Adcroft
in order to re-generate river runoff maps when the land-sea mask changes. This tool can be found here:
https://github.com/adcroft/regrid_runoff

## Current functionality
Currently, the plan is to allow the tool to be used as part of the PISM-MOM
coupling framework of Kreuzer et al. 21, as well as standalone (e.g. when
you want sea level/ land-sea masks to change without dynamically coupling in
an ice sheet or solid earth model). 