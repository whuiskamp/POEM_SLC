#!/usr/bin/env python3

# For standalone ocean-ice runs, use this script to calculate 
# the virtual mass flux to change the sea level by a given
# amount (m) over a given time (years). This assumes density of 
# fresh water is 1000kg/m3 and that model is using a noleap 
# calendar.
# 
# vprec = calc_mass_flux(slc,time)
# Input:
# - slc:   The required sea level change you wish to calculate,
#          e.g. m/100years
# - time:  Time over which the slc should occur (years)
# 
# Output:
# - vprec: The virtual mass flux applied to the model in
#          kg/m-2/s-1

import argparse

########### Argument parsing ########### 
parser = argparse.ArgumentParser(
            description=
            "Calculate the virtual mass flux (in kg/m-2/s-1) the model requires to \
                change global mean sea level by a given rate (m/yr)"
            )
parser.add_argument('-s', '--slc', type = float, action="store", dest="slc",
                        required=True, 
                        help="Total change in sea level")
parser.add_argument('-t', '--time', type = float, action="store", dest="time",
                        required=True, 
                        help="Time over which change in sea level occurs")
args = parser.parse_args()

########### Calculate the flux #########
def calc_vprec():
    t_year = 60*60*24*365 # No. of seconds in a year (noleap)

    # Calculate chosen sea level change as a mass flux (assuming fresh water flux)...
    mass = args.slc * 1000 # kg/m-2 
    print('Calculating: ')
    print(str(args.slc)+'m of sea level change over '+str(args.time)+' years')
     
    vprec = mass/(t_year*args.time) # kg/m-2/s-1
    print('Resulting in a flux of '+str(vprec)+' kg/m-2/s-1')
    return vprec

if __name__ == "__main__":

    calc_vprec()