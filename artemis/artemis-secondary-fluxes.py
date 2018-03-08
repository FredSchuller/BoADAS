# artemis-secondary-fluxes.py
# Fluxes for ArTeMiS secondary calibrators
#
# This table is used by the BoA function redcal defined in artemis.boa

# June 2017: created (copy of saboca-secondary-fluxes.py)
# 2017-06-09 - fschuller  : updated values for 350, added values for 450
# 2018-01-19 - fschuller  : updated all values, based on all 2017 data

calibF350 = {
    'G10.47'	: 397.26, # +/- 31.93
    'B13134'	: 163.67, # +/- 22.80
    'G34.3'	: 424.40, # +/- 38.87
    'G45.1'	:  96.31, # +/-  9.33
    'G10.62'	: 365.69, # +/- 41.60
    'G5.89'	: 247.90, # +/- 26.68 
    'I16293'    : 149.73, # +/- 11.61
    'NGC253'	:  41.59, # +/-  5.03
    'ARP220'	:   8.45, # +/-  1.08
    'NGC4945'	:  65.31, # +/- 11.06
    'Carina'    :  87.96, # +/- 12.85
    'VYCMA'     :  11.76, # +/-  1.65
    'CRL618'    :  17.89, # +/-  2.78
    'K3-50A'	:  96.30, # +/-  4.20
    'N2071IR'   :  55.25, # +/-  3.75
    'CRL2688'   :  34.05, # +/-  2.75
    'CW-LEO'    :  23.13  # +/-  3.32
    }

calibF450 = {
    'G10.47'	: 198.99, # +/- 14.60
    'B13134'	:  83.19, # +/-  5.89
    'G34.3'	: 236.14, # +/- 18.88
    'G45.1'	:  44.83, # +/-  2.97
    'G10.62'	: 186.46, # +/- 18.10
    'G5.89'	: 121.89, # +/- 10.26
    'I16293'    :  83.92, # +/-  7.78
    'NGC253'	:  19.34, # +/-  2.29
    'ARP220'	:   4.13, # +/-  0.76
    'NGC4945'	:  33.85, # +/-  3.77
    'Carina'    :  71.87, # +/-  8.99
    'VYCMA'     :   6.81, # +/-  0.85
    'CRL618'    :  12.91, # +/-  1.36
    'K3-50A'    :  59.40, # +/-  1.00
    'N2071IR'   :  36.05, # +/-  1.85
    'CRL2688'   :  30.85, # +/-  7.25
    'CW-LEO'    :  12.99  # +/-  2.24
    }

# Add a few synonyms
calibF350['IRC+10216'] = calibF350['CW-LEO']
calibF350['N253']      = calibF350['NGC253']
calibF350['NGC-253']   = calibF350['NGC253']
calibF350['G10.47B1']  = calibF350['G10.47']
calibF350['IRAS16293'] = calibF350['I16293']
calibF350['VY-CMA']    = calibF350['VYCMA']

calibF450['IRC+10216'] = calibF450['CW-LEO']
calibF450['N253']      = calibF450['NGC253']
calibF450['NGC-253']   = calibF450['NGC253']
calibF450['G10.47B1']  = calibF450['G10.47']
calibF450['IRAS16293'] = calibF450['I16293']
calibF450['VY-CMA']    = calibF450['VYCMA']


