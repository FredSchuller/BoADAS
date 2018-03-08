# saboca-secondary-fluxes.py
# Fluxes for secondary calibrators, version 2010-May-05
# Created 2010-May-05 - G. Siringo
#
# This table is used by the BoA script scalib.boa - called by "redscal" from BoA
# Note: use ONLY UPPERCASE in the source names, the reduction script relies on that

# NOTE: as of May 2010 no reliable determination of 
# calibrators yet available - values for testing
# to be replaced soon
# Jun 2010: 1st iteration values entered - MDU

# May 2013: update list and flux values - FSc

calibFluxes = {
    'G34.3'	: 386.9,
    'G10.47'	: 335.5,
    'G10.62'	: 308.0,
    'G5.89'	: 216.5,
    'B13134'	: 155.0,
    'IRAS16293' : 113.7,
    'K3-50A'	: 104.4,
    'G45.1'	:  92.8,
    'NGC4945'	:  68.5,
    'CRL2688'	:  45.2,
    'NGC253'	:  39.6,
    'CRL618'	:  22.5,
    'HL_TAU'	:  15.8,
    'V883-ORI'	:  12.0,
    'ARP220'	:  10.2,
    'CW-LEO'    :  22.4,
    'VYCMA'     :  11.1
    }


# Add a few synonyms
calibFluxes['VY-CMA']    = calibFluxes['VYCMA']
calibFluxes['J1800-241'] = calibFluxes['G5.89']
calibFluxes['IRC+10216'] = calibFluxes['CW-LEO']
calibFluxes['N253'] = calibFluxes['NGC253']
calibFluxes['NGC-253'] = calibFluxes['NGC253']
#calibFluxes['IRAS16342-38'] = calibFluxes['IRAS16342']
calibFluxes['HLTAU'] = calibFluxes['HL_TAU']
