# secondary-fluxes.boa
# Fluxes for secondary calibrators, version 2008 Nov.
# Created 2008/11/21 - F. Schuller
#         2009/5/07: added Carina

# This dictionary is used by the script: reduce-calibrator.boa
# Note: use ONLY UPPERCASE in the source names, the reduction script
# relies on that

calibFluxes = {
    'CRL618'   : 4.8,
    'V883-ORI' : 1.36,
    'N2071IR'  : 9.1,
    'VYCMA'    : 1.53,
    'B13134'   : 12.9,
    'IRAS16293': 16.1,
    'G10.62'   : 33.0,
    'G34.3'    : 55.3,
    'G45.1'    : 8.0,
    'K3-50A'   : 14.7,
    'CRL2688'  : 5.5,
    'CW-LEO'   : 4.1,
    'HLTAU'    : 2.3,
    'G5.89'    : 27.6,
    'CARINA'   : 34.1
    }

# Note: still using Dec. 2007 values for: CW-LEO, HLTAU and G5.89

# Add a few synonyms
calibFluxes['VY-CMA']    = calibFluxes['VYCMA']
calibFluxes['J1800-241'] = calibFluxes['G5.89']
calibFluxes['IRC+10216'] = calibFluxes['CW-LEO']
