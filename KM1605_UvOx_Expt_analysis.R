# KM1605_UvOx_Expt_analysis.R

# Analysis of an ultraviolet radiation oxiation experiment conducted with large-volume whole seawater samples during a SCOPE cruise (KM1605) aboard the R/V Kilo Moana

# Created 11/29/16 by J.R.C.

# Follows same general initial path as PAL1314_PAL1516_environmental_samples.R

### workspace preparation ####

# load necessary libraries

library(LOBSTAHS)
library(chemCal)
library(stats)
library(stringr)

# set wd

setwd("/Users/jrcollins/Code/LipidPhotoOxBox/")

### load some necessary metadata; define functions ###

# some chemical data

DNPPE_mg_mL_BD_extracts_2016 = 0.051 # concentration of the DNPPE added during B&D extractions
# in Van Mooy Lab in 2016
DNPPE_MW = 875.081 # MW DNPPE, g/mol
