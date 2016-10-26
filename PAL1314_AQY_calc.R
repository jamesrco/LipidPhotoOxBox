# PAL1314_AQY_calc.R

# Purpose: Calculate apparent quantum yields for photolysis of PUFA-containing
# PC species examined in the PAL1314 liposome experiments 

# ****** Assumes some variables created by PAL1314_liposome_expts.R and
# UV_TS_analysis_PAL1314.R are already in user's workspace; these files can be
# found at https://github.com/jamesrco/LipidPhotoOxBox

# Created 10/26/16 by J.R.C.

# first, pull out necessary change in concentration data for PC 22:6, 22:6 from
# Exp 13

Exp_13_PC.samp.pmol.mL.norm.mean