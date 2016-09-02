# ContainerUV_VISTransPlots.R
#
# Purpose: Read, format, plot data from UV-VIS spectrophotometer transmittance
# profiles for various materials used for incubations in lipid & whole seawater
# community photo-oxidation experiments
#
# Created: 8/23/16 by James Collins, james.r.collins@aya.yale.edu
# Released under MIT License

# load in % transmittance data from separate raw data files
# these were acquired using a Thermo Evolution 300 UV-VIS spectrophotometer at
# 0.2 nm increments from 190 to 800 nm

# get current working directory so we can reset it later
current.wd = getwd()

# set wd to data location
setwd("data/raw/Evolution_300/container_transmittance_profiles") 

# read in data, concatenate into single data frame
quartz = read.csv("Quartz_%T_190_to_800.csv") # one wall of a 50 mL fused quartz vial
colnames(quartz) = c("lambda","percent_transmittance_quartz")
borosilicate = read.csv("Borosilicate_%T_190_to_800.csv") # one wall of a borosilicate ("regular" lab glass) EPA vial 
colnames(borosilicate) = c("lambda","percent_transmittance_borosilicate")
tedlar_PVF = read.csv("Tedlar_%T_190_to_800_b.csv") # 2 mil thickness Tedlar (PVF), one wall of Tedlar VOC sampling bag
colnames(tedlar_PVF) = c("lambda","percent_transmittance_tedlar_PVF")
PVF_PET = read.csv("Tedlar+mylar_%T_190_to_800_b.csv") # 4 mil thickness Mylar sheet plus Tedlar bag
colnames(PVF_PET) = c("lambda","percent_transmittance_PVF_PET")
mylar_PET_4mil = read.csv("Mylar_only_%T_190_to_800.csv") # 4 mil thickness Mylar (PET)
colnames(mylar_PET_4mil) = c("lambda","percent_transmittance_mylar_PET_4mil")
PC_bottle = read.csv("PC_bottle_%T_190_to_800.csv") # polycarbonate (PC) bottle, this one only at 1 nm resolution
colnames(PC_bottle) = c("lambda","percent_transmittance_PC_bottle")