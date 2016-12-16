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

# assuming initial wd is that of the LipidPhotoOxBox repository, get current
# working directory so we can reset it later

# get current working directory so we can reset it later
initial.wd = getwd()

# set wd to data location
setwd("/Users/jrcollins/Code/LipidPhotoOxBox/data/raw/Evolution_300/container_transmittance_profiles") 

# read in data, concatenate into single data frame
quartz = read.csv("Quartz_%T_190_to_800.csv", 
                  stringsAsFactors = FALSE) # one wall of a 50 mL fused quartz vial
borosilicate = read.csv("Borosilicate_%T_190_to_800.csv", 
                        stringsAsFactors = FALSE) # one wall of a borosilicate ("regular" lab glass) EPA vial 
tedlar_PVF = read.csv("Tedlar_%T_190_to_800_b.csv", 
                      stringsAsFactors = FALSE) # 2 mil thickness Tedlar (PVF), one wall of Tedlar VOC sampling bag
PVF_PET = read.csv("Tedlar+mylar_%T_190_to_800_b.csv", 
                   stringsAsFactors = FALSE) # 4 mil thickness Mylar sheet plus Tedlar bag
mylar_PET_4mil = read.csv("Mylar_only_%T_190_to_800.csv", 
                          stringsAsFactors = FALSE) # 4 mil thickness Mylar (PET)
PC_bottle = read.csv("PC_bottle_%T_190_to_800.csv", 
                     stringsAsFactors = FALSE) # polycarbonate (PC) bottle, this one only at 1 nm resolution

PctTransData = cbind(as.numeric(as.character(quartz[,1])),
                     as.numeric(as.character(quartz[,2])),
                     as.numeric(as.character(borosilicate[,2])),
                     as.numeric(as.character(tedlar_PVF[,2])),
                     as.numeric(as.character(PVF_PET[,2])),
                     as.numeric(as.character(mylar_PET_4mil[,2])))

PctTransData = as.data.frame(PctTransData[-c(1:8),]) # get rid of some unneeded header info

# convert fields to numeric

colnames(PctTransData) = c("lambda_nm",
                           "transmittance_quartz_pct",
                           "transmittance_borosilicate_pct",
                           "transmittance_tedlar_PVF_pct",
                           "transmittance_PVF_PET_pct",
                           "transmittance_mylar_PET_4mil_pct")

# save the data file (we'll need it later)

save(PctTransData, file = "Container_percent_transmittance.RData")

# make main plot for liposome experiment chapter

transPlotCol = c("darkblue","skyblue","darkred","coral1") # define colors
transPlotLty = c("solid","dashed","dotted","dotdash") # define lty

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "GlassUV_VISTrans.pdf",
    width = 8, height = 6, pointsize = 12,
    bg = "white")

par(mar=c(5,5,1,1))
plot(PctTransData$lambda_nm,PctTransData$transmittance_quartz_pct,"l",
     col = transPlotCol[1], lty = transPlotLty[1], lwd = "1",
     ylim = c(0,100), xlim = c(200,700), ylab = "Percent transmittance",
     xlab = "Wavelength (nm)")
lines(PctTransData$lambda_nm,PctTransData$transmittance_borosilicate_pct,
      col = transPlotCol[2], lty = transPlotLty[2], lwd = "2")
# lines(PctTransData$lambda_nm,PctTransData$transmittance_tedlar_PVF_pct,
#       col = transPlotCol[3], lty = transPlotLty[3], lwd = "2")
# lines(PctTransData$lambda_nm,PctTransData$transmittance_PVF_PET_pct,
#       col = transPlotCol[4], lty = transPlotLty[4], lwd = "2")

legend(x = 450, y = 40, bty = "n",
#        legend = c("Quartz glass vial","Borosilicate glass vial","PVF incubation bag",
#          "PVF bag w/PET screen"),
       legend = c("Quartz glass vial","Borosilicate glass vial"),
       col = transPlotCol, lty = transPlotLty, lwd = 2)

dev.off()

# make plot for inset (290-315 nm)

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "GlassUV_VISTrans_inset_290_315.pdf",
    width = 4, height = 4, pointsize = 12,
    bg = "white")

par(mar=c(5,5,1,1))
plot(PctTransData$lambda_nm,PctTransData$transmittance_quartz_pct,"l",
     col = transPlotCol[1], lty = transPlotLty[1], lwd = "1",
     ylim = c(0,100), xlim = c(290,315), ylab = "Percent transmittance",
     xlab = "Wavelength (nm)")
lines(PctTransData$lambda_nm,PctTransData$transmittance_borosilicate_pct,
      col = transPlotCol[2], lty = transPlotLty[2], lwd = "2")
# lines(PctTransData$lambda_nm,PctTransData$transmittance_tedlar_PVF_pct,
#       col = transPlotCol[3], lty = transPlotLty[3], lwd = "2")
# lines(PctTransData$lambda_nm,PctTransData$transmittance_PVF_PET_pct,
#       col = transPlotCol[4], lty = transPlotLty[4], lwd = "2")

dev.off()

# make main plot for thesis Appendix B

transPlotCol = c("darkblue","skyblue","darkred","coral1") # define colors
transPlotLty = c("solid","dashed","dotted","dotdash") # define lty

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "Glass_bags_UV_VISTrans.pdf",
    width = 8, height = 6, pointsize = 12,
    bg = "white")

par(mar=c(5,5,1,1))
plot(PctTransData$lambda_nm,PctTransData$transmittance_quartz_pct,"l",
     col = transPlotCol[1], lty = transPlotLty[1], lwd = "1",
     ylim = c(0,100), xlim = c(200,700), ylab = "Percent transmittance",
     xlab = "Wavelength (nm)")
lines(PctTransData$lambda_nm,PctTransData$transmittance_borosilicate_pct,
      col = transPlotCol[2], lty = transPlotLty[2], lwd = "2")
 lines(PctTransData$lambda_nm,PctTransData$transmittance_tedlar_PVF_pct,
       col = transPlotCol[3], lty = transPlotLty[3], lwd = "2")
 lines(PctTransData$lambda_nm,PctTransData$transmittance_PVF_PET_pct,
       col = transPlotCol[4], lty = transPlotLty[4], lwd = "2")

legend(x = 450, y = 40, bty = "n",
               legend = c("Quartz glass vial","Borosilicate glass vial","PVF incubation bag",
                 "PVF bag w/PET screen"),
       legend = c("Quartz glass vial","Borosilicate glass vial"),
       col = transPlotCol, lty = transPlotLty, lwd = 2)

dev.off()

# make plot for inset (290-315 nm)

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "Glass_bags_UV_VISTrans_inset_290_315.pdf",
    width = 4, height = 4, pointsize = 12,
    bg = "white")

par(mar=c(5,5,1,1))
plot(PctTransData$lambda_nm,PctTransData$transmittance_quartz_pct,"l",
     col = transPlotCol[1], lty = transPlotLty[1], lwd = "1",
     ylim = c(0,100), xlim = c(290,315), ylab = "Percent transmittance",
     xlab = "Wavelength (nm)")
lines(PctTransData$lambda_nm,PctTransData$transmittance_borosilicate_pct,
      col = transPlotCol[2], lty = transPlotLty[2], lwd = "2")
 lines(PctTransData$lambda_nm,PctTransData$transmittance_tedlar_PVF_pct,
       col = transPlotCol[3], lty = transPlotLty[3], lwd = "2")
 lines(PctTransData$lambda_nm,PctTransData$transmittance_PVF_PET_pct,
       col = transPlotCol[4], lty = transPlotLty[4], lwd = "2")

dev.off()