# LipidUV_VISAbsPlots.R
#
# Purpose: Read, format, plot data from UV-VIS spectrophotometer absorbance
# profiles of the various lipids for the lipid photo-ox project; make some
# calculations of decadic molar extinction coefficients and other properties
#
# Created: 9/6/2016 by James Collins, james.r.collins@aya.yale.edu
# Released under MIT License

# load in % transmittance data from separate raw data files
# these were acquired using a Thermo Evolution 300 UV-VIS spectrophotometer at
# 0.2 nm increments from 190 to 800 nm

# assuming initial wd is that of the LipidPhotoOxBox repository, get current
# working directory so we can reset it later

# get current working directory so we can reset it later
initial.wd = getwd()

# set wd to data location
setwd("data/raw/Evolution_300/lipids_in_MeOH") 

# read in data, concatenate into single data frame
PC22_6 = read.csv("PC_22-6_in_MeOH.csv", 
                  stringsAsFactors = FALSE)
PC22_1 = read.csv("PC_22-1_in_MeOH.csv", 
                        stringsAsFactors = FALSE)
PC22_0 = read.csv("PC_22-0_in_MeOH,heated.csv", 
                  stringsAsFactors = FALSE)

# convert fields to numeric

LipidAbsData = cbind(as.numeric(as.character(PC22_6[,1])),
                     as.numeric(as.character(PC22_6[,2])),
                     as.numeric(as.character(PC22_1[,2])),
                     as.numeric(as.character(PC22_0[,2])))

LipidAbsData = as.data.frame(LipidAbsData[-c(1:8),]) # get rid of some unneeded header info

colnames(LipidAbsData) = c("lambda_nm",
                           "abs_PC22_6",
                           "abs_PC22_1",
                           "abs_PC22_0")

# calculate decadic molar extinction coefficients (epsilon, in per M per cm)

# specify some constants
pathlength_cm = 10 # pathlength of cuvette used, in cm
molarity_mM = 1.113 # molarity of lipids in MeOH (all 3 measured at same molarity)
molarity_M = molarity_mM/1000 # molarity, in mol/L

# calculations
LipidAbsData$epsilon_M_cm_PC22_6 = LipidAbsData$abs_PC22_6/(molarity_M*pathlength_cm)
LipidAbsData$epsilon_M_cm_PC22_1 = LipidAbsData$abs_PC22_1/(molarity_M*pathlength_cm)
LipidAbsData$epsilon_M_cm_PC22_0 = LipidAbsData$abs_PC22_0/(molarity_M*pathlength_cm)

# save the data file (we'll need it later)

save(LipidAbsData, file = "PC_lipid_abs_data.RData")

# make plots for liposome experiment chapter

# absorbance

absPlotCol = hsv(c(0.2, 0.57, 0.95), 1, 1, 0.8) # define colors
absPlotLty = c("solid","dashed","dotdash") # define lty

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "PCLipidAbs.pdf",
    width = 8, height = 6, pointsize = 12,
    bg = "white")

par(mar=c(5,5,1,1))
plot(LipidAbsData$lambda_nm,LipidAbsData$abs_PC22_0,"l",
     col = absPlotCol[1], lty = absPlotLty[1], lwd = "1.5",
     ylim = c(0,3.5), xlim = c(225,700), ylab = "Absorbance",
     xlab = "Wavelength (nm)")
lines(LipidAbsData$lambda_nm,LipidAbsData$abs_PC22_1,
      col = absPlotCol[2], lty = absPlotLty[2], lwd = "1.5")
lines(LipidAbsData$lambda_nm,LipidAbsData$abs_PC22_6,
      col = absPlotCol[3], lty = absPlotLty[3], lwd = "1.5")

legend(x = 450, y = 2, bty = "n",
       legend = c("22:0/22:0 PC","22:1/22:1 PC","22:6/22:6 PC"),
       col = absPlotCol, lty = absPlotLty, lwd = 2)

dev.off()

# make plot for inset (290-315 nm)

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "PCLipidAbs_inset_290_315.pdf",
    width = 4, height = 4, pointsize = 12,
    bg = "white")

par(mar=c(5,5,1,1))
plot(LipidAbsData$lambda_nm,LipidAbsData$abs_PC22_0,"l",
     col = absPlotCol[1], lty = absPlotLty[1], lwd = "1.5",
     ylim = c(0,3.5), xlim = c(290,315), ylab = "Absorbance",
     xlab = "Wavelength (nm)")
lines(LipidAbsData$lambda_nm,LipidAbsData$abs_PC22_1,
      col = absPlotCol[2], lty = absPlotLty[2], lwd = "1.5")
lines(LipidAbsData$lambda_nm,LipidAbsData$abs_PC22_6,
      col = absPlotCol[3], lty = absPlotLty[3], lwd = "1.5")

dev.off()

# decadic molar extinction coefficients

absPlotCol = hsv(c(0.2, 0.57, 0.95), 1, 1, 0.8) # define colors
absPlotLty = c("solid","dashed","dotdash") # define lty

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "PCLipidMolExtCoeff.pdf",
    width = 8, height = 6, pointsize = 12,
    bg = "white")

par(mar=c(5,5,1,1))
plot(LipidAbsData$lambda_nm,log(LipidAbsData$epsilon_M_cm_PC22_0),"l",
     col = absPlotCol[1], lty = absPlotLty[1], lwd = "1.5",
     ylim = c(0,6), xlim = c(225,500),
     ylab = expression(paste("log ",epsilon[i]," (",M^-1," ",cm^-1,")")),
     xlab = "Wavelength (nm)")
lines(LipidAbsData$lambda_nm,log(LipidAbsData$epsilon_M_cm_PC22_1),
      col = absPlotCol[2], lty = absPlotLty[2], lwd = "1.5")
lines(LipidAbsData$lambda_nm,log(LipidAbsData$epsilon_M_cm_PC22_6),
      col = absPlotCol[3], lty = absPlotLty[3], lwd = "1.5")

legend(x = 400, y = 4, bty = "n",
       legend = c("22:0/22:0 PC","22:1/22:1 PC","22:6/22:6 PC"),
       col = absPlotCol, lty = absPlotLty, lwd = 2)

dev.off()

# some data

# wavelength of 22:6 absorbance maximum
LipidAbsData$lambda_nm[LipidAbsData$abs_PC22_6==max(LipidAbsData$abs_PC22_6)]