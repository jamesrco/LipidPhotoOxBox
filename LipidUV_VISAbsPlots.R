# LipidUV_VISAbsPlots.R
#
# Purpose: Read, format, plot data from UV-VIS spectrophotometer absorbance
# profiles of the various lipids for the lipid photo-ox project; make some
# calculations of molar absorption ("extinction") coefficients and other properties
#
# Calculations of coefficients in both decadic and Napierian units
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

# load in "older" (fall 2016) data
# set wd to data location
setwd("/Users/jamesrco/Code/LipidPhotoOxBox/data/raw/Evolution_300/lipids_in_MeOH_initial") 

# read in data, concatenate into single data frame
# original profiles made at 1.113 mM concentration on 10 Aug 2016
PC22_6_init = read.csv("PC_22-6_in_MeOH_1.113mM.csv", 
                  stringsAsFactors = FALSE)
PC22_1_init = read.csv("PC_22-1_in_MeOH_1.113mM.csv", 
                        stringsAsFactors = FALSE)
PC22_0_init = read.csv("PC_22-0_in_MeOH,heated_1.113mM.csv", 
                  stringsAsFactors = FALSE)

# also, load in another 22:6 profile made on 26 Oct 2016 at much diluted 
# concentration (spec appeared to be saturated in original run)
# note that these were only collected at unit wavelength resolution

# **** these initial PC 22:6 data are highly suspect **** 
PC22_6_dil_init = read.csv("PC_22-6_in_MeOH_diluted_0.0027825mM.csv", 
                  stringsAsFactors = FALSE)
PC22_6_dil_init = PC22_6_dil_init[-c(1:8),]

# convert fields to numeric

LipidAbsData_init = cbind(as.numeric(as.character(PC22_6_init[,1])),
                     as.numeric(as.character(PC22_6_init[,2])),
                     as.numeric(as.character(PC22_1_init[,2])),
                     as.numeric(as.character(PC22_0_init[,2])))

LipidAbsData_init = as.data.frame(LipidAbsData_init[-c(1:8),]) # get rid of some unneeded header info

colnames(LipidAbsData_init) = c("lambda_nm",
                           "abs_PC22_6",
                           "abs_PC22_1",
                           "abs_PC22_0")

# append diluted 22:6 data, requires jerry-rigging since collected at different interval

LipidAbsData_init$abs_PC22_6_dil_init = NA

for (i in 1:nrow(LipidAbsData_init)) {
  
  if (length(PC22_6_dil_init$scan026[PC22_6_dil_init$Batch==LipidAbsData_init$lambda_nm[i]])>0) {
    
    LipidAbsData_init$abs_PC22_6_dil_init[i] = 
      as.numeric(PC22_6_dil_init$scan026[PC22_6_dil_init$Batch==LipidAbsData_init$lambda_nm[i]])
    
  }
  
}
 
# calculate decadic molar absorption coefficients (epsilon, in per M per cm)

# specify some constants
pathlength_cm = 10 # pathlength of cuvette used, in cm (same for all dates)
molarity_mM_20160810 = 1.113 # molarity of lipids in MeOH (original conc. on 10 Aug 16)
molarity_mM_20161026 = 0.0027825 # molarity of diluted 22:6 sol'n evaluated on 26 Oct 16
molarity_M_20160810 = molarity_mM_20160810/1000 # molarity, in mol/L
molarity_M_20161026 = molarity_mM_20161026/1000 # molarity, in mol/L

# calculations
LipidAbsData_init$epsilon_M_cm_PC22_6 = LipidAbsData_init$abs_PC22_6/(molarity_M_20160810*pathlength_cm)
LipidAbsData_init$epsilon_M_cm_PC22_1 = LipidAbsData_init$abs_PC22_1/(molarity_M_20160810*pathlength_cm)
LipidAbsData_init$epsilon_M_cm_PC22_0 = LipidAbsData_init$abs_PC22_0/(molarity_M_20160810*pathlength_cm)
# require a small baseline correction for the second batch of 22:6 data
LipidAbsData_init$epsilon_M_cm_PC22_6_dil = LipidAbsData_init$abs_PC22_6_dil/(molarity_M_20161026*pathlength_cm)

# also, calculate some Napierian molar absorption coefficients (kappa, in per M per cm); requires scaling factor of ln(10)

LipidAbsData_init$kappa_M_cm_PC22_6 = LipidAbsData_init$epsilon_M_cm_PC22_6*log(10)
LipidAbsData_init$kappa_M_cm_PC22_1 = LipidAbsData_init$epsilon_M_cm_PC22_1*log(10)
LipidAbsData_init$kappa_M_cm_PC22_0 = LipidAbsData_init$epsilon_M_cm_PC22_0*log(10)
LipidAbsData_init$kappa_M_cm_PC22_6_dil = LipidAbsData_init$epsilon_M_cm_PC22_6_dil*log(10)

# reset the working directory
setwd(initial.wd)

# load in new (verification run) Jan 2017 data
# set wd to data location
setwd("/Users/jamesrco/Code/LipidPhotoOxBox/data/raw/Evolution_300/lipids_in_MeOH_rerun_Jan2017") 

# read in data, concatenate into single data frame
# verification run, including DHA, of new standards from Avanti and Cayman Chem run on 17 Jan 2017
# didn't run PC 22:0
# included a dilution series for PC 22-6
PC22_6_011387mM_Jan17 = read.csv("PC_22-6_in_MeOH_0.11387mM_second_run.csv", 
                       stringsAsFactors = FALSE)
PC22_6_01576055mM_Jan17 = read.csv("PC_22-6_in_MeOH_0.1576055mM.csv", 
                                 stringsAsFactors = FALSE)
PC22_6_005679mM_Jan17 = read.csv("PC_22-6_in_MeOH_0.05679mM.csv", 
                                   stringsAsFactors = FALSE)
PC22_6_00285mM_Jan17 = read.csv("PC_22-6_in_MeOH_0.0285mM.csv", 
                                   stringsAsFactors = FALSE)
PC22_1_Jan17 = read.csv("PC_22-1_in_MeOH_1.113mM.csv", 
                       stringsAsFactors = FALSE)
DHA_01032mM_Jan17 = read.csv("DHA_in_MeOH_0.1032mM.csv", 
                       stringsAsFactors = FALSE)

# convert fields to numeric

LipidAbsData_Jan17 = cbind(as.numeric(as.character(PC22_6_011387mM_Jan17[,1])),
                          as.numeric(as.character(PC22_6_011387mM_Jan17[,2])),
                          as.numeric(as.character(PC22_6_01576055mM_Jan17[,2])),
                          as.numeric(as.character(PC22_6_005679mM_Jan17[,2])),
                          as.numeric(as.character(PC22_6_00285mM_Jan17[,2])),
                          as.numeric(as.character(PC22_1_Jan17[,2])),
                          as.numeric(as.character(DHA_01032mM_Jan17[,2])))

LipidAbsData_Jan17 = as.data.frame(LipidAbsData_Jan17[-c(1:8),]) # get rid of some unneeded header info

colnames(LipidAbsData_Jan17) = c("lambda_nm",
                                "abs_PC22_6_011387mM",
                                "abs_PC22_6_01576055mM",
                                "abs_PC22_6_005679mM",
                                "abs_PC22_6_00285mM",
                                "abs_PC22_1_1113mM",
                                "abs_DHA_01032mM")

# calculate decadic molar absorption coefficients (epsilon, in per M per cm)

# calculations
# concentrations must be in mol/L
LipidAbsData_Jan17$epsilon_M_cm_PC22_6_011387mM = LipidAbsData_Jan17$abs_PC22_6_011387mM/((0.11387/1000)*pathlength_cm)
LipidAbsData_Jan17$epsilon_M_cm_PC22_6_01576055mM = LipidAbsData_Jan17$abs_PC22_6_01576055mM/((0.1576055/1000)*pathlength_cm)
LipidAbsData_Jan17$epsilon_M_cm_PC22_6_005679mM = LipidAbsData_Jan17$abs_PC22_6_005679mM/((0.05679/1000)*pathlength_cm)
LipidAbsData_Jan17$epsilon_M_cm_PC22_6_00285mM = LipidAbsData_Jan17$abs_PC22_6_00285mM/((0.0285/1000)*pathlength_cm)
LipidAbsData_Jan17$epsilon_M_cm_PC22_1_1113mM = LipidAbsData_Jan17$abs_PC22_1_1113mM/((1.113/1000)*pathlength_cm)
LipidAbsData_Jan17$epsilon_M_cm_DHA_01032mM = LipidAbsData_Jan17$abs_DHA_01032mM/((0.1032/1000)*pathlength_cm)

# also, calculate some Napierian molar absorption coefficients (kappa, in per M per cm); requires scaling factor of ln(10)
LipidAbsData_Jan17$kappa_M_cm_PC22_6_011387mM = LipidAbsData_Jan17$epsilon_M_cm_PC22_6_011387mM*log(10)
LipidAbsData_Jan17$kappa_M_cm_PC22_6_01576055mM = LipidAbsData_Jan17$epsilon_M_cm_PC22_6_01576055mM*log(10)
LipidAbsData_Jan17$kappa_M_cm_PC22_6_005679mM = LipidAbsData_Jan17$epsilon_M_cm_PC22_6_005679mM*log(10)
LipidAbsData_Jan17$kappa_M_cm_PC22_6_00285mM = LipidAbsData_Jan17$epsilon_M_cm_PC22_6_00285mM*log(10)
LipidAbsData_Jan17$kappa_M_cm_PC22_1_1113mM = LipidAbsData_Jan17$epsilon_M_cm_PC22_1_1113mM*log(10)
LipidAbsData_Jan17$kappa_M_cm_DHA_01032mM = LipidAbsData_Jan17$epsilon_M_cm_DHA_01032mM*log(10)

# make plots for liposome experiment chapter/paper

# first, using older fall 2016 data

# absorbance, 22:0 and 22:1

absPlotCol = hsv(c(0.2, 0.57, 0.95), 1, 1, 0.8) # define colors
absPlotLty = c("solid","dashed","dotdash") # define lty

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "PCLipidAbs_22-0,22-1.pdf",
    width = 8, height = 6, pointsize = 12,
    bg = "white")

par(mar=c(5,5,1,1))
plot(LipidAbsData_init$lambda_nm,LipidAbsData_init$abs_PC22_0,"l",
     col = absPlotCol[1], lty = absPlotLty[1], lwd = "1.5",
     ylim = c(0,1.75), xlim = c(225,500), ylab = "Decadic absorbance",
     xlab = "Wavelength (nm)")
lines(LipidAbsData_init$lambda_nm,LipidAbsData_init$abs_PC22_1,
      col = absPlotCol[2], lty = absPlotLty[2], lwd = "1.5")

legend(x = 450, y = 1, bty = "n",
       legend = c("22:0/22:0 PC","22:1/22:1 PC","22:6/22:6 PC (1:400 dilution"),
       col = absPlotCol, lty = absPlotLty, lwd = 2)

dev.off()

# absorbance, 22:6

absPlotCol = hsv(c(0.2, 0.57, 0.95), 1, 1, 0.8) # define colors
absPlotLty = c("solid","dashed","dotdash") # define lty

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "PCLipidAbs_22-6.pdf",
    width = 8, height = 6, pointsize = 12,
    bg = "white")

par(mar=c(5,5,1,1))
plot(LipidAbsData_init$lambda_nm[seq(1,nrow(LipidAbsData_init),5)],
     LipidAbsData_init$abs_PC22_6_dil[seq(1,nrow(LipidAbsData_init),5)],"l",
     col = absPlotCol[3], lty = absPlotLty[3], lwd = "1.5",
     ylim = c(0,0.75), xlim = c(225,500), ylab = "Decadic absorbance",
     xlab = "Wavelength (nm)")

legend(x = 450, y = 2, bty = "n",
       legend = c("22:0/22:0 PC","22:1/22:1 PC","22:6/22:6 PC (1:400 dilution"),
       col = absPlotCol, lty = absPlotLty, lwd = 2)

dev.off()

# make plot for inset (290-315 nm)

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "PCLipidAbs_inset_290_315.pdf",
    width = 4, height = 4, pointsize = 12,
    bg = "white")

par(mar=c(5,5,1,1))
plot(LipidAbsData_init$lambda_nm,LipidAbsData_init$abs_PC22_0,"l",
     col = absPlotCol[1], lty = absPlotLty[1], lwd = "1.5",
     ylim = c(0,3.5), xlim = c(290,315), ylab = "Decadic absorbance",
     xlab = "Wavelength (nm)")
lines(LipidAbsData_init$lambda_nm,LipidAbsData_init$abs_PC22_1,
      col = absPlotCol[2], lty = absPlotLty[2], lwd = "1.5")
lines(LipidAbsData_init$lambda_nm[seq(1,nrow(LipidAbsData_init),5)],
      LipidAbsData_init$abs_PC22_6_dil[seq(1,nrow(LipidAbsData_init),5)],
      col = absPlotCol[3], lty = absPlotLty[3], lwd = "1.5")

dev.off()

# molar absorption coefficients
# plot of decadic values

absPlotCol = hsv(c(0.2, 0.57, 0.95), 1, 1, 0.8) # define colors
absPlotLty = c("solid","dashed","dotdash") # define lty

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "PCLipidMolExtCoeff_decadic.pdf",
    width = 8, height = 6, pointsize = 12,
    bg = "white")

par(mar=c(5,5,1,1))
plot(LipidAbsData_init$lambda_nm,log(LipidAbsData_init$epsilon_M_cm_PC22_0),"l",
     col = absPlotCol[1], lty = absPlotLty[1], lwd = "1.5",
     ylim = c(0,10), xlim = c(225,500),
     ylab = expression(paste("log ",epsilon[i]," (",M^-1," ",cm^-1,")")),
     xlab = "Wavelength (nm)")
lines(LipidAbsData_init$lambda_nm,log(LipidAbsData_init$epsilon_M_cm_PC22_1),
      col = absPlotCol[2], lty = absPlotLty[2], lwd = "1.5")
lines(LipidAbsData_init$lambda_nm[seq(1,nrow(LipidAbsData_init),5)],
      log(LipidAbsData_init$epsilon_M_cm_PC22_6_dil[seq(1,nrow(LipidAbsData_init),5)]),
      col = absPlotCol[3], lty = absPlotLty[3], lwd = "1.5")

legend(x = 400, y = 8, bty = "n",
       legend = c("22:0/22:0 PC","22:1/22:1 PC","22:6/22:6 PC"),
       col = absPlotCol, lty = absPlotLty, lwd = 2)

dev.off()

# molar absorption coefficients
# plot of Napierian values

absPlotCol = hsv(c(0.2, 0.57, 0.95), 1, 1, 0.8) # define colors
absPlotLty = c("solid","dashed","dotdash") # define lty

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "PCLipidMolExtCoeff_Napierian.pdf",
    width = 8, height = 6, pointsize = 12,
    bg = "white")

par(mar=c(5,5,1,1))
plot(LipidAbsData_init$lambda_nm,log(LipidAbsData_init$kappa_M_cm_PC22_0),"l",
     col = absPlotCol[1], lty = absPlotLty[1], lwd = "1.5",
     ylim = c(0,10), xlim = c(225,500),
     ylab = expression(paste("log ",kappa[i]," (",M^-1," ",cm^-1,")")),
     xlab = "Wavelength (nm)")
lines(LipidAbsData_init$lambda_nm,log(LipidAbsData_init$kappa_M_cm_PC22_1),
      col = absPlotCol[2], lty = absPlotLty[2], lwd = "1.5")
lines(LipidAbsData_init$lambda_nm[seq(1,nrow(LipidAbsData_init),5)],
      log(LipidAbsData_init$kappa_M_cm_PC22_6_dil[seq(1,nrow(LipidAbsData_init),5)]),
      col = absPlotCol[3], lty = absPlotLty[3], lwd = "1.5")

legend(x = 400, y = 8, bty = "n",
       legend = c("22:0/22:0 PC","22:1/22:1 PC","22:6/22:6 PC"),
       col = absPlotCol, lty = absPlotLty, lwd = 2)

dev.off()

# now, using Jan 2017 data
# will still have to use PC 22:0 data from the fall 2016 run since we didn't re-run PC 22:0

# absorbance, 22:0 and 22:1

absPlotCol = hsv(c(0.1, 0.35, 0.6, 0.85), 1, 1, 0.8) # define colors
absPlotLty = c("solid","dashed","dotdash","dotted") # define lty

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "PCLipidAbs_22-0,22-1.pdf",
    width = 8, height = 6, pointsize = 12,
    bg = "white")

par(mar=c(5,5,1,1))
plot(LipidAbsData_init$lambda_nm,LipidAbsData_init$abs_PC22_0,"l",
     col = absPlotCol[1], lty = absPlotLty[1], lwd = "1.5",
     ylim = c(0,1.75), xlim = c(225,500), ylab = "Decadic absorbance",
     xlab = "Wavelength (nm)")
lines(LipidAbsData_init$lambda_nm,LipidAbsData_init$abs_PC22_1,
      col = absPlotCol[2], lty = absPlotLty[2], lwd = "1.5")

legend(x = 450, y = 1, bty = "n",
       legend = c("22:0/22:0 PC","22:1/22:1 PC","22:6/22:6 PC","DHA"),
       col = absPlotCol, lty = absPlotLty, lwd = 2)

dev.off()

# absorbance, 22:6 and DHA together (new measurements)

# abs_PC22_6_011387mM
# abs_PC22_6_01576055mM
# abs_PC22_6_005679mM
# abs_PC22_6_00285mM
# abs_PC22_1_1113mM
# abs_DHA_01032mM

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "PCLipidAbs_22-6_DHA.pdf",
    width = 8, height = 6, pointsize = 12,
    bg = "white")

par(mar=c(5,5,1,1))
plot(LipidAbsData_Jan17$lambda_nm,LipidAbsData_Jan17$abs_PC22_6_011387mM,"l",
     col = absPlotCol[3], lty = absPlotLty[3], lwd = "1.5",
     ylim = c(0,1.1), xlim = c(225,500), ylab = "Decadic absorbance",
     xlab = "Wavelength (nm)")
lines(LipidAbsData_Jan17$lambda_nm,LipidAbsData_Jan17$abs_DHA_01032mM,
      col = absPlotCol[4], lty = absPlotLty[4], lwd = "1.5")

legend(x = 450, y = 1, bty = "n",
       legend = c("22:0/22:0 PC","22:1/22:1 PC","22:6/22:6 PC","DHA"),
       col = absPlotCol, lty = absPlotLty, lwd = 2)

dev.off()

# # absorbance, DHA
# 
# # abs_PC22_6_011387mM
# # abs_PC22_6_01576055mM
# # abs_PC22_6_005679mM
# # abs_PC22_6_00285mM
# # abs_PC22_1_1113mM
# # abs_DHA_01032mM
# 
# par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this
# 
# pdf(file = "PCLipidAbs_DHA.pdf",
#     width = 8, height = 6, pointsize = 12,
#     bg = "white")
# 
# par(mar=c(5,5,1,1))
# plot(LipidAbsData_Jan17$lambda_nm,LipidAbsData_Jan17$abs_DHA_01032mM,"l",
#      col = absPlotCol[3], lty = absPlotLty[3], lwd = "1.5",
#      ylim = c(0,1.1), xlim = c(225,500), ylab = "Absorbance",
#      xlab = "Wavelength (nm)")
# 
# dev.off()


# decadic molar absorption coefficients

# epsilon_M_cm_PC22_6_011387mM
# epsilon_M_cm_PC22_6_01576055mM
# epsilon_M_cm_PC22_6_005679mM
# epsilon_M_cm_PC22_6_00285mM
# epsilon_M_cm_PC22_1_1113mM
# epsilon_M_cm_DHA_01032mM

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "PCLipidMolExtCoeff_with_DHA_decadic.pdf",
    width = 8, height = 6, pointsize = 12,
    bg = "white")

par(mar=c(5,5,1,1))
plot(LipidAbsData_init$lambda_nm,log(LipidAbsData_init$epsilon_M_cm_PC22_0),"l",
     col = absPlotCol[1], lty = absPlotLty[1], lwd = "1.5",
     ylim = c(0,7), xlim = c(225,500),
     ylab = expression(paste("log ",epsilon[i]," (",M^-1," ",cm^-1,")")),
     xlab = "Wavelength (nm)")
lines(LipidAbsData_init$lambda_nm,log(LipidAbsData_init$epsilon_M_cm_PC22_1),
      col = absPlotCol[2], lty = absPlotLty[2], lwd = "1.5")
 lines(LipidAbsData_Jan17$lambda_nm,log(LipidAbsData_Jan17$epsilon_M_cm_PC22_6_011387mM),
       col = absPlotCol[3], lty = absPlotLty[3], lwd = "1.5")
lines(LipidAbsData_Jan17$lambda_nm,log(LipidAbsData_Jan17$epsilon_M_cm_DHA_01032mM),
      col = absPlotCol[4], lty = absPlotLty[4], lwd = "1.5")
legend(x = 350, y = 6, bty = "n",
       legend = c("22:0/22:0 PC","22:1/22:1 PC","22:6/22:6 PC","DHA"),
       col = absPlotCol, lty = absPlotLty, lwd = 2)

dev.off()

# Napierian molar absorption coefficients

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "PCLipidMolExtCoeff_with_DHA_Napierian.pdf",
    width = 8, height = 6, pointsize = 12,
    bg = "white")

par(mar=c(5,5,1,1))
plot(LipidAbsData_init$lambda_nm,log(LipidAbsData_init$kappa_M_cm_PC22_0),"l",
     col = absPlotCol[1], lty = absPlotLty[1], lwd = "1.5",
     ylim = c(0,7), xlim = c(225,500),
     ylab = expression(paste("log ",kappa[i]," (",M^-1," ",cm^-1,")")),
     xlab = "Wavelength (nm)")
lines(LipidAbsData_init$lambda_nm,log(LipidAbsData_init$kappa_M_cm_PC22_1),
      col = absPlotCol[2], lty = absPlotLty[2], lwd = "1.5")
lines(LipidAbsData_Jan17$lambda_nm,log(LipidAbsData_Jan17$kappa_M_cm_PC22_6_011387mM),
      col = absPlotCol[3], lty = absPlotLty[3], lwd = "1.5")
lines(LipidAbsData_Jan17$lambda_nm,log(LipidAbsData_Jan17$kappa_M_cm_DHA_01032mM),
      col = absPlotCol[4], lty = absPlotLty[4], lwd = "1.5")
legend(x = 350, y = 6, bty = "n",
       legend = c("22:0/22:0 PC","22:1/22:1 PC","22:6/22:6 PC","DHA"),
       col = absPlotCol, lty = absPlotLty, lwd = 2)

dev.off()

# save a data object containing the combined, validated absorption results reported in thesis/manuscript
LipidAbsData =
  LipidAbsData_Jan17[,-c(2:7,9:12,15:18)]
LipidAbsData$epsilon_M_cm_PC22_1 = LipidAbsData_init$epsilon_M_cm_PC22_1[1:1551]
LipidAbsData$epsilon_M_cm_PC22_0 = LipidAbsData_init$epsilon_M_cm_PC22_0[1:1551]
LipidAbsData$kappa_M_cm_PC22_1 = LipidAbsData_init$kappa_M_cm_PC22_1[1:1551]
LipidAbsData$kappa_M_cm_PC22_0 = LipidAbsData_init$kappa_M_cm_PC22_0[1:1551]

colnames(LipidAbsData)[2:5] = c("epsilon_M_cm_PC22_6","epsilon_M_cm_DHA",
                                    "kappa_M_cm_PC22_6","kappa_M_cm_DHA")

# some reordering

LipidAbsData = LipidAbsData[,c(1:3,6,7,4,5,8,9)]

save(LipidAbsData,file = "/Users/jamesrco/Code/LipidPhotoOxBox/data/nice/container_and_lipid_absorbances/PC_lipid_UV-VIS_molar_absorption_data.RData")

# some data

# wavelength of 22:6 absorbance maximum
  LipidAbsData_Jan17$lambda_nm[LipidAbsData_Jan17$abs_PC22_6_011387mM==max(LipidAbsData_Jan17$abs_PC22_6_011387mM, na.rm =T)]