# PAL1314_AQY_calc.R

# Purpose: Calculate apparent quantum yields for photolysis of PUFA-containing
# PC species examined in the PAL1314 liposome experiments 

# Formulae have been updated since time of thesis publication; now consistent
# with those found in manuscript in prep for submission

# ****** Assumes some variables created by PAL1314_liposome_expts.R,
# LipidUV_VISAbsPlots.R, ContainerUV_VISTransPlots.R, and
# UV_TS_analysis_PAL1314.R are already in user's workspace; these files can be
# found at https://github.com/jamesrco/LipidPhotoOxBox

# Created 10/26/16 by J.R.C.

# set working directory to parent LipidPhotoOxBox repo
setwd("/Users/jamesrco/Code/LipidPhotoOxBox")

##### seawater contribution to absorbance ##### 

# first, pull in and process some spectrophotometric profiles of seawater from Palmer
# need these to calculate the contribution to absorbance of the seawater
# matrix over the actionable wavelength range

# set wd to data location
setwd("data/raw/Evolution_300/PAL1314_LMG1401_UV_VIS_seawater_absorbances") 

# read in absorbance data
PAL1314_LMG1401_UV_VIS_SW_abs = read.csv("PAL1314_LMG1401_UV_VIS_seawater_absorbance_Evolution_300.csv", 
                           stringsAsFactors = FALSE, skip = 9)

# exploratory plot

color_ramp.SW.abs = colorRampPalette(c("green", "red")) # create a color ramp

SW.abs.plot.colors = rev(color_ramp.SW.abs(ncol(PAL1314_LMG1401_UV_VIS_SW_abs)))

plot(PAL1314_LMG1401_UV_VIS_SW_abs$Wavelength_nm,
     PAL1314_LMG1401_UV_VIS_SW_abs[,2],"l",
     col = SW.abs.plot.colors[1], lty = 1, lwd = "0.5",
     ylim = c(0,0.2), xlim = c(290,720),
     ylab = "Absorbance",
     xlab = "Wavelength (nm)")

text(453,
     PAL1314_LMG1401_UV_VIS_SW_abs[25,2],
     labels = colnames(PAL1314_LMG1401_UV_VIS_SW_abs)[2],
     offset = -.25, pos = 3, cex = 0.5)

# overlay other depths

for (i in 3:ncol(PAL1314_LMG1401_UV_VIS_SW_abs)) {
  
  lines(PAL1314_LMG1401_UV_VIS_SW_abs$Wavelength_nm,
        PAL1314_LMG1401_UV_VIS_SW_abs[,i],
        col = SW.abs.plot.colors[i], lty = 1, lwd = "0.5")
  
#  if (i %in% c(7,9,11,13,15,16,17)) {
    text(453,
         PAL1314_LMG1401_UV_VIS_SW_abs[25,i],
         labels = colnames(PAL1314_LMG1401_UV_VIS_SW_abs)[i],
         offset = -.25, pos = 3, cex = 0.5)
    
#  }
  
}

# some more data visualization, using just the PAL1314 data (not LMG1401)

# load in some metadata

PAL1314_DOC_log_metadata = read.csv("PAL1314_DOC_log_extract.csv", 
                                         stringsAsFactors = FALSE, skip = 1, header = T)


plot(PAL1314_LMG1401_UV_VIS_SW_abs$Wavelength_nm,
     PAL1314_LMG1401_UV_VIS_SW_abs[,2],
     "l", col = SW.abs.plot.colors[36], lwd = 0.5,
     lty = 1,
     ylim = c(-.05,0.6),
     xlim = c(230,395))

# overlay others, by color according to depth and date

for (i in c(3:14)) {
  
  depth = as.character(PAL1314_DOC_log_metadata$Depth_m[PAL1314_DOC_log_metadata$Sample_ID==colnames(PAL1314_LMG1401_UV_VIS_SW_abs)[i]])
  stn = PAL1314_DOC_log_metadata$PAL_LTER_Station_ID[PAL1314_DOC_log_metadata$Sample_ID==colnames(PAL1314_LMG1401_UV_VIS_SW_abs)[i]]
  date = PAL1314_DOC_log_metadata$Date[PAL1314_DOC_log_metadata$Sample_ID==colnames(PAL1314_LMG1401_UV_VIS_SW_abs)[i]]
  
  linetype = switch(depth,
                    "0" = 1,
                    "5" = 2,
                    "10" = 3,
                    "20" = 4)
  
  line.col = switch(date,
                     "12-Dec-13" = 1,
                     "24-Dec-13" = 25,
                    "27-Dec-13" = 50)
  
  lines(PAL1314_LMG1401_UV_VIS_SW_abs$Wavelength_nm,
        PAL1314_LMG1401_UV_VIS_SW_abs[,i],
        col = SW.abs.plot.colors[line.col], lty = linetype, lwd = "0.5")
  
  text(453,
       PAL1314_LMG1401_UV_VIS_SW_abs[25,i],
       labels = colnames(PAL1314_LMG1401_UV_VIS_SW_abs)[i],
       offset = -.25, pos = 3, cex = 0.5)
  
}

# it appears depth doesn't matter nearly so much as date (an earlier examination
# showed that the station (E vs B) didn't matter too much either)

# so, can calculate some averages based on date groupings, then display centroids
# and uncertainties; will convert to absorbtion coefficient (alpha, units of per m)

# first, define our pathlength

pathlength_cm.PAL1314_SW_abs = 10 # pathlength of cuvette used, in cm

# create a 12 Dec 13 subset and a 27 Dec 13 subset

PAL1314_UV_VIS_SW_abs_profile_means = as.data.frame(matrix(NA, nrow(PAL1314_LMG1401_UV_VIS_SW_abs),
                                                   5))

PAL1314_UV_VIS_SW_abs_profile_means[,1] =
  PAL1314_LMG1401_UV_VIS_SW_abs[,1]

colnames(PAL1314_UV_VIS_SW_abs_profile_means) = 
  c("Wavelength_nm","PAL1314_12Dec13_SW_alpha_per_m.mean",
    "PAL1314_12Dec13_SW_alpha_per_m.sd","PAL1314_27Dec13_SW_alpha_per_m.mean",
    "PAL1314_27Dec13_SW_alpha_per_m.sd")

# 12 Dec 13 calcs

PAL1314.SW.abs.subset.12Dec13 = PAL1314_LMG1401_UV_VIS_SW_abs[,colnames(PAL1314_LMG1401_UV_VIS_SW_abs) %in% PAL1314_DOC_log_metadata$Sample_ID[PAL1314_DOC_log_metadata$Date=="12-Dec-13"]]

PAL1314_UV_VIS_SW_abs_profile_means$PAL1314_12Dec13_SW_alpha_per_m.mean = apply(
  PAL1314.SW.abs.subset.12Dec13,1,mean)*log(10)/(pathlength_cm.PAL1314_SW_abs/100)

PAL1314_UV_VIS_SW_abs_profile_means$PAL1314_12Dec13_SW_alpha_per_m.sd = apply(
  PAL1314.SW.abs.subset.12Dec13,1,sd)*log(10)/(pathlength_cm.PAL1314_SW_abs/100)

# apply a baseline correction

PAL1314_UV_VIS_SW_abs_profile_means$PAL1314_12Dec13_SW_alpha_per_m.mean = 
  PAL1314_UV_VIS_SW_abs_profile_means$PAL1314_12Dec13_SW_alpha_per_m.mean+
  abs(min(PAL1314_UV_VIS_SW_abs_profile_means$PAL1314_12Dec13_SW_alpha_per_m.mean))

# 27 Dec 13 calcs

PAL1314.SW.abs.subset.27Dec13 = PAL1314_LMG1401_UV_VIS_SW_abs[,colnames(PAL1314_LMG1401_UV_VIS_SW_abs) %in% PAL1314_DOC_log_metadata$Sample_ID[PAL1314_DOC_log_metadata$Date %in% c("27-Dec-13")]]

PAL1314_UV_VIS_SW_abs_profile_means$PAL1314_27Dec13_SW_alpha_per_m.mean = apply(
  PAL1314.SW.abs.subset.27Dec13,1,mean)*log(10)/(pathlength_cm.PAL1314_SW_abs/100)

PAL1314_UV_VIS_SW_abs_profile_means$PAL1314_27Dec13_SW_alpha_per_m.sd = apply(
  PAL1314.SW.abs.subset.27Dec13,1,sd)*log(10)/(pathlength_cm.PAL1314_SW_abs/100)

# apply a baseline correction

PAL1314_UV_VIS_SW_abs_profile_means$PAL1314_27Dec13_SW_alpha_per_m.mean = 
  PAL1314_UV_VIS_SW_abs_profile_means$PAL1314_27Dec13_SW_alpha_per_m.mean+
  abs(min(PAL1314_UV_VIS_SW_abs_profile_means$PAL1314_27Dec13_SW_alpha_per_m.mean))

# plot mean abs from one of the two date groupings, to start

par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

pdf(file = "PAL1314_SW_alpha_means.pdf",
    width = 8, height = 6, pointsize = 12,
    bg = "white")

par(mar=c(5,5,1,5))

plot(PAL1314_UV_VIS_SW_abs_profile_means$Wavelength_nm,
     PAL1314_UV_VIS_SW_abs_profile_means$PAL1314_12Dec13_SW_alpha_per_m.mean,
     "l",
     ylim = c(-.05,15),
     xlim = c(230,395),
     xlab = "Wavelength (nm)",
     ylab = "Alpha of seawater(lambda) (m-1)")

# uncertainties (±sd), 12 Dec 13 data

polygon(c(PAL1314_UV_VIS_SW_abs_profile_means$Wavelength_nm,
          rev(PAL1314_UV_VIS_SW_abs_profile_means$Wavelength_nm)),
        c(PAL1314_UV_VIS_SW_abs_profile_means$PAL1314_12Dec13_SW_alpha_per_m.mean+
            PAL1314_UV_VIS_SW_abs_profile_means$PAL1314_12Dec13_SW_alpha_per_m.sd,
          rev(PAL1314_UV_VIS_SW_abs_profile_means$PAL1314_12Dec13_SW_alpha_per_m.mean-
                PAL1314_UV_VIS_SW_abs_profile_means$PAL1314_12Dec13_SW_alpha_per_m.sd)),
        border = NA, col = "lightcyan")

# uncertainties (±sd), 27 Dec 13 data

polygon(c(PAL1314_UV_VIS_SW_abs_profile_means$Wavelength_nm,
          rev(PAL1314_UV_VIS_SW_abs_profile_means$Wavelength_nm)),
        c(PAL1314_UV_VIS_SW_abs_profile_means$PAL1314_27Dec13_SW_alpha_per_m.mean+
            PAL1314_UV_VIS_SW_abs_profile_means$PAL1314_27Dec13_SW_alpha_per_m.sd,
          rev(PAL1314_UV_VIS_SW_abs_profile_means$PAL1314_27Dec13_SW_alpha_per_m.mean-
                PAL1314_UV_VIS_SW_abs_profile_means$PAL1314_27Dec13_SW_alpha_per_m.sd)),
        border = NA, col = "mistyrose")

# plot means again over top of the polygons

lines(PAL1314_UV_VIS_SW_abs_profile_means$Wavelength_nm,
      PAL1314_UV_VIS_SW_abs_profile_means$PAL1314_12Dec13_SW_alpha_per_m.mean,
      col = "turquoise3")

lines(PAL1314_UV_VIS_SW_abs_profile_means$Wavelength_nm,
      PAL1314_UV_VIS_SW_abs_profile_means$PAL1314_27Dec13_SW_alpha_per_m.mean,
      col = "red")

dev.off()

##### pulling in tranmissivity data for the quartz and borosilicate glass ##### 
# from object PctTransData, created by ContainerUV_VISTransPlots.R

FracTrans = PctTransData[,c(1:3)]
FracTrans[,c(2:3)] = FracTrans[,c(2:3)]/100 # change to fraction

##### d[22:6]/dt calculations, in pmol/L/hr, from Exp 13 ##### 

# first, calculate d[22:6]/dt (change in concentration) using data
# data for PC 22:6, 22:6 from Exp 13; will need a few figures — one for apparent
# dark oxidation (autooxidation component), one for UVB+UVA oxidation,
# one for UVB oxidation

# define our dt, in hrs
dt_Exp13_hr = 8.2

# preallocate (second element will be uncertainty)
d22_6_dt.Exp13.pmol_mL_hr.dark = vector(mode = "double", length = 2)
d22_6_dt.Exp13.pmol_mL_hr.noUVB = vector(mode = "double", length = 2)
d22_6_dt.Exp13.pmol_mL_hr.plusUVB = vector(mode = "double", length = 2)
d22_6_dt.Exp13.pmol_mL_hr.UVB.photo_ox = vector(mode = "double", length = 2)
d22_6_dt.Exp13.pmol_mL_hr.UVA_UVB.photo_ox = vector(mode = "double", length = 2)
d22_6_dt.Exp13.pmol_mL_hr.UVA.photo_ox = vector(mode = "double", length = 2)

# calculations
# dark
d22_6_dt.Exp13.pmol_mL_hr.dark[1] = 
  (Exp_13_PC.samp.pmol.mL.norm.mean[c("PC 44:12"),c("Dark_control_no_HB_2013-12-14 17:50:00.mean")]-
  Exp_13_PC.samp.pmol.mL.norm.mean[c("PC 44:12"),c("Dark_control_2013-12-14 09:30:00.mean")])/
  dt_Exp13_hr

d22_6_dt.Exp13.pmol_mL_hr.dark[2] = (sqrt((Exp_13_PC.samp.pmol.mL.norm.mean[c("PC 44:12"),c("Dark_control_no_HB_2013-12-14 17:50:00.se")])^2+
                                           (Exp_13_PC.samp.pmol.mL.norm.mean[c("PC 44:12"),c("Dark_control_2013-12-14 09:30:00.se")])^2))/
  dt_Exp13_hr

# -UVB
d22_6_dt.Exp13.pmol_mL_hr.noUVB[1] = 
  (Exp_13_PC.samp.pmol.mL.norm.mean[c("PC 44:12"),c("EPA_no_HB_2013-12-14 17:50:00.mean")]-
  Exp_13_PC.samp.pmol.mL.norm.mean[c("PC 44:12"),c("Dark_control_2013-12-14 09:30:00.mean")])/
  dt_Exp13_hr

d22_6_dt.Exp13.pmol_mL_hr.noUVB[2] = (sqrt((Exp_13_PC.samp.pmol.mL.norm.mean[c("PC 44:12"),c("EPA_no_HB_2013-12-14 17:50:00.se")])^2+
                                          (Exp_13_PC.samp.pmol.mL.norm.mean[c("PC 44:12"),c("Dark_control_2013-12-14 09:30:00.se")])^2))/
  dt_Exp13_hr

# +UVB
d22_6_dt.Exp13.pmol_mL_hr.plusUVB[1] = 
  (Exp_13_PC.samp.pmol.mL.norm.mean[c("PC 44:12"),c("Quartz_no_HB_2013-12-14 17:50:00.mean")]-
  Exp_13_PC.samp.pmol.mL.norm.mean[c("PC 44:12"),c("Dark_control_2013-12-14 09:30:00.mean")])/
  dt_Exp13_hr

d22_6_dt.Exp13.pmol_mL_hr.plusUVB[2] = (sqrt((Exp_13_PC.samp.pmol.mL.norm.mean[c("PC 44:12"),c("Quartz_no_HB_2013-12-14 17:50:00.se")])^2+
                                             (Exp_13_PC.samp.pmol.mL.norm.mean[c("PC 44:12"),c("Dark_control_2013-12-14 09:30:00.se")])^2))/
  dt_Exp13_hr

# now, define control-corrected rates loss due to all photooxidation (direct+indirect),
# direct photooxidation, and indirect photo-oxidation

d22_6_dt.Exp13.pmol_mL_hr.UVA_UVB.photo_ox[1] = d22_6_dt.Exp13.pmol_mL_hr.plusUVB[1]-
  d22_6_dt.Exp13.pmol_mL_hr.dark[1]
d22_6_dt.Exp13.pmol_mL_hr.UVA_UVB.photo_ox[2] = sqrt((d22_6_dt.Exp13.pmol_mL_hr.plusUVB[2])^2+
  (d22_6_dt.Exp13.pmol_mL_hr.dark[2])^2)

d22_6_dt.Exp13.pmol_mL_hr.UVA.photo_ox[1] = d22_6_dt.Exp13.pmol_mL_hr.noUVB[1]-
  d22_6_dt.Exp13.pmol_mL_hr.dark[1]
d22_6_dt.Exp13.pmol_mL_hr.UVA.photo_ox[2] = sqrt((d22_6_dt.Exp13.pmol_mL_hr.noUVB[2])^2+
                                                  (d22_6_dt.Exp13.pmol_mL_hr.dark[2])^2)

d22_6_dt.Exp13.pmol_mL_hr.UVB.photo_ox[1] = d22_6_dt.Exp13.pmol_mL_hr.UVA_UVB.photo_ox[1]-
  d22_6_dt.Exp13.pmol_mL_hr.UVA.photo_ox[1]
d22_6_dt.Exp13.pmol_mL_hr.UVB.photo_ox[2] = sqrt((d22_6_dt.Exp13.pmol_mL_hr.UVA_UVB.photo_ox[2])^2+
                                                       (d22_6_dt.Exp13.pmol_mL_hr.UVA.photo_ox[2])^2)

##### calculation of total radiation received at each wavelength during the experiment #####
# produces the E_n,p,sigma defined in Equation 7 in the manuscript

# i.e., the radiation received by the JAZ sensor, in the deck tank, 
# at the same depth as the vials, over course of Exp 13

# integrated figures in units of photons/cm2/wavelength

# extract subset of data for the experimental time interval
PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.sub = 
  PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec[
    PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec$Timestamp_GMT>=
      as.POSIXct('2013-12-14 12:30:00', tz = "GMT") &  # 9:30 local time
      PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec$Timestamp_GMT<=
      as.POSIXct('2013-12-14 20:50:00', tz = "GMT") # 17:50 local time
    ,]

# convert to units of photons/cm2, step by step
PAL1314_JAZ_subsurf_hires_full_spectrum_W_m2.14Dec.sub = 
  PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.sub[,5:(ncol(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.sub)-1)]/100

lambda_nm_JAZ = as.numeric(colnames(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.sub)[5:(ncol(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.sub)-1)])

PAL1314_JAZ_subsurf_hires_full_spectrum_umol_photons_m2_s.14Dec.sub =
  matrix(NA, ncol = length(lambda_nm_JAZ), nrow = nrow(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.14Dec.sub))

for (i in 1:ncol(PAL1314_JAZ_subsurf_hires_full_spectrum_umol_photons_m2_s.14Dec.sub)) {
  
  PAL1314_JAZ_subsurf_hires_full_spectrum_umol_photons_m2_s.14Dec.sub[,i] =
    PAL1314_JAZ_subsurf_hires_full_spectrum_W_m2.14Dec.sub[,i]*lambda_nm_JAZ[i]*0.836*(10^-2)
  
}

# integrate over time, by wavelength

# first, time integration at each wavelength

# need package caTools, if not loaded already

detach("package:RSEIS", unload=TRUE)
library(caTools)

# preallocate vector
PAL1314.JAZ.14Dec.E_n_p_sigma_umol_photons_m2 = vector(mode = "double",
                                                     length = ncol(PAL1314_JAZ_subsurf_hires_full_spectrum_umol_photons_m2_s.14Dec.sub))

for (i in 1:length(PAL1314.JAZ.14Dec.E_n_p_sigma_umol_photons_m2)) {
  
  # requires vector of times in second from UV_TS_analysis_PAL1314.R, with t = 0 being timepoint at beginning of experiment
  
  PAL1314.JAZ.14Dec.E_n_p_sigma_umol_photons_m2[i] =
    caTools::trapz(Exp13_timeint.s[1:length(Exp13_timeint.s)],PAL1314_JAZ_subsurf_hires_full_spectrum_umol_photons_m2_s.14Dec.sub[,i])
  
}

##### workup of data from Experiment 13 ##### 

# define a few parameters and retrieve some data that won't change with wavelength

# retrieve initial concentrations of lipids in experiments, and (for PC 22:6) change in concentration data

# PC 22:6
Init_conc_PC22_6_pmol_mL = 
  Exp_13_PC.samp.pmol.mL.norm.mean[c("PC 44:12"),c("Dark_control_2013-12-14 09:30:00.mean")]

d22_6_dt.UVA = d22_6_dt.Exp13.pmol_mL_hr.UVA.photo_ox[1]

d22_6_dt.UVB = d22_6_dt.Exp13.pmol_mL_hr.UVB.photo_ox[1]

d22_6_dt.UVA_UVB = d22_6_dt.Exp13.pmol_mL_hr.UVA_UVB.photo_ox[1]

# PC 22:6
Init_conc_PC22_6_pmol_mL = 
  Exp_13_PC.samp.pmol.mL.norm.mean[c("PC 44:12"),c("Dark_control_2013-12-14 09:30:00.mean")]


# effective pathlength through minor axis of quartz vials used in experiment
# since vials were left on their side throughout experiment = diameter (2 cm)
Quartz_vial_pathlength_cm = 2
m3_per_L = 0.001
cm_per_m = 100

##### (optional) code for calculation of uncertainities using a monte carlo analysis ##### 

# ------ optional monte carlo code
# define number of simulations
numSim = 5000

# preallocate matrix to hold values simulated during monte carlo run
Theta_Exp13_PC_22_6.sim = matrix(nrow = numSim, ncol = 3)

for (i in 1:numSim) {
  
  # user feedback
  cat("Simulation number:",as.character(i),"\n")
  
  # ------ end optional monte carlo code
  
# preallocate matrix to hold values of term under integral in Eq. 6 in manuscript
# first column will be for quartz vial, second for EPA vial
  
# term under integral will be in units of mol photons/volume/time (time is implied; it is the
# 8.2 hr sampling interval of the experiment)


Integral_Exp13_PC_22_6 = matrix(NA, ncol = 2, nrow =
                          length(PAL1314.JAZ.14Dec.E_n_p_sigma_umol_photons_m2))

colnames(Integral_Exp13_PC_22_6) = c("Integral_PC_22_6_quartz","Integral_PC_22_6_EPA")

for (j in 1:nrow(Integral_Exp13_PC_22_6)) { # iterate by wavelength
  
  if (lambda_nm_JAZ[j]<500) { # since don't have any lipid absorbance data above 500 nm
    
    # retrieve necessary values for this wavelength; some manipulations since
    # in some cases, variables were sampled at different intervals from each other
    
    # time-integrated photon flux at this wavelength (over sampling interval) 
    E_n_p_sigma_photons_m2 = PAL1314.JAZ.14Dec.E_n_p_sigma_umol_photons_m2[j]/1000000
    
    # vessel transmittances for this wavelength
    
    T_quartz = FracTrans$transmittance_quartz_pct[abs(FracTrans$lambda_nm-lambda_nm_JAZ[j])==min(abs(FracTrans$lambda_nm-lambda_nm_JAZ[j]))]
    
    T_EPA = FracTrans$transmittance_borosilicate_pct[abs(FracTrans$lambda_nm-lambda_nm_JAZ[j])==min(abs(FracTrans$lambda_nm-lambda_nm_JAZ[j]))]
    
    # total alpha of all elements in system; Eq. 8 in manuscript
    
    alpha_total_per_m = 
      
    # Nap molar absorption coefficient from LipidAbsData$epsilon_M_cm_PC22_6
    # calculated in LipidUV_VISAbsPlots.R
    epsilon_per_M_per_cm = LipidAbsData$epsilon_M_cm_PC22_6[LipidAbsData$lambda_nm==round(lambda_nm_JAZ[j])]

    # absorbance coefficient of Palmer seawater (the matrix) at this wavelength (and uncertainty)
    alpha_per_m = PAL1314_UV_VIS_SW_abs_profile_means$PAL1314_12Dec13_SW_alpha_per_m.mean[abs(PAL1314_UV_VIS_SW_abs_profile_means$Wavelength_nm-lambda_nm_JAZ[j])==min(abs(PAL1314_UV_VIS_SW_abs_profile_means$Wavelength_nm-lambda_nm_JAZ[j]))]
#     
    # ------ optional monte carlo code
    alpha_per_m.sd = PAL1314_UV_VIS_SW_abs_profile_means$PAL1314_12Dec13_SW_alpha_per_m.sd[abs(PAL1314_UV_VIS_SW_abs_profile_means$Wavelength_nm-lambda_nm_JAZ[j])==min(abs(PAL1314_UV_VIS_SW_abs_profile_means$Wavelength_nm-lambda_nm_JAZ[j]))]

    alpha_per_m = rnorm(1, mean = alpha_per_m, sd = alpha_per_m.sd)
    # ------ end optional monte carlo code
    
    # finally, calculate Ks for this wavelength
    # have to make some unit conversions
    
    # this is Ch 4, Equation 7
    
    Ks_Exp13_PC_22_6[j,1] =
      (E_mol_photons_m2*T_quartz*epsilon_per_M_per_cm*cm_per_m*m3_per_L*(1-10^(-alpha_per_m*(Quartz_vial_pathlength_cm/100))))/
      (alpha_per_m*(Quartz_vial_pathlength_cm/100))
    
    Ks_Exp13_PC_22_6[j,2] =
      (E_mol_photons_m2*T_EPA*epsilon_per_M_per_cm*cm_per_m*m3_per_L*(1-10^(-alpha_per_m*(Quartz_vial_pathlength_cm/100))))/
      (alpha_per_m*(Quartz_vial_pathlength_cm/100))
    
  }
  
}

##### calculation of polychromatic quantum yields, UVB, UVB, and total UVR (UVA+UVB) ##### 

# UVA and total UVR yields only use data out to 395.5 nm, since don't have good
# extinction data after that point

# retrieve initial concentration of PC 22:6 used in experiments, and change in concentration data
Init_conc_PC22_6_pmol_mL = 
  Exp_13_PC.samp.pmol.mL.norm.mean[c("PC 44:12"),c("Dark_control_2013-12-14 09:30:00.mean")]

d22_6_dt.UVA = d22_6_dt.Exp13.pmol_mL_hr.UVA.photo_ox[1]

d22_6_dt.UVB = d22_6_dt.Exp13.pmol_mL_hr.UVB.photo_ox[1]

d22_6_dt.UVA_UVB = d22_6_dt.Exp13.pmol_mL_hr.UVA_UVB.photo_ox[1]

# ------ optional monte carlo code

Init_conc_PC22_6_pmol_mL = rnorm(1, mean = Exp_13_PC.samp.pmol.mL.norm.mean[c("PC 44:12"),c("Dark_control_2013-12-14 09:30:00.mean")],
       sd = Exp_13_PC.samp.pmol.mL.norm.mean[c("PC 44:12"),c("Dark_control_2013-12-14 09:30:00.se")])

d22_6_dt.UVA = rnorm(1, mean = d22_6_dt.Exp13.pmol_mL_hr.UVA.photo_ox[1], sd = d22_6_dt.Exp13.pmol_mL_hr.UVA.photo_ox[2])

d22_6_dt.UVB = rnorm(1, mean = d22_6_dt.Exp13.pmol_mL_hr.UVB.photo_ox[1], sd = d22_6_dt.Exp13.pmol_mL_hr.UVB.photo_ox[2])

d22_6_dt.UVA_UVB = rnorm(1, mean = d22_6_dt.Exp13.pmol_mL_hr.UVA_UVB.photo_ox[1], sd = d22_6_dt.Exp13.pmol_mL_hr.UVA_UVB.photo_ox[2])
# ------ end optional monte carlo code

Theta_Exp13_PC_22_6 = vector(mode = "double", length = 3)

names(Theta_Exp13_PC_22_6) =
  c("Theta_PC_22_6_UVA_UVB","Theta_PC_22_6_UVA","Theta_PC_22_6_UVB")

# total UVR, use Ks for quartz vials

Ks_int_TUVR = caTools::trapz(lambda_nm_JAZ[(lambda_nm_JAZ>=290 & lambda_nm_JAZ<=395.5)],Ks_Exp13_PC_22_6[(lambda_nm_JAZ>=290 & lambda_nm_JAZ<=395.5),1])

Theta_Exp13_PC_22_6[1] =
  (-d22_6_dt.UVA_UVB*8.2)/
    (Init_conc_PC22_6_pmol_mL*
       Ks_int_TUVR)

# UVA, use Ks for EPA vials

Ks_int_UVA = caTools::trapz(lambda_nm_JAZ[(lambda_nm_JAZ>315 & lambda_nm_JAZ<=395.5)],Ks_Exp13_PC_22_6[(lambda_nm_JAZ>315 & lambda_nm_JAZ<=395.5),2])

Theta_Exp13_PC_22_6[2] =
  (-d22_6_dt.UVA*8.2)/
  (Init_conc_PC22_6_pmol_mL*
     Ks_int_UVA)

# UVB

Ks_int_UVB = caTools::trapz(lambda_nm_JAZ[(lambda_nm_JAZ>=290 & lambda_nm_JAZ<=315)],Ks_Exp13_PC_22_6[(lambda_nm_JAZ>=290 & lambda_nm_JAZ<=315),1])-
  caTools::trapz(lambda_nm_JAZ[(lambda_nm_JAZ>=290 & lambda_nm_JAZ<=315)],Ks_Exp13_PC_22_6[(lambda_nm_JAZ>=290 & lambda_nm_JAZ<=315),2])

Theta_Exp13_PC_22_6[3] =
  (-d22_6_dt.UVB*8.2)/
  (Init_conc_PC22_6_pmol_mL*
     Ks_int_UVB)

# ------ optional monte carlo code
Theta_Exp13_PC_22_6.sim[i,] = Theta_Exp13_PC_22_6

}

# return standard deviation of result
apply(Theta_Exp13_PC_22_6.sim,2,sd)

# store in a variable

Theta_Exp13_PC_22_6_sigma =
  apply(Theta_Exp13_PC_22_6.sim,2,sd)

# ------ end optional monte carlo code

# ##### d[22:6]/dt calculations, in pmol/L/hr, from Exp 03a #####
# # using the re-analyzed dataset from Nov 2016
# 
# # ***** the data underlying these calculations was not statisically significant; just calculating to obtain figures to report in the manuscript *****
# 
# # first, calculate d[22:6]/dt (change in concentration) using data
# # data for PC 22:6, 22:6 from Exp 03a; will need a few figures — one for apparent
# # dark oxidation (autooxidation component), one for UVB+UVA oxidation,
# # one for UVB oxidation
# 
# # define our dt, in hrs
# dt_Exp03a_hr = 9.8
# 
# # preallocate (second element will be uncertainty)
# d22_6_dt.Exp03a.pmol_mL_hr.dark = vector(mode = "double", length = 2)
# d22_6_dt.Exp03a.pmol_mL_hr.noUVB = vector(mode = "double", length = 2)
# d22_6_dt.Exp03a.pmol_mL_hr.plusUVB = vector(mode = "double", length = 2)
# d22_6_dt.Exp03a.pmol_mL_hr.UVB.photo_ox = vector(mode = "double", length = 2)
# d22_6_dt.Exp03a.pmol_mL_hr.UVA_UVB.photo_ox = vector(mode = "double", length = 2)
# d22_6_dt.Exp03a.pmol_mL_hr.UVA.photo_ox = vector(mode = "double", length = 2)
# 
# # calculations
# # dark
# d22_6_dt.Exp03a.pmol_mL_hr.dark[1] =
#   (Exp_03a_PC.samp.pmol.mL.norm.rerun.mean[c("PC 44:12"),c("Dark_control_no_HB_2013-12-02 18:28:00.mean")]-
#      Exp_03a_PC.samp.pmol.mL.norm.rerun.mean[c("PC 44:12"),c("Dark_control_no_HB_2013-12-02 08:40:00.mean")])/
#   dt_Exp03a_hr
# 
# d22_6_dt.Exp03a.pmol_mL_hr.dark[2] = (sqrt((Exp_03a_PC.samp.pmol.mL.norm.rerun.mean[c("PC 44:12"),c("Dark_control_no_HB_2013-12-02 18:28:00.se")])^2+
#                                              (Exp_03a_PC.samp.pmol.mL.norm.rerun.mean[c("PC 44:12"),c("Dark_control_no_HB_2013-12-02 08:40:00.se")])^2))/
#   dt_Exp03a_hr
# 
# # -UVB
# d22_6_dt.Exp03a.pmol_mL_hr.noUVB[1] =
#   (Exp_03a_PC.samp.pmol.mL.norm.rerun.mean[c("PC 44:12"),c("EPA_no_HB_2013-12-02 18:28:00.mean")]-
#      Exp_03a_PC.samp.pmol.mL.norm.rerun.mean[c("PC 44:12"),c("Dark_control_no_HB_2013-12-02 08:40:00.mean")])/
#   dt_Exp03a_hr
# 
# d22_6_dt.Exp03a.pmol_mL_hr.noUVB[2] = (sqrt((Exp_03a_PC.samp.pmol.mL.norm.rerun.mean[c("PC 44:12"),c("EPA_no_HB_2013-12-02 18:28:00.se")])^2+
#                                               (Exp_03a_PC.samp.pmol.mL.norm.rerun.mean[c("PC 44:12"),c("Dark_control_no_HB_2013-12-02 08:40:00.se")])^2))/
#   dt_Exp03a_hr
# 
# # +UVB
# d22_6_dt.Exp03a.pmol_mL_hr.plusUVB[1] =
#   (Exp_03a_PC.samp.pmol.mL.norm.rerun.mean[c("PC 44:12"),c("Quartz_no_HB_2013-12-02 18:28:00.mean")]-
#      Exp_03a_PC.samp.pmol.mL.norm.rerun.mean[c("PC 44:12"),c("Dark_control_no_HB_2013-12-02 08:40:00.mean")])/
#   dt_Exp03a_hr
# 
# d22_6_dt.Exp03a.pmol_mL_hr.plusUVB[2] = (sqrt((Exp_03a_PC.samp.pmol.mL.norm.rerun.mean[c("PC 44:12"),c("Quartz_no_HB_2013-12-02 18:28:00.se")])^2+
#                                                 (Exp_03a_PC.samp.pmol.mL.norm.rerun.mean[c("PC 44:12"),c("Dark_control_no_HB_2013-12-02 08:40:00.se")])^2))/
#   dt_Exp03a_hr
# 
# # # now, define control-corrected rates loss due to all photooxidation (direct+indirect),
# # # direct photooxidation, and indirect photo-oxidation
# # 
# # d22_6_dt.Exp03a.pmol_mL_hr.UVA_UVB.photo_ox[1] = d22_6_dt.Exp03a.pmol_mL_hr.plusUVB[1]-
# #   d22_6_dt.Exp03a.pmol_mL_hr.dark[1]
# # d22_6_dt.Exp03a.pmol_mL_hr.UVA_UVB.photo_ox[2] = sqrt((d22_6_dt.Exp03a.pmol_mL_hr.plusUVB[2])^2+
# #                                                         (d22_6_dt.Exp03a.pmol_mL_hr.dark[2])^2)
# # 
# # d22_6_dt.Exp03a.pmol_mL_hr.UVA.photo_ox[1] = d22_6_dt.Exp03a.pmol_mL_hr.noUVB[1]-
# #   d22_6_dt.Exp03a.pmol_mL_hr.dark[1]
# # d22_6_dt.Exp03a.pmol_mL_hr.UVA.photo_ox[2] = sqrt((d22_6_dt.Exp03a.pmol_mL_hr.noUVB[2])^2+
# #                                                     (d22_6_dt.Exp03a.pmol_mL_hr.dark[2])^2)
# # 
# # d22_6_dt.Exp03a.pmol_mL_hr.UVB.photo_ox[1] = d22_6_dt.Exp03a.pmol_mL_hr.UVA_UVB.photo_ox[1]-
# #   d22_6_dt.Exp03a.pmol_mL_hr.UVA.photo_ox[1]
# # d22_6_dt.Exp03a.pmol_mL_hr.UVB.photo_ox[2] = sqrt((d22_6_dt.Exp03a.pmol_mL_hr.UVA_UVB.photo_ox[2])^2+
# #                                                     (d22_6_dt.Exp03a.pmol_mL_hr.UVA.photo_ox[2])^2)
# 
# # given that rate of removal appears to have been about the same (i.e., all the way to zero concentration) in both the treatments and controls, assume alternatively that the rate represents a minimum rate of removal:
# 
# d22_6_dt.Exp03a.pmol_mL_hr.UVA_UVB.photo_ox[1] = d22_6_dt.Exp03a.pmol_mL_hr.plusUVB[1]
# 
# d22_6_dt.Exp03a.pmol_mL_hr.UVA_UVB.photo_ox[2] = sqrt((d22_6_dt.Exp03a.pmol_mL_hr.plusUVB[2])^2+
#                                                         (d22_6_dt.Exp03a.pmol_mL_hr.dark[2])^2)
# 
# d22_6_dt.Exp03a.pmol_mL_hr.UVA.photo_ox[1] = d22_6_dt.Exp03a.pmol_mL_hr.noUVB[1]
# 
# d22_6_dt.Exp03a.pmol_mL_hr.UVA.photo_ox[2] = sqrt((d22_6_dt.Exp03a.pmol_mL_hr.noUVB[2])^2+
#                                                     (d22_6_dt.Exp03a.pmol_mL_hr.dark[2])^2)
# 
# d22_6_dt.Exp03a.pmol_mL_hr.UVB.photo_ox[1] = d22_6_dt.Exp03a.pmol_mL_hr.UVA_UVB.photo_ox[1]
# 
# d22_6_dt.Exp03a.pmol_mL_hr.UVB.photo_ox[2] = sqrt((d22_6_dt.Exp03a.pmol_mL_hr.UVA_UVB.photo_ox[2])^2+
#                                                     (d22_6_dt.Exp03a.pmol_mL_hr.UVA.photo_ox[2])^2)
# 
# ##### calculation of total radiation received #####
# # i.e., the radiation received by the JAZ sensor, in the deck tank,
# # at the same depth as the vials, over course of Exp 03a
# 
# # integrated figures in units of photons/cm2/wavelength
# 
# # extract subset of data for the experimental time interval
# PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.02Dec.sub =
#   PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.02Dec[
#     PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.02Dec$Timestamp_GMT>=
#       as.POSIXct('2013-12-02 11:40:00', tz = "GMT") &  # 8:40 local time
#       PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.02Dec$Timestamp_GMT<=
#       as.POSIXct('2013-12-02 21:28:00', tz = "GMT") # 18:28 local time
#     ,]
# 
# # convert to units of photons/cm2, step by step
# PAL1314_JAZ_subsurf_hires_full_spectrum_W_m2.02Dec.sub =
#   PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.02Dec.sub[,5:(ncol(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.02Dec.sub)-1)]/100
# 
# lambda_nm_JAZ = as.numeric(colnames(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.02Dec.sub)[5:(ncol(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.02Dec.sub)-1)])
# 
# PAL1314_JAZ_subsurf_hires_full_spectrum_umol_photons_m2_s.02Dec.sub =
#   matrix(NA, ncol = length(lambda_nm_JAZ), nrow = nrow(PAL1314_JAZ_subsurf_hires_full_spectrum_uW_cm2.02Dec.sub))
# 
# for (i in 1:ncol(PAL1314_JAZ_subsurf_hires_full_spectrum_umol_photons_m2_s.02Dec.sub)) {
# 
#   PAL1314_JAZ_subsurf_hires_full_spectrum_umol_photons_m2_s.02Dec.sub[,i] =
#     PAL1314_JAZ_subsurf_hires_full_spectrum_W_m2.02Dec.sub[,i]*lambda_nm_JAZ[i]*0.836*(10^-2)
# 
# }
# 
# # integrate over time, by wavelength
# 
# # first, time integration at each wavelength
# 
# # need package caTools, if not loaded already
# 
# detach("package:RSEIS", unload=TRUE)
# library(caTools)
# 
# # preallocate vector
# PAL1314.JAZ.02Dec.sub_umol_photons_m2_total = vector(mode = "double",
#                                                      length = ncol(PAL1314_JAZ_subsurf_hires_full_spectrum_umol_photons_m2_s.02Dec.sub))
# 
# for (i in 1:length(PAL1314.JAZ.02Dec.sub_umol_photons_m2_total)) {
# 
#   # requires vector of times in second from UV_TS_analysis_PAL1314.R, with t = 0 being timepoint at beginning of experiment
# 
#   PAL1314.JAZ.02Dec.sub_umol_photons_m2_total[i] =
#     caTools::trapz(Exp03a_timeint.s[1:length(Exp03a_timeint.s)],PAL1314_JAZ_subsurf_hires_full_spectrum_umol_photons_m2_s.02Dec.sub[,i])
# 
# }
# 
# ##### calculation of wavelength-specific absorbances, from Exp 03a #####
# 
# # define a few parameters
# 
# # effective pathlength through minor axis of quartz vials used in experiment
# # since vials were left on their side throughout experiment
# Quartz_vial_pathlength_cm = 2
# m3_per_L = 0.001
# cm_per_m = 100
# 
# ##### (optional) code for calculation of uncertainities using a monte carlo analysis #####
# 
# # # ------ optional monte carlo code
# # # define number of simulations
# # numSim = 5000
# #
# # # preallocate matrix to hold values simulated during monte carlo run
# # Theta_Exp03a_PC_22_6.sim = matrix(nrow = numSim, ncol = 3)
# #
# # for (i in 1:numSim) {
# #
# #   # user feedback
# #   cat("Simulation number:",as.character(i),"\n")
# #
# #   # ------ end optional monte carlo code
# 
# # preallocate matrix to hold Ks values
# # first column will be for quartz vial, second for EPA vial
# 
# # Ks will be in units of mol photons/mol rxn/time (time is implied; it is the
# # 9.8 hr sampling interval of the experiment)
# 
# Ks_Exp03a_PC_22_6 = matrix(NA, ncol = 2, nrow =
#                              length(PAL1314.JAZ.02Dec.sub_umol_photons_m2_total))
# 
# colnames(Ks_Exp03a_PC_22_6) = c("Ks_PC_22_6_quartz","Ks_PC_22_6_EPA")
# 
# for (j in 1:nrow(Ks_Exp03a_PC_22_6)) { # iterate by wavelength
# 
#   if (lambda_nm_JAZ[j]<800) { # since don't have any lipid absorbance data above 800 nm
# 
#     # retrieve necessary values for this wavelength; some manipulations since
#     # in some cases, variables were sampled at different intervals from each other
# 
#     # time-integrated photon flux at this wavelength (over sampling interval)
#     E_mol_photons_m2 = PAL1314.JAZ.02Dec.sub_umol_photons_m2_total[j]/1000000
# 
#     # decadic molar absorption coefficient from LipidAbsData$epsilon_M_cm_PC22_6_dil
#     # calculated in LipidUV_VISAbsPlots.R
#     epsilon_per_M_per_cm = LipidAbsData$epsilon_M_cm_PC22_6_dil[LipidAbsData$lambda_nm==round(lambda_nm_JAZ[j])]
# 
#     # absorbance coefficient of Palmer seawater (the matrix) at this wavelength (and uncertainty)
#     alpha_per_m = PAL1314_UV_VIS_SW_abs_profile_means$PAL1314_12Dec13_SW_alpha_per_m.mean[abs(PAL1314_UV_VIS_SW_abs_profile_means$Wavelength_nm-lambda_nm_JAZ[j])==min(abs(PAL1314_UV_VIS_SW_abs_profile_means$Wavelength_nm-lambda_nm_JAZ[j]))]
#     #
#     #       # ------ optional monte carlo code
#     #       alpha_per_m.sd = PAL1314_UV_VIS_SW_abs_profile_means$PAL1314_12Dec13_SW_alpha_per_m.sd[abs(PAL1314_UV_VIS_SW_abs_profile_means$Wavelength_nm-lambda_nm_JAZ[j])==min(abs(PAL1314_UV_VIS_SW_abs_profile_means$Wavelength_nm-lambda_nm_JAZ[j]))]
#     #
#     #       alpha_per_m = rnorm(1, mean = alpha_per_m, sd = alpha_per_m.sd)
#     #       # ------ end optional monte carlo code
# 
#     # vessel transmittances for this wavelength
# 
#     T_quartz = FracTrans$transmittance_quartz_pct[abs(FracTrans$lambda_nm-lambda_nm_JAZ[j])==min(abs(FracTrans$lambda_nm-lambda_nm_JAZ[j]))]
# 
#     T_EPA = FracTrans$transmittance_borosilicate_pct[abs(FracTrans$lambda_nm-lambda_nm_JAZ[j])==min(abs(FracTrans$lambda_nm-lambda_nm_JAZ[j]))]
# 
#     # finally, calculate Ks for this wavelength
#     # have to make some unit conversions
# 
#     Ks_Exp03a_PC_22_6[j,1] =
#       (E_mol_photons_m2*T_quartz*epsilon_per_M_per_cm*cm_per_m*m3_per_L*(1-10^(-alpha_per_m*(Quartz_vial_pathlength_cm/100))))/
#       (alpha_per_m*(Quartz_vial_pathlength_cm/100))
# 
#     Ks_Exp03a_PC_22_6[j,2] =
#       (E_mol_photons_m2*T_EPA*epsilon_per_M_per_cm*cm_per_m*m3_per_L*(1-10^(-alpha_per_m*(Quartz_vial_pathlength_cm/100))))/
#       (alpha_per_m*(Quartz_vial_pathlength_cm/100))
# 
#   }
# 
# }
# 
# ##### calculation of polychromatic quantum yields, UVB, UVB, and total UVR (UVA+UVB) #####
# 
# # UVA and total UVR yields only use data out to 395.5 nm, since don't have good
# # extinction data after that point
# 
# # retrieve initial concentration of PC 22:6 used in experiments, and change in concentration data
# Init_conc_PC22_6_pmol_mL =
#   Exp_03a_PC.samp.pmol.mL.norm.rerun.mean[c("PC 44:12"),c("Dark_control_no_HB_2013-12-02 08:40:00.mean")]
# 
# d22_6_dt.UVA = d22_6_dt.Exp03a.pmol_mL_hr.UVA.photo_ox[1]
# 
# d22_6_dt.UVB = d22_6_dt.Exp03a.pmol_mL_hr.UVB.photo_ox[1]
# 
# d22_6_dt.UVA_UVB = d22_6_dt.Exp03a.pmol_mL_hr.UVA_UVB.photo_ox[1]
# 
# #   # ------ optional monte carlo code
# #
# #   Init_conc_PC22_6_pmol_mL = rnorm(1, mean = Exp_03a_PC.samp.pmol.mL.norm.rerun.mean[c("PC 44:12"),c("Dark_control_no_HB_2013-12-02 08:40:00.mean")],
# #                                    sd = Exp_03a_PC.samp.pmol.mL.norm.rerun.mean[c("PC 44:12"),c("Dark_control_no_HB_2013-12-02 08:40:00.se")])
# #
# #   d22_6_dt.UVA = rnorm(1, mean = d22_6_dt.Exp03a.pmol_mL_hr.UVA.photo_ox[1], sd = d22_6_dt.Exp03a.pmol_mL_hr.UVA.photo_ox[2])
# #
# #   d22_6_dt.UVB = rnorm(1, mean = d22_6_dt.Exp03a.pmol_mL_hr.UVB.photo_ox[1], sd = d22_6_dt.Exp03a.pmol_mL_hr.UVB.photo_ox[2])
# #
# #   d22_6_dt.UVA_UVB = rnorm(1, mean = d22_6_dt.Exp03a.pmol_mL_hr.UVA_UVB.photo_ox[1], sd = d22_6_dt.Exp03a.pmol_mL_hr.UVA_UVB.photo_ox[2])
# #   # ------ end optional monte carlo code
# 
# Theta_Exp03a_PC_22_6 = vector(mode = "double", length = 3)
# 
# names(Theta_Exp03a_PC_22_6) =
#   c("Theta_PC_22_6_UVA_UVB","Theta_PC_22_6_UVA","Theta_PC_22_6_UVB")
# 
# # total UVR, use Ks for quartz vials
# 
# Ks_int_TUVR = caTools::trapz(lambda_nm_JAZ[(lambda_nm_JAZ>=290 & lambda_nm_JAZ<=395.5)],Ks_Exp03a_PC_22_6[(lambda_nm_JAZ>=290 & lambda_nm_JAZ<=395.5),1])
# 
# Theta_Exp03a_PC_22_6[1] =
#   (-d22_6_dt.UVA_UVB*8.2)/
#   (Init_conc_PC22_6_pmol_mL*
#      Ks_int_TUVR)
# 
# # UVA, use Ks for EPA vials
# 
# Ks_int_UVA = caTools::trapz(lambda_nm_JAZ[(lambda_nm_JAZ>315 & lambda_nm_JAZ<=395.5)],Ks_Exp03a_PC_22_6[(lambda_nm_JAZ>315 & lambda_nm_JAZ<=395.5),2])
# 
# Theta_Exp03a_PC_22_6[2] =
#   (-d22_6_dt.UVA*8.2)/
#   (Init_conc_PC22_6_pmol_mL*
#      Ks_int_UVA)
# 
# # UVB
# 
# Ks_int_UVB = caTools::trapz(lambda_nm_JAZ[(lambda_nm_JAZ>=290 & lambda_nm_JAZ<=315)],Ks_Exp03a_PC_22_6[(lambda_nm_JAZ>=290 & lambda_nm_JAZ<=315),1])-
#   caTools::trapz(lambda_nm_JAZ[(lambda_nm_JAZ>=290 & lambda_nm_JAZ<=315)],Ks_Exp03a_PC_22_6[(lambda_nm_JAZ>=290 & lambda_nm_JAZ<=315),2])
# 
# Theta_Exp03a_PC_22_6[3] =
#   (-d22_6_dt.UVB*8.2)/
#   (Init_conc_PC22_6_pmol_mL*
#      Ks_int_UVB)
# 
# # ------ optional monte carlo code
# Theta_Exp03a_PC_22_6.sim[i,] = Theta_Exp03a_PC_22_6
# 
# }
# 
# # return standard deviation of result
# apply(Theta_Exp03a_PC_22_6.sim,2,sd)
# 
# # ------ end optional monte carlo code
