# PAL1314_AQY_calc.R

# Purpose: Calculate apparent quantum yields for photolysis of PUFA-containing
# PC species examined in the PAL1314 liposome experiments 

# ****** Assumes some variables created by PAL1314_liposome_expts.R and
# UV_TS_analysis_PAL1314.R are already in user's workspace; these files can be
# found at https://github.com/jamesrco/LipidPhotoOxBox

# Created 10/26/16 by J.R.C.

# set working directory to parent LipidPhotoOxBox repo
setwd("/Users/jrcollins/Code/LipidPhotoOxBox")

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

# 27 Dec 13 calcs

PAL1314.SW.abs.subset.27Dec13 = PAL1314_LMG1401_UV_VIS_SW_abs[,colnames(PAL1314_LMG1401_UV_VIS_SW_abs) %in% PAL1314_DOC_log_metadata$Sample_ID[PAL1314_DOC_log_metadata$Date %in% c("27-Dec-13")]]

PAL1314_UV_VIS_SW_abs_profile_means$PAL1314_27Dec13_SW_alpha_per_m.mean = apply(
  PAL1314.SW.abs.subset.27Dec13,1,mean)*log(10)/(pathlength_cm.PAL1314_SW_abs/100)

PAL1314_UV_VIS_SW_abs_profile_means$PAL1314_27Dec13_SW_alpha_per_m.sd = apply(
  PAL1314.SW.abs.subset.27Dec13,1,sd)*log(10)/(pathlength_cm.PAL1314_SW_abs/100)

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

##### AQY calculation, after Kieber at al, Env. Sci P&I, 2014 ##### 

# first, calculate d[22:6]/dt (change in concentration) using data
# data for PC 22:6, 22:6 from Exp 13; will need a few figures — one for apparent
# dark oxidation (autooxidation component), one for -UVB (indirect photolysis),
# one for +UVB

d22_6_dt.pmol_L_hr.autoox
d22_6_dt.pmol_L_hr.autoox
d22_6_dt.pmol_L_hr.autoox

Exp_13_PC.samp.pmol.mL.norm.mean$

