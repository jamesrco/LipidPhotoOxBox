# PAL1718_UV_data_analysis.R

# for PAL1718 data

# We have actinometry this year. So, checking on some older PAL UV data to see
# what the daily fluxes look like in the active wavebands for the actinometers\\

# Ref for actinometry is Jankowski, Kieber, Mopper, Photochemistry and
# Photobiology, 1999, 70(3): 319–328

# Calculate some time-integrated doses from the NOAA data
# *** requires the object PAL1314_est_spectra_at_depth_uW_cm2 from
# UV_TS_analysis_PAL1314.R ***

# preallocate an object to hold the data, will be equiv to E_n,p,sigma in Eq. 10
# in the now-rejected GCA manuscript

PAL1314_E_n_p_sigma_umol_photons_m2_d_nm_incident =
  as.data.frame(matrix(NA,length(lipidox.calcdates),length(NOAA_AntUV_lambdas)))
rownames(PAL1314_E_n_p_sigma_umol_photons_m2_d_nm_incident) = lipidox.calcdates
colnames(PAL1314_E_n_p_sigma_umol_photons_m2_d_nm_incident) = NOAA_AntUV_lambdas

for (i in 1:length(lipidox.calcdates)) {
  
  # iterate through dates on which we want an estimate, from lipidox.calcdates
  
  # integrated figures in units of photons/cm2/wavelength
  
  # extract subset of data for this date
  
  PAL1314_est_spectra.today.uW_cm2 =  
    PAL1314_est_spectra_at_depth_uW_cm2[[1]][
      PAL1314_est_spectra_at_depth_uW_cm2$Timestamp_GMT>
        lipidox.calcdates[i]-60*60*12 & 
        PAL1314_est_spectra_at_depth_uW_cm2$Timestamp_GMT<=
        lipidox.calcdates[i]+60*60*12,]
  
  Todays.timestamps = 
    PAL1314_est_spectra_at_depth_uW_cm2$Timestamp_GMT[
      PAL1314_est_spectra_at_depth_uW_cm2$Timestamp_GMT>
        lipidox.calcdates[i]-60*60*12 & 
        PAL1314_est_spectra_at_depth_uW_cm2$Timestamp_GMT<=
        lipidox.calcdates[i]+60*60*12]
  
  # create a vector of times in seconds
  
  Timevector.s = Todays.timestamps-(lipidox.calcdates[i]-60*60*12)
  
  units(Timevector.s) = "secs"
  
  Timevector.s = as.numeric(Timevector.s)
  
  # convert to units of photons/cm2, step by step
  PAL1314_est_spectra.today.W_m2 = 
    PAL1314_est_spectra.today.uW_cm2/100
  
  PAL1314_est_spectra.today.umol_photons_m2_s =
    matrix(NA, ncol = length(NOAA_AntUV_lambdas), nrow = nrow(PAL1314_est_spectra.today.uW_cm2))
  
  for (k in 1:ncol(PAL1314_est_spectra.today.umol_photons_m2_s)) {
    
    PAL1314_est_spectra.today.umol_photons_m2_s[,k] =
      PAL1314_est_spectra.today.W_m2[,k]*NOAA_AntUV_lambdas[k]*0.836*(10^-2)
    
  }
  
  # integrate over time, by wavelength
  
  # first, time integration at each wavelength
  
  # preallocate vector
  PAL1314_E_n_p_sigma_umol_photons_m2_d_nm.today = vector(mode = "double",
                                                          length = ncol(PAL1314_est_spectra.today.umol_photons_m2_s))
  
  for (l in 1:length(PAL1314_E_n_p_sigma_umol_photons_m2_d_nm.today)) {
    
    # requires vector of times in seconds, with t = 0 being beginning of day
    
    PAL1314_E_n_p_sigma_umol_photons_m2_d_nm.today[l] =
      caTools::trapz(Timevector.s[1:length(Timevector.s)],PAL1314_est_spectra.today.umol_photons_m2_s[,l])
    
  }
  
  # store data in appropriate location
  
  PAL1314_E_n_p_sigma_umol_photons_m2_d_nm_incident[i,] =
    PAL1314_E_n_p_sigma_umol_photons_m2_d_nm.today
  
}

# now, can calculate some integrals for certain wavebands

# nitrite actinometer response bw is 330-380 nm
# nitrate is 311-333

PAL1314_E_n_p_sigma_umol_photons_m2_d_330_380nm =
  as.data.frame(matrix(NA,nrow(PAL1314_E_n_p_sigma_umol_photons_m2_d_nm_incident),1))
rownames(PAL1314_E_n_p_sigma_umol_photons_m2_d_330_380nm) =
  rownames(PAL1314_E_n_p_sigma_umol_photons_m2_d_nm_incident)

for (i in 1:nrow(PAL1314_E_n_p_sigma_umol_photons_m2_d_330_380nm)) {
  PAL1314_E_n_p_sigma_umol_photons_m2_d_330_380nm[i,1] =
    caTools::trapz(NOAA_AntUV_lambdas[(NOAA_AntUV_lambdas>=330 & NOAA_AntUV_lambdas<=380)],
                   as.numeric(PAL1314_E_n_p_sigma_umol_photons_m2_d_nm_incident[i,NOAA_AntUV_lambdas>=330 & 
                                                                                  NOAA_AntUV_lambdas<=380]))
}

PAL1314_E_n_p_sigma_umol_photons_m2_d_311_333nm =
  as.data.frame(matrix(NA,nrow(PAL1314_E_n_p_sigma_umol_photons_m2_d_nm_incident),1))
rownames(PAL1314_E_n_p_sigma_umol_photons_m2_d_311_333nm) =
  rownames(PAL1314_E_n_p_sigma_umol_photons_m2_d_nm_incident)

for (i in 1:nrow(PAL1314_E_n_p_sigma_umol_photons_m2_d_311_333nm)) {
  PAL1314_E_n_p_sigma_umol_photons_m2_d_311_333nm[i,1] =
    caTools::trapz(NOAA_AntUV_lambdas[(NOAA_AntUV_lambdas>=311 & NOAA_AntUV_lambdas<=333)],
                   as.numeric(PAL1314_E_n_p_sigma_umol_photons_m2_d_nm_incident[i,NOAA_AntUV_lambdas>=311 & 
                                                                                  NOAA_AntUV_lambdas<=333]))
}

