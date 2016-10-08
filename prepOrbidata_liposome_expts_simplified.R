
# ******************************************************************
################ Basic user begin editing here #############
# ******************************************************************

################ User: define locations of data files and database(s) #############

working_dir = "/Volumes/Lab/Jamie Collins/mzXML/PAL1314/Exp_13/" # specify working directory
setwd(working_dir) # set working directory to working_dir

# specify directories subordinate to the working directory in which the .mzXML files for xcms can be found; per xcms documentation, use subdirectories within these to divide files according to treatment/primary environmental variable (e.g., station number along a cruise transect) and file names to indicate timepoint/secondary environmental variable (e.g., depth)

mzXMLdirs = c("neg/","pos/")

# specify which of the directories above you wish to analyze this time through 

chosenFileSubset = "neg/"

################# Load in mzXML files, get xcms settings from IPO or user input #############

# load selected subset for processing

mzXMLfiles.raw = list.files(chosenFileSubset, recursive = TRUE, full.names = TRUE)

# verify the ion mode of the data in these files

# subset.polarity = getSubsetIonMode(mzXMLfiles.raw)

# provide some feedback to user

print(paste0("Loaded ",length(mzXMLfiles.raw)," mzXML files. These files contain ",subset.polarity," ion mode data. Raw dataset consists of:"))

print(mzXMLfiles.raw)

# check whether user has elected to exclude any files, and exclude them if they happen to be in this subset

if (exists("excluded.mzXMLfiles") & length("excluded.mzXMLfiles")>0) {
  
  excludedfiles = getFNmatches(IDnumlist = excluded.mzXMLfiles, filelist = mzXMLfiles.raw) # index files to be excluded
  
  print(paste0("The following files will be excluded from processing based on user's input:"))
  print(mzXMLfiles.raw[excludedfiles])
        
  mzXMLfiles = mzXMLfiles.raw[-excludedfiles] # exclude the files from mzXMLfiles
  
} else {
  
  mzXMLfiles = mzXMLfiles.raw
  
}

################# Create xcmsSet using selected settings #############

print(paste0("Creating xcmsSet object from ",length(mzXMLfiles)," mzXML files remaining in dataset using specified settings..."))

# create xcms xset object; runs WAY faster with multicore tasking enabled; 

xset_centWave = xcmsSet(mzXMLfiles,
                        method = "centWave",
                        profparam = centW.profparam, 
                        ppm = centW.ppm,
                        peakwidth = c(centW.min_peakwidth,centW.max_peakwidth),
                        fitgauss = centW.fitgauss,
                        noise = centW.noise,
                        mzdiff = centW.mzdiff,
                        verbose.columns = centW.verbose.columns,
                        snthresh = centW.snthresh,
                        integrate = centW.integrate,
                        prefilter = centW.prefilter,
                        mzCenterFun = centW.mzCenterFun,
                        #                 sleep = centW.sleep
                        nSlaves = centW.nSlaves
)

print(paste0("xcmsSet object xset_centWave created:"))

print(xset_centWave)

# Some notes:
#
#  1. If using massifquant or centWave and you are sure your input data are centroided, can ignore warning message "It looks like this file is in profile mode. [method] can process only centroid mode data !" since this is just based on a heuristic. That is, you can ignore the message if you are certain data are in centroid mode. You can verify this by opening one of your converted .mzXML files in a text reader. You should see: <dataProcessing centroided="1"></dataProcessing> (a "0" is bad)
# 
#     For more on this error, see http://metabolomics-forum.com/viewtopic.php?f=8&t=267 or https://groups.google.com/forum/#!topic/xcms/xybDDQTaQiY
#
#  2. So long as the number of peak data insertion problems is relatively low (i.e., < 100), you can safely ignore the error. Otherwise, might try lowering the ppm
#
#  3. On-the-fly plotting features (i.e., with sleep ≥ 0.001 enabled) don't appear to function properly in Mac RStudio

#####################################################################################
##### Grouping and retention time correction using xcms (and IPO, if desired) #######
#####################################################################################

################# Perform grouping and retention time correction on dataset #############

print(paste0("Performing grouping and retention time correction on dataset"))
print(paste0("Using group.density and retcor.",retcor.meth))

# initial grouping

# # method "nearest" with settings below seems to work better than method = "density," but takes absolutely forever; however, it seems to take less time crunching centWave picked data than massifquant picked data

# xset_centWave = group(xset_centWave,
# method = "nearest",
# mzVsRTbalance=10,
# mzCheck=0.2,
# rtCheck=30,
# kNN=10
# ) 

# using method = "density" with settings from above

xset_gr = group(xset_centWave,
                method = "density",
                bw = density.bw,
                minfrac = density.minfrac,
                minsamp = density.minsamp,
                mzwid = density.mzwid,
                max = density.max,
                sleep = density.sleep
)

# chromatographic alignment (retention time correction)

if (retcor.meth=="loess") {
  
  xset_gr.ret = retcor(xset_gr,
#                        method = "loess", # this appears unnecessary
                        missing = loess.missing,
                        extra = loess.extra,
                        smooth = "loess",
                        span = loess.span,
                        family = loess.family,
                        plottype = loess.plottype,
                        col = NULL,
                        ty = NULL
  )
    
} else if (retcor.meth=="obiwarp") {
  
  xset_gr.ret = retcor.peakgroups(xset_gr,
                        method = "obiwarp",
                        plottype = obiwarp.plottype,
                        profStep = obiwarp.profStep,
                        center = obiwarp.center,
                        response = obiwarp.response,
                        distFunc = obiwarp.distFunc,
                        gapInit = obiwarp.gapInit,
                        gapExtend = obiwarp.gapInit,
                        factorDiag = obiwarp.factorDiag,
                        factorGap = obiwarp.factorGap,
                        localAlignment = obiwarp.localAlignment,
                        initPenalty = 0
  )
  
}

# perform grouping again

print(paste0("Performing second peak grouping after application of retcor..."))

# using method = "density" with settings from above

xset_gr.ret.rg = group(xset_gr.ret,
                      method = "density",
                      bw = density.bw,
                      minfrac = density.minfrac,
                      minsamp = density.minsamp,
                      mzwid = density.mzwid,
                      max = density.max,
                      sleep = density.sleep
)

# fill missing peaks

print(paste0("Filling missing peaks..."))

xset_gr.ret.rg.fill = fillPeaks.chrom(xset_gr.ret.rg, nSlaves = 4)

#####################################################################################
##### Isotope peak identification, creation of xsAnnotate object using CAMERA #######
#####################################################################################

print(paste0("Applying CAMERA to identify isotopic peaks, create xsAnnotate object, and create CAMERA pseudospectra using correlation of xcms peak groups between and within samples. These pseudospectra are the groups within which the adduct hierarchy and retention time screening criteria will be applied using LOBSTAHS"))

# first, a necessary workaround to avoid a import error; see https://support.bioconductor.org/p/69414/
imports = parent.env(getNamespace("CAMERA"))
unlockBinding("groups", imports)
imports[["groups"]] = xcms::groups
lockBinding("groups", imports)

# create annotated xset using wrapper annotate(), allowing us to perform all CAMERA tasks at once

xset_a = annotate(xset_gr.ret.rg.fill,
                  
                  quick=FALSE, # set to FALSE because we want to run groupCorr; will also cause CAMERA to run adduct annotation. while LOBSTAHS will do its own adduct identification later, it doesn't hurt to do this now if it lets CAMERA create better pseudospectra  
                  sample=NA, # use all samples
                  nSlaves=4, # use 4 sockets
                  
                  # group FWHM settings
                  # using defaults for now
                  
                  sigma=6,
                  perfwhm=0.6,
                  
                  # groupCorr settings
                  # using defaults for now
                  
                  cor_eic_th=0.75,
                  graphMethod="hcs",
                  pval=0.05,
                  calcCiS=TRUE,
                  calcIso=TRUE,
                  calcCaS=FALSE, # weird results with this set to TRUE
                  
                  # findIsotopes settings
                  
                  maxcharge=4,
                  maxiso=4,
                  minfrac=0.5, # 0.25?
                  
                  # adduct annotation settings
                  
                  psg_list=NULL,
                  rules=NULL,
                  polarity=subset.polarity,
                  multiplier=3,
                  max_peaks=100,
                  
                  # common to multiple tasks
                  
                  intval="into",
                  ppm=2.5,
                  mzabs=0.0015
                  
                  )

cleanParallel(xset_a) # kill sockets

# at this point, should have an xsAnnotate object called "xset_a" in hand, which will serve as the primary input to the main screening and annotation function "doLOBscreen" in LOBSTAHS

print(paste0("xsAnnotate object 'xset_a' has been created. User can now use LOBSTAHS to perform screening..."))

print(xset_a)

library(LOBSTAHS)
LOBset = doLOBscreen(xset_a, polarity="negative",match.ppm=2.5)
getLOBpeaklist(LOBset,gen.csv=TRUE)

Exp_13_neg = LOBset
save(file = "Exp_13_neg.RData", Exp_13_neg)
