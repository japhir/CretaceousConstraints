#### BEGIN ANALYSIS
# Load the R package Astrochron
# This analysis uses version 0.6.6. Please use versions >= 0.6.6.

library(astrochron);

##################################################################
# ANALYZE THE LIBSACK FMI DATA
##################################################################
# (1) Read the Libsack FMI data from Locklair & Sageman (2008).
#     This should be a comma-separated-value file (.csv), with first column as depth (meters)
# and second column as FMI.

FMI=read("data.csv");

# (2) This data set has a sampling interval that ranges from 0.030478 to 0.030479 m.
# Interpolate the data to a sampling interval of 0.03 m.

FMI_0.03=linterp(FMI, dt=0.03);

# (3) Tune the FMI record using the long-eccentricity cycle.
# (a) Extract Locklair & Sageman’s (2008) long-eccentricity cycle using bandpass filtering.

longEcc=bandpass(FMI_0.03,flow=0.14,fhigh=0.26,xmax=0.5,padfac=5,win=2,p=0.66);

# (b) Find the peak associated with each long eccentricity cycle.

eccMax=peak(longEcc);

# (c) Construct the floating time (elapsed time) vs. core depth map for tuning.

timeControl=cb(eccMax[,2],(0:18)*405);

plot(timeControl,type="l",lwd=2,xlab="Depth (m)",ylab="Elapsed Time (ka)");

# (d) Tune (the original FMI data) using the time vs. core depth map.

tuned=tune(FMI,timeControl,extrapolate=T);

# (4) This tuned data set has a sampling interval that ranges from 1.714387 to 3.774922 ka.
# Interpolate the tuned data to the median sampling interval of ~2.5 ka.

tuned_2.5=linterp(tuned, dt=2.5);

# (5) Convert result from floating (elapsed) time to radioisotopically-anchored time,
# using the nominal radioisotopic anchoring (S.p. ammonite biozone).
# The radioisotopic age, and its depth in the Libsack core, come from Table 1.

anchorAt=resample(timeControl,xout=2147.62,genplot=F)[,2];
anchored=anchorTime(tuned_2.5,time=anchorAt,age=89370,timeDir=2);

# Now create a plot of radioisotopically-anchored time vs. depth.
timeDepth=tuned;
timeDepth[2]=FMI[1];
anchoredTimeDepth=anchorTime(timeDepth,time=anchorAt,age=89370,timeDir=2,genplot=F);

pl(1);
plot(anchoredTimeDepth,type="l",lwd=2,col="red",ylim=c(max(FMI[1]),min(FMI[1])),xlab="Time (ka)",ylab="Depth (m)",cex.lab=1.2);

# (6) Conduct evolutive power spectral analysis (EPSA) and evolutive harmonic analysis (EHA)
# for the tuned & anchored FMI data using a 500-ka moving window (with linear trend
# removal), and three 2pi prolate tapers. Plot amplitude normalized to unity
# (for each window) to reveal changes in relative strength.

pwr=eha(anchored,win=500,fmax=.1,output=2,pl=1,pad=5000,genplot=3,ydir=-1, xlab="Frequency (cycles/ka)",ylab="Age (ka)");

# (7) Determine power modulation of obliquity terms using EPSA results.
# Integrate the obliquity power from 0.018 to 0.037 cycles/ka.

integrate_obl=integratePower(pwr,flow=0.018,fhigh=0.037,npts=201,pad=5000,ln=T,ydir=-1);

# (8) Determine power modulation of short eccentricity terms.
# Integrate the short eccentricity power from 0.007 to 0.012 cycles/ka.

integrate_ecc=integratePower(pwr,flow=0.007,fhigh=0.012,npts=201,pad=5000,ln=T,ydir=-1);

# (9) Evaluate amplitude modulation of the long eccentricity term (405 ka), following
# removal of bias associated with long-term (>1 Myr) variance, using Lowess.

longEcc2=bandpass(noLow(anchored,0.1,genplot=F),flow=.002,fhigh=.0035,win=2,p=0.66, padfac=5,xmax=.02);

hilEcc2=hilbert(longEcc2,addmean=T);

# (10) Plot summary figures.

xlim1=c(82707.41,89932.41);
pl(r=3,c=1);
plot(cb(integrate_obl,c(1,2)),type="l",lwd=2,col="red",ylab="Obliquity Band Power",xlab="Time (ka)", cex.lab=1.2,xlim=xlim1);
plot(cb(integrate_ecc,c(1,2)),type="l",lwd=2,col="red",ylab="Short-Eccentricity Band Power",xlab="Time (ka)",cex.lab=1.2,xlim=xlim1);
plot(longEcc2,type="l",lwd=2,col="red",ylab="Long-Eccentricity Bandpass",xlab="Time (ka)", cex.lab=1.2,xlim=xlim1);
lines(hilEcc2);
pl(r=3,c=1);
plot(cb(integrate_obl,c(1,4)),type="l",lwd=2,col="red",ylab="Obliquity/Total Power",xlab="Time (ka)", cex.lab=1.2,xlim=xlim1);
plot(cb(integrate_ecc,c(1,4)),type="l",lwd=2,col="red",ylab="Short-Eccentricity/Total Power",xlab="Time (ka)",cex.lab=1.2,xlim=xlim1);
plot(longEcc2,type="l",lwd=2,col="red",ylab="Long-Eccentricity Bandpass",xlab="Time (ka)",
cex.lab=1.2,xlim=xlim1);
lines(hilEcc2);
pl(r=3,c=1);
plot(cb(integrate_ecc,c(1,2)),type="l",lwd=2,col="red",ylab="Short-Eccentricity Band Power",xlab="Time (ka)",cex.lab=1.2,xlim=xlim1);
plot(cb(integrate_ecc,c(1,4)),type="l",lwd=2,col="red",ylab="Short-Eccentricity/Total Power",xlab="Time (ka)",cex.lab=1.2,xlim=xlim1);
plot(longEcc2,type="l",lwd=2,col="red",ylab="Long-Eccentricity Bandpass",xlab="Time    (ka)",cex.lab=1.2,xlim=xlim1);
lines(hilEcc2);
