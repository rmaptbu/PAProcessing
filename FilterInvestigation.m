clear;
corrmean=2001;
for corrsize=500:500:4000
    corrmin=corrmean-corrsize/2;
    corrmax=corrmean+corrsize/2;
    lowpass=250;
    highpass=10;
    wallfilter=0;
    figname=['WallOffNocorrsize',num2str(corrsize)];
    ReadFiles(highpass,lowpass,wallfilter,corrmin,corrmax,figname);
end

corrsize=2000;
corrmin=corrmean-corrsize/2;
corrmax=corrmean+corrsize/2;
for highpass=2:2:14
    lowpass=250;
    wallfilter=1;
    figname=['WallOnHighpass',num2str(highpass)];
    ReadFiles(highpass,lowpass,wallfilter,corrmin,corrmax,figname);
end

for highpass=2:2:14
    lowpass=0;
    wallfilter=1;
    figname=['WallOnNoLowpassHighpass',num2str(highpass)];
    ReadFiles(highpass,lowpass,wallfilter,corrmin,corrmax,figname);
end