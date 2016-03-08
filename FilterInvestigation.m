clear;
corrmean=2001;
% for corrsize=500:500:500
%     corrmin=corrmean-corrsize/2;
%     corrmax=corrmean+corrsize/2;
%     lowpass=250;
%     highpass=10;
%     wallfilter=0;
%     figname=['test',num2str(corrsize)];
%     ReadFiles(highpass,lowpass,wallfilter,corrmin,corrmax,figname);
% end

corrmin=1;
corrmax=5000;
highpass=10;
lowpass=0;
wallfilter=0;
figname=['WallOffLowpass0Highpass',num2str(highpass)];
ReadFiles(highpass,lowpass,wallfilter,corrmin,corrmax,figname);

corrsize=500;
corrmin=1;
corrmax=5000;
for highpass=2:2:14
    lowpass=0;
    wallfilter=0;
    figname=['WallOffLowpass0Highpass',num2str(highpass)];
    ReadFiles(highpass,lowpass,wallfilter,corrmin,corrmax,figname);
end

for highpass=2:2:14
    lowpass=250;
    wallfilter=1;
    figname=['WallOnLowpass250Highpass',num2str(highpass)];
    ReadFiles(highpass,lowpass,wallfilter,corrmin,corrmax,figname);
end

for highpass=2:2:14
    lowpass=250;
    wallfilter=0;
    figname=['WallOffLowpass250Highpass',num2str(highpass)];
    ReadFiles(highpass,lowpass,wallfilter,corrmin,corrmax,figname);
end