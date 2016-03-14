clearvars -except emdd;

corrmin=1;
corrmax=5000;
highpass=00;
lowpass=250;
wallfilter=0;
emd_high=0;
emd_low=0;
for i=1:2
    disp(i);
    order=i;
    figname=['IMFXcorrs_Mode',num2str(order)];
    ReadFiles(highpass,lowpass,wallfilter,corrmin,corrmax,figname, emdd, emd_low, emd_high,order);
end
for i=5:6
    disp(i);
    order=i;
    figname=['IMFXcorrs_Mode',num2str(order)];
    ReadFiles(highpass,lowpass,wallfilter,corrmin,corrmax,figname, emdd, emd_low, emd_high,order);
end