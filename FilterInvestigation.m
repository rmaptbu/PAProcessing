clearvars -except emdd;

corrmin=1;
corrmax=5000;
highpass=00;
lowpass=250;
wallfilter=0;
emd_high=0;
for i=6:7
    disp(i);
    emd_low=i-1;
    figname=['emd_low_restricted',num2str(emd_low)];
    ReadFiles(highpass,lowpass,wallfilter,corrmin,corrmax,figname, emdd, emd_low, emd_high);
end
