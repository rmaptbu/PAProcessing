%% General Setup
clear;
options=containers.Map;
%Oscilloscope options
options('dt')=0.25; %sampling interval in ns
options('sampling_rate') = 4000; %samples per microsecond
%Filtering options
options('highpass') = 10; %high pass filter?(freq in Mhz) 0=no filter
options('lowpass') = 250; %low pass filter?(freq in Mhz) 0=no filter
options('wallfilter') = 0; %remove mean signal from all signals?
%Xcorr of entire waveform options
options('corrmin')=1750;
options('corrmax')=2250;
%Time Gating options
options('time_gating')=1;
options('draw')=0; %only works with time_gating=1
options('remove_outliers')=1; %currently not implemented
options('N')=1000;%starting poing
options('q')=30;%stepsize
options('w')=200; %interrogation windows (actual length is w-1)
options('W')=2000; %walking length
N=options('N');
q=options('q');
w=options('w');
W=options('W');
jmax=ceil(W/q);

%Plotting options
varying_separation = 0; %pulse sep the controlled variable?
varying_flowrate = 1; %or is the flow rate being changed?

%% Read files
% path for PA data
basepath=...
    '/Users/Thore/Documents/MATLAB/PAProcessing/JosData/raw/';
cd(basepath)
clearvars basepath

%Read files that were accepted
fileID = fopen('velocitymeasurements_250MHz_peak0_Interpolant.csv','r');
formatSpec = '%f %f %f %f %f %f %f %f %f %f %f %f Good result';
% 1   Dial/rate (ml/h)
% 2   Record number
% 3   Number of files
% 4   Known velocity (mm/s)
% 5   Known resolution (mm/s)
% 6   Measured velocity from fit to xcorr peak (mm/s)
% 7   Measured resolution from fit to xcorr peak (mm/s)
% 8   Measured pulse separation (ms)
% 9   Theoretical resolution based on oscilloscope sampling interval (mm/s)
% 10  Amplitude of selected peak
% 11  Amplitude of selected peak relative to the next largest peak
% 12  Amplitude of selected peak relative to the RMS of the xcorr amplitude
% 13  Quality of result (options: "Good result", Default; "Bad result
% (incorrect pulse sep - laser misfiring)"; "Bad
size = [12 Inf];
files_raw=fscanf(fileID,formatSpec, size); clearvars size;
fclose(fileID);

%Convert files_raw into filenames
names={};
for i=1:size(files_raw,2)
    if files_raw(1,i)
        minus='-';
    else
        minus='';
    end
    if files_raw(2,i)==1
        str=['rate_',minus,...
            num2str(files_raw(1,i)),'_0.5sep_250Mhz_12.764trig_',...
            num2str(files_raw(3,i)),'files_',...
            num2str(files_raw(2,i)),'.csv'];
        names=[names,str];
    end
end
clearvars files_raw str minus formatSpec fileID;

%% analyse each file
number_of_files=length(names);
for i=1:number_of_files;
    display(i);
    filename=names{i};
    %find flow rate from filename
    [shift, profile, pressure]=PAF(filename, options);
    shift_all(:,i)=shift;
    profile_all(:,i)=profile;
end

%% create summary plots..
%shift_all contains the measured shift obtained by
%correlating the entire%waveform
shift_all_E=padarray(shift_all(2,:)', [0 size(profile_all,1)-1], ...
    'symmetric', 'post');
figure;
%adjust the number in for loop to equal number of files anaylsed (n)
%make sure subplot has enough space: 
sbX=5;
sbY=5;
%make sure that you increase the subplot no as the no of files increase
for i=1:number_of_files; subplot(sbX,sbY,i);hold on; box on
    for j=1:jmax-1;
        plot(j:j+1,profile_all(j:j+1,i),...
            'LineWidth',2,'Color',[j/jmax,.5-j/(jmax*2),1-j/jmax]);
    end;
    hold on
    plot(shift_all_E(i,:)', 'Color', [1 0 0]);
    S=strrep(names(i),'_','\_');
    S=S{1};
    S=[S(1:21),'xcorr', num2str(shift_all_E(i,1))];
    title(S);
    ylim([0 5])
    hold off;
end
subplot(sbX,sbY,i+1:i+2)
plot(pressure(:,1),'Color',[0 0 0]);hold on
for j=1:jmax-1;
    box on
    plot(N+j*q:N+j*q+49,pressure(N+j*q:N+j*q+49,1),...
        'LineWidth',2,'Color',[j/jmax,.5-j/(jmax*2),1-j/jmax]);
end;hold off;

if varying_separation
    %string match separation
    expr='rate_-[\d.]+_[\d.]+';
    [a, index]=regexp(names,expr);
    for i=1:size(index,1)
        sep(i) = str2double(names{i}(12:index{i,1}));
    end
    %sep in ms, xcorr peak in ns
    %dT=T*1.04E-7
    sep_exp=sep*1.04E-1;

    subplot(sbX,sbY,i+3)
    scatter(sep, shift_all,'Marker', '*','MarkerFaceColor', [0 0 0],...
        'MarkerEdgeColor', [0 0 0])
    hold on
    box on
    %plot(sep_exp,sep_exp)
    set(gca, 'XGrid', 'on', 'YGrid', 'on')
    xlabel('Pulse Separation (ms)')
    ylabel('Measured Time Shift (ns)')
end
if varying_flowrate
    subplot(sbX,sbY,i+3)
    scatter(shift_all(1,:),shift_all(2,:))
    xlabel('Known Rate (ml/h)');
    ylabel('Measured Time Shift (ns)');
    box on;
end
cd('/Users/Thore/Documents/MATLAB/PAProcessing/');