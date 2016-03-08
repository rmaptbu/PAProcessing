clear;
% path for PA data
%% General Setup
NumP = 5000; %number of points
%% load files
basepath=...
'/Users/Thore/Documents/MATLAB/PAProcessing/Speckle/0226_blood0.1Hct/';
cd(basepath)

%Read files that were accepted
fileID = fopen('velocitymeasurements_full_peak0_Interpolant_rejected.csv','r');
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
sz = [12 Inf];
files_raw=fscanf(fileID,formatSpec, sz); clearvars sz;
fclose(fileID);
%Convert files_raw into filenames
names={};
for i=1:size(files_raw,2)
    str=['rate_-',...
        num2str(files_raw(1,i)),'_1.0sep_full_12.340trig_',...
        num2str(files_raw(3,i)),'files_',...
        num2str(files_raw(2,i)),'.csv'];
    
    names=[names,str];
end
clear('files_raw', 'str', 'minus', 'formatSpec', 'fileID');

dt=0.25; %sampling interval in ns
number_of_files=length(names);

%% do carpet plots

for u=1:number_of_files;
    display(u);
    filename=names{u};
    
    %load all files into "all_points" (just some variable)
    all_points = csvread([filename,'']); 
    %whatever it says here needs to add up to the full name of the file
%     all_points=cat(1,all_points,csvread([filename,'2.csv']));
    %
    i=0;
    %take 1000 point, split them off, filter attach to "pressure"
    %step by 1000, repeat, until end of file
    %split all data inte separate acquistions (indext in dimension 2)
    for x=NumP:2*NumP:size(all_points,1)-NumP;
        i=i+1;     
        pressure(u,1:NumP,i)= all_points(x-(NumP-1):x,2);
    end 
    clear('all_points');
end
p=mean(pressure,3);
for i=1:number_of_files-1
    xc(i,:)=xcorr(p(i,:),p(i+1,:),'unbiased');
end
%% Plots
%titles={'0','25','50','75','100','125'};
figure;
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 11)
for i=1:number_of_files
    subplot(number_of_files,1,i);
    imagesc(squeeze(pressure(i,:,:))');
    %title(titles{i});
end
figure;hold on;
p=mean(pressure,3);
for i=1:1:number_of_files
        plot(p(i,:),'.-','LineWidth', 2, 'Color',...
            [i/number_of_files, .5-i/(number_of_files*2),...
            1-i/number_of_files]);
end