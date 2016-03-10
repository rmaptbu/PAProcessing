function []= ReadFiles(highpass,lowpass,wallfilter,corrmin,corrmax,figname)
%% General Setup
options=containers.Map;
%Oscilloscope options
dt = 0.25; %sampling interval in ns
sampling_rate = 4000; %samples per microsecond
%timegating options
%Start at N,N+w -> N+q,N+w+q ->...->N+(jmax-1)q,N+(jmax-1)q+w
%==>> jmax is number of windows
%Ensure walking length needs to be integer multiple of stepisze
N=1000; %Starting point
q=50; %stepsize
w=200; %Interrogation windows
W=2000; %walking length
W=W-rem(W,q);
jmax=W/q;

%Plotting options
varying_separation = 0; %pulse sep the controlled variable?
varying_flowrate = 1; %or is the flow rate being changed?

%% Read files
% path for PA data
warning('Ensure filename for valid velocity measurements is correct');
basepath=...
    '/Users/Thore/Documents/MATLAB/PAProcessing/JosData/610nm/';
cd(basepath)

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
sz = [12 Inf];
files_raw=fscanf(fileID,formatSpec, sz); clearvars sz;
fclose(fileID);
%Convert files_raw into filenames
names={};
other_files={};
for i=1:size(files_raw,2)
    if files_raw(1,i) %flow rate is non zero
        minus='-';
    else
        minus='';
    end
    str=['rate_',minus,...
        num2str(files_raw(1,i)),'_0.5sep_250Mhz_12.764trig_',...
        num2str(files_raw(3,i)),'files_',...
        num2str(files_raw(2,i)),'.csv'];
    if i>1;
        %check if previous entries have same flowrate
        prevmatch1=[files_raw(1,1:i-1)==files_raw(1,i)];
        prevmatch2=[files_raw(3,1:i-1)==files_raw(3,i)];
        if any(prevmatch1);
            f1=find(prevmatch1); %index of non-zero elements            
            if any(prevmatch2);
                f2=find(prevmatch2); %index of non-zero elements 
                f3=find(ismember(f1,f2));                
                idx=f1(f3(1)); %idx is first occurence of same prev entry
                idx=map(idx); %index corresp. to pos in names{}
                str=[];
                %list of indexes for other files with from same flowrate
                other_files{idx}(length(f3))=files_raw(2,i);                  
            end
        end        
    end
    names=[names,str];
    map(i)=length(names);
end
clear('files_raw', 'str', 'minus', 'formatSpec', 'fileID',...
    'map', 'prevmatch1', 'prevmatch2', 'idx', 'f1', 'f2', 'f3', 'str');

%% analyse each file
number_of_files=length(names);
paf=PAF(dt, sampling_rate);
h = waitbar(0, 'Initialising Waitbar');
emd={};
for i=1:number_of_files;
    msg=['Analysing files...',num2str(number_of_files-i),' files left'];
    waitbar(i/number_of_files,h,msg);
    filename=names{i};
    %find flow rate from filename
    paf.ReadData(filename,other_files{i});
    
    if lowpass; paf.lowpass(lowpass);end
    if highpass; paf.highpass(highpass);end
    if wallfilter; paf.wallfilter();end
%     paf.emd(0,0);
%     emd{i}=paf.imf; %empirical mode decomposition
    paf.pressure=sum(emdd{i}{:}
    shift = paf.xcorr(1750, 2250);
    profile = paf.TimeGating(N, q, w, W, 0, 0);
    pressure = paf.pressure;
    
    shift_all(:,i)=shift';
    profile_all(:,i)=profile;
end
clos(h);
[shift_all I]=sortrows(shift_all');
shift_all=shift_all';

%% find mean of selected region of flowprofile
fmin=25;
fmax=35;
for i=1:number_of_files
    shift_profile(:,i) = [shift_all(1,i),mean(profile_all(fmin:fmax,i))];
    %fit parbolic flow profile
    p = polyfit(fmin:fmax,squeeze(profile_all(fmin:fmax,i))',2);
    xmin=-p(2)/(2*p(1));
    shift_poly(:,i) = [shift_all(1,i),p(1)*xmin^2+p(2)*xmin+p(3)];
end
%% create summary plots..
%shift_all contains the measured shift obtained by
%correlating the entire%waveform
shift_all_E=padarray(shift_all(2,:)', [0 size(profile_all,1)-1], ...
    'symmetric', 'post');
fig=figure('Visible','on');
%adjust the number in for loop to equal number of files anaylsed (n)
%make sure subplot has enough space: 
sbY=5;
sbX=5;
%make sure that you increase the subplot no as the no of files increase
for i=1:number_of_files; subplot(sbY,sbX,i);hold on; box on
    for j=1:jmax-1;
        plot(j:j+1,profile_all(j:j+1,I(i)),...
            'LineWidth',2,'Color',[j/jmax,.5-j/(jmax*2),1-j/jmax]);
    end;
    hold on
    plot(shift_all_E(i,:)', 'Color', [1 0 0]);
%     S=strrep(names(i),'_','\_');
%     S=S{1};
%     S=[S(1:21),'xcorr', num2str(shift_all_E(i,1))];
    title([num2str(shift_all(1,i)),' ml/s ',...
        num2str(length(other_files{i})),' file(s)']);
    ylim([-5 5])
    xlim([1 jmax])
    ax=gca;
    ax.YGrid = 'on';
    hold off;
end
%All profiles sumamry
subplot(sbY,sbX,i+1);hold on; box on
for j=1:number_of_files;
    plot(profile_all(:,I(j)),...
        'LineWidth',1,'Color',[j/number_of_files,...
        .5-j/(number_of_files*2),1-j/number_of_files]);
end;
title('all');
ylim([-5 5])
xlim([20 45])
ax=gca;
ax.YGrid = 'on';
hold off;

i=(ceil(number_of_files/sbX))*sbX; %go to start of a new row
subplot(sbY,sbX,i+1:i+2)
plot(pressure(:,1),'Color',[0 0 0]);hold on
for j=1:jmax-1;
    box on
    plot(N+j*q:N+j*q+49,pressure(N+j*q:N+j*q+49,1),...
        'LineWidth',2,'Color',[j/jmax,.5-j/(jmax*2),1-j/jmax]);
end;hold off;

subplot(sbY,sbX,i+3); box on; hold on;
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
    scatter(sep, shift_all,'Marker', '*','MarkerFaceColor', [0 0 0],...
        'MarkerEdgeColor', [0 0 0])
    %plot(sep_exp,sep_exp)
    set(gca, 'XGrid', 'on', 'YGrid', 'on')
    xlabel('Pulse Separation (ms)')
    ylabel('Measured Time Shift (ns)')
else %varying_flowrate
    plot(shift_all(1,:),shift_all(2,:),'k.-')
    xlabel('Known Rate (ml/h)');
    ylabel('Measured Time Shift (ns)');
    xlim([min(shift_all(1,:)) max(shift_all(1,:))]);
    
%     plot(shift_profile(1,:),shift_profile(2,:),'b.-')
%     plot(shift_poly(1,:),shift_poly(2,:),'r.-')
end
hold off
    str={['highpass = ', num2str(highpass),'Mhz'],...
        ['lowpass = ', num2str(lowpass),'Mhz'],...
        ['wallfilter = ', num2str(wallfilter)],...
        ['corrmin = ', num2str(paf.corrmin)],...
        ['corrmax = ', num2str(paf.corrmax)]};
    annotation('textbox',[ 1-1/sbX 0.025 1/sbX 1/sbY],...
        'String',str, 'FitBoxToText', 'on');
%set size
set(fig, 'Position', [100 100 256*sbX 256*sbY]);
%% save files, return
%save file
disp(['saving figure ', figname]);
try
    cd([basepath,'figures/']);
catch
    warning('Creating directory');
    mkdir([basepath,'figures/']);
    cd([basepath,'figures/']);
end
set(gcf,'PaperPositionMode','auto')
print(fig,figname,'-dpng','-r300')
close(fig);
cd('/Users/Thore/Documents/MATLAB/PAProcessing/');