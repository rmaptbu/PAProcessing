function [shift, profile, pressure] = PAF(filename, options)
%unpack options into namespace
Keys=keys(options);
for key=Keys
    value=options(key{1});
    evalc([key{1},'=',num2str(value)]);
end
clearvars Keys key value options

%shift(flowrate, xcorr_peak_pos in ns)
%profile is timegating xcorr peak pos plotted against windows position
expr='-?[\d.]+_0\.5sep';
[start, fin]=regexp(filename,expr);
flow_rate = str2double(filename(start:fin-7));
clearvars start fin expr

%load all files into "all_points" (just some variable)
all_points = csvread([filename,'']);
%whatever it says here needs to add up to the full name of the file
%     all_points=cat(1,all_points,csvread([filename,'2.csv']));
%     all_points=cat(1,all_points,csvread([filename,'3.csv']));
%     all_points=cat(1,all_points,csvread([filename,'4.csv']));
%     all_points=cat(1,all_points,csvread([filename,'5.csv']));
%
i=0;
%take 5000 point, split them off, filter attach to "pressure"
%step by 5000, repeat, until end of file
%split all data inte separate acquistions (indext in dimension 2)
for x=5000:5000:size(all_points,1);
    i=i+1;    
    %Split file
    pressure(1:5000,i)=all_points(x-4999:x,2);    
end
clear('all_points','x');

%% Filtering
%lowpass
if lowpass
    %Filter
    Fp=lowpass; %enter in MHz
    Fp=Fp*2*pi/sampling_rate; %sampling in samples per microscond
    %stop band frequency in rad/sample
    Fst=lowpass*1.1; %in Mhz
    Fst=Fst*2*pi/sampling_rate;
    d = fdesign.lowpass('Fp,Fst,Ap,Ast',Fp,Fst,1,60);
    % % Ap ? amount of ripple allowed in the pass band in decibels
    % (the default units). Also called Apass.
    % % Ast ? attenuation in the stop band in decibels
    % (the default units). Also called Astop.
    % % F3db ? cutoff frequency for the point 3 dB point below the
    % passband value. Specified in normalized frequency units.
    % % Fc ? cutoff frequency for the point 6 dB point below the
    % passband value. Specified in normalized frequency units.
    % % Fp ? frequency at the start of the pass band. Specified in
    % normalized frequency units. Also called Fpass.
    % % Fst ? frequency at the end of the stop band. Specified in
    % normalized frequency units. Also called Fstop.
    % % N ? filter order.
    % % Na and Nb are the order of the denominator and numerator.
    hd = design(d,'butter');
    for i=1:size(pressure,2)
        pressure(:,i) = filter(hd,pressure(:,i));
        %             d1r = wrev(d1);
        %             d2 = filter(hd,d1r);
        %             pressure_lowp(:,i) = wrev(d2);
    end
end 
%highpass
if highpass
    %pass band frequency in rad/sample
    Fp=highpass; %enter in MHz
    Fp=Fp*2*pi/sampling_rate; %sampling in samples per microscond
    %stop band frequency in rad/sample
    Fst=highpass/2; %in Mhz
    Fst=Fst*2*pi/sampling_rate;
    d = fdesign.highpass('Fst,Fp,Ast,Ap',Fst,Fp,60,1);
    hd = design(d,'butter');
    for i=1:size(pressure,2)
        pressure(:,i) = filter(hd,pressure(:,i));
        %             d1r = wrev(d1);
        %             d2 = filter(hd,d1r);
        %             pressure_lowp(:,i) = wrev(d2);
    end
end
%wallfilter
if wallfilter
    meansignal=squeeze(mean(pressure,2));
    pressure=pressure-repmat(meansignal,[1 size(pressure,2)]);
end
clearvars hd d Fst Fp highpass lowpass wallfilter
%% XCORR
%cross correlate all signal and normalise by maximum

width=corrmax-corrmin;
for i=1:2:size(pressure,2)-1;
    xcorrs(:,(i+1)/2)=xcorr(pressure(corrmin:corrmax,i+1),...
        pressure(corrmin:corrmax,i),'biased');
    xcorrs_norm(:,(i+1)/2)=xcorrs(:,(i+1)/2)./max(xcorrs(:,(i+1)/2));
    %xcorrs_norm: normalised xcorr of all pairs
    %xcorrs_norm(xcorr_value, Pair_index)
end
clearvars xcorrs corrmin corrmax

%% find maximum of x-correlation of entire waveform
%take average of all normalised cross correlations and find position of
%maximum
%"shift" is position of maximum
[~, shift_i]=max(xcorrs_norm);
shift_std=std(shift_i); clear('shift_i');
%ensemble correlation
meanxcorr=squeeze(mean(xcorrs_norm,2));
[~, shift]=max(meanxcorr);

%Interpolate location of maximum
t = -width:width; 
t = t * dt;%map of time points
%create 100+100 points around maximum data point
%if I=1 then peak is at centre. i.e. need to shift by one
%to get rid of offset (peak at I-1 in absolute time)
tI = (shift-width-2):1.0/100:(shift-width);
tI = tI * dt;
xcorr_interp = ... %100 point interpolation
    interp1(t,meanxcorr,tI,'spline');
[~,max_pos]=max(xcorr_interp');
xcorr_peak_pos=tI(max_pos);


shift=(shift-width)*dt;
shift=[flow_rate, shift];
%SHIFT_ALL=cat(1,SHIFT_ALL,xcorr_peak_pos);

clear('max_shift_i','meanxcorr', 'xcorrs');
%% Time gating
jmax=W/q;
assert(~rem(W,q),'W/q not integer');
if time_gating
    %starting point for xcorr
    corr_bounds=ceil(0.5*w); 
    %between 0and1.defines size of search for corr peak
    %assume first pair correlates
    for i=1:2:size(pressure,2)-1;
        for j=1:jmax
            xcorrwindows(1:2*w+1,(i+1)/2,j)=xcorr(...
                pressure((N+(j-1)*q):(N+(j-1)*q+w),i+1),...
                pressure((N+(j-1)*q):(N+(j-1)*q+w),i),'biased');
        end
    end
    %Normalise
    if normalise
        for i=1:length(xcorrwindows(1,:,1))
            for j=1:length(xcorrwindows(1,1,:))
                xcorrwindows(:,i,j)=...
                    xcorrwindows(:,i,j)./max(xcorrwindows(:,i,j));
            end
        end
    end
    
    xcorrwindows_ensemble=squeeze(mean(xcorrwindows,2));
    %w+1 is centre
    [M, I]=max(xcorrwindows_ensemble(w+1-corr_bounds:w+1+corr_bounds,:),...
        [],1);
    %adjust I to work with indexing of xcorrwindows_ensemble
    %coor_bounds+1 is at centre. needs to be shifted to w+1
    I=I+w-corr_bounds;
    clearvars xcorrwindows
%     if remove_outliers
%         for i=1:length(I);
%             if I(i)<1 || I(i)>w*0.6 %peak out of bounds
%                 %play with 1 and 0.6 to adjust bound for outliers
%                 if i~=1
%                     %replace outlier by previous data point
%                     I_0(i)=I_0(i-1);
%                 else
%                     I_0(i)=0;
%                 end
%                 display(['changing ',num2str(i)])
%             else
%                 I_0(i)=I(i);
%             end
%         end
%         I=I_0;
%     end
    offset=0;%w-1-corr_bounds; %legacy bit. ignore for now
    I=I';
    
    %Interpolate location of maximum
    t=-w:w;
    t=t*dt;
    %map of time points
    for i=1:size(I,1)
        %create 100+100 points around maximum data point
        %if I=1 then peak is at centre. i.e. need to shift by one
        %to get rid of offset (peak at I-1 in absolute time)
        tI(i,:) = (I(i)-w-2):1.0/100:(I(i)-w);
        tI(i,:) = dt*tI(i,:);
        xcorr_interp(i,:) = ... %100 point interpolation
            interp1(t,xcorrwindows_ensemble(:,i),tI(i,:),'spline');
    end
    [~,max_pos]=max(xcorr_interp');
    for i=1:size(max_pos,2)
        xcorr_peak_pos(i)=tI(i,max_pos(i));
    end
    
    profile=xcorr_peak_pos';
else
    profile = 0;
end

%% PLOTS
if draw
    figure; hold on
    set(0,'DefaultAxesFontName', 'Times New Roman')
    set(0,'DefaultAxesFontSize', 11)
    for i=1:jmax-1;
        plot(i:i+1,I(i:i+1),'LineWidth',2,'Color',...
            [i/jmax,.5-i/(jmax*2),1-i/jmax]);
    end ;
    %     normal xcorr plot
    % overview plot
    figure
    nplots=4;
    ppg=floor(jmax/nplots); %plots per graph
    for j=0:nplots-1
        subplot(4,2,5+j)
        %         subplot(3,5,j+1)
        hold on;
        xlim([w w+100]);
        for i=j*ppg+1:(j+1)*ppg;
            plot(squeeze(xcorrwindows_ensemble(:,i)),'black');
        end
        %max points
        
        for i=j*ppg+1:(j+1)*ppg-1;
            plot(I(i:i+1)+offset,M(i:i+1),'b-o',...
                'LineWidth',2,'Color',[i/jmax,.5-i/(jmax*2),1-i/jmax]);
        end
        hold off
    end
    
    %Mean signal
    subplot(4,2,1:2)
    meansignal=squeeze(mean(pressure,2));
    plot(meansignal);
    hold on;
    for i=1:jmax-1;
        plot(N+i*q:N+i*q+49,meansignal(N+i*q:N+i*q+49),...
            'LineWidth',2,'Color',[i/jmax,.5-i/(jmax*2),1-i/jmax]);
    end ;
    
    %all xcorrs
    subplot(4,2,3:4);
    %%
    figure;hold on
    for i=1:jmax-1;
        fig1=plot(squeeze(xcorrwindows_ensemble(:,i)./...
            max(abs(xcorrwindows_ensemble(:,i)))),...
            'Color',[.2,.2,.2]);
        fig1.Color(4)=0.1;
    end
    hold off
    %%
    %         figure;plot(maxprofile);
    %     savefig([filename,'.fig']);
    
    figure;
    subplot(1,2,1);hold on
    plot(pressure(:,1));
    for i=1:jmax-1;
        plot(N+i*q:N+i*q+49,pressure(N+i*q:N+i*q+49,1),...
            'LineWidth',2,'Color',[i/jmax,.5-i/(jmax*2),1-i/jmax]);
    end ;hold off
    subplot(1,2,2);hold on
    plot(pressure(:,2));
    for i=1:jmax-1;
        plot(N+i*q:N+i*q+49,pressure(N+i*q:N+i*q+49,2),...
            'LineWidth',2,'Color',[i/jmax,.5-i/(jmax*2),1-i/jmax]);
    end ;
end