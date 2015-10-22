%% DESCRIPTION

% Select files to evaluate
% Calculate cross-correlation for all waveform pairs, and mean xcorr for
% each file, and the time shift for each file
% Calculate FFT for each waveform
% Take the mean FFT for all waveforms, and for each velocity
% Calculate the mean time shift for each velocity: 
%           this is the mean of multiple files recorded at the same
%           velocity
%           It is also possible just to analyse a single file, but then
%           there will be zero error bar (the error bar is the standard
%           deviation of multiple measurements at the same velocity)
% Convert to mm/s and write the data to a file

% NB: requires the lbox2 function
% lbox should be saved in the same folder as this script
% The Matlab current directory should contain the data to be analysed
% When running this script, click "Add to path" when prompted

clear
close all

%% OPTIONS AND INPUTS

framesPerFile = 100;
extreme = [Inf,250]; %two cases of extremes to exclude: Inf = exclude none; 250 = exclude measured velocities>250 mm/s
vlim = 50; %calculate means for all values below vlim

% Filtering
filterData=true;
lows=250;%[5,10,20,30];
highs=250;%[5,10,20,30,50,80];%[5:5:100]; %
filterStructure = 'lowpass'; %choose lowpass, highpass, bandpass
filterType = 'doubleFiltering'; %choose filter (i.e. conventional), filtfilt (zero-phase), doubleFiltering (also zero-phase)

runAutomatically = false; %analyses files one after the other with no user input
split = false; %splits files into smaller parts e.g. if working with a 50 frame file, split=2 gives 2x 25 frame files (only the first 24 will be analysed because need an even number of frames)
plotXcorr = true; %option to plot the individual cross-correlation functions
analyse_velocities=true; %analyse the calculated velocities
plotVspectra = true; %option to plot the mean frequency spectrum for each velocity
plotAnalysis = true; %option to plot the results of the analysis: mean fractional error and mean resolution
savePlots = false; %option to save all the plots


%% START CALCULATIONS

for Fi=1:length(lows)
    
    %%%%%%%%%% Prepare data for analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if runAutomatically
        reply = 'Y';
    else
        reply = input([num2str(framesPerFile),' frames per file? Type "N" to cancel, or any key to continue\n'],'s');
    end
    
    if isempty(reply)
        reply = 'Y';
    end
    
    if reply~='N'
        
        if runAutomatically
            reply = 'N';
        else
            reply = input('Select files to evaluate (via lbox2 function)?\nType "Y" to accept, or any key to continue with current file selection\n','s');
        end
        
        if isempty(reply)
            reply = 'N';
        end
        if reply=='Y'
            lbox2 %select files to evaluate
            replyC = input('Select files from Directory List, then press any key to continue.\n(Alternatively, type N to cancel)\n', 's');
            if isempty(replyC)
                replyC = 'Y';
            end
        else replyC = 'Y'; %don't cancel
        end
        
        mac=1; %using mac

        try

            if mac

                filenames=readtable('files_meanFFT.txt');

                filenames=table2cell(filenames);

            else

                [N,filenames]=xlsread('files_meanFFT');

            end           

        catch

            file_error=lasterr

        end
        
        if split
            splits = 2; %number of parts to split file into
        else
            splits = 1;
        end
        splitNo = 1;
        filenames = repmat(filenames,1,splits);
        filenames = reshape(filenames',size(filenames,1)*splits,1);
        
        if rem(framesPerFile,2)~=0 %odd number of frames: exclude last frame
            evenFramesPerFile = framesPerFile-1;
        else evenFramesPerFile = framesPerFile; %(framesPerFile-4)./splits;
            %usually set evenFramesPerFile = framesPerFile to include all frames
            %If splitting file: set evenFramesPerFile to number of frames
            %in each file part e.g. framesPerFile./splits
        end
        
        ERRno=1; col=1;
        numberFiles = length(filenames)
        scrsz=get(0,'ScreenSize');
        if plotXcorr
            figXcorr = figure('Position',[0 0 scrsz(3) scrsz(4)]); sub=1; figNO=1;
        end
        ts = zeros(numberFiles,2);
        for fileNo=1:length(filenames)
            
            if replyC == 'N'
                break
            end
            
            fileNo
            try
                filename=filenames(fileNo);
                clear data
                data=dlmread(filename{1});
            catch
                error_number=ERRno
                ERRno=ERRno+1;
                error_message=lasterr
                continue
            end
            
            %%%%%%%%%% Filter data %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            dt = data(2,1)-data(1,1);
            str_dt = num2str(dt,'%1.4e'); %set to 4 d.p.
            dt = str2num(str_dt); %sampling interval in microseconds
            nopoints = length(data)./framesPerFile; %no. of points per frame
            f = 1/(dt*10^-6)*(0:nopoints/2)/(nopoints); %frequency axis in Hz
            plotHz = 10^6; %plot Hz / kHz / MHz e.g. 1000=kHz, 10^6=MHz
            if plotHz==1
                unit='Hz';
            elseif plotHz==1000
                unit='kHz';
            elseif plotHz==1000000
                unit='MHz';
            end
            plotprop = 10; %proportion of points to plot (must be greater than or equal to 2)
            
            t = dt*(nopoints-1):-dt:-dt*(nopoints-1); %vector of time shifts for xcorr (microseconds)
            
            data = data(evenFramesPerFile*(splitNo-1)*nopoints+1:evenFramesPerFile*splitNo*nopoints,:);
            if splitNo==splits
                splitNo = 1; %start again
            else splitNo = splitNo+1;
            end
            
            d = zeros(nopoints,2);
            c = zeros(length(t),evenFramesPerFile/2);
            for frame=1:evenFramesPerFile
                if rem(frame,2)==0 %second frame in a pair
                    xcol = 2;
                else xcol = 1;
                    d = zeros(nopoints,2);
                end
                dd = data((frame-1)*nopoints+1:frame*nopoints,2);
                if filterData
                    sampling_rate = 1./dt; %sampling rate in MHz
                    order=2;
                    switch filterStructure
                        case 'bandpass'
                            lowF = lows(Fi); %lower cutoff in MHz
                            highF = highs(Fi); %upper cutoff in MHz
                            [b,a]=butter(order,[lowF/sampling_rate*2,highF/sampling_rate*2]); %bandpass Butterworth filter
                            
                        case 'lowpass'
                            lowF = lows(Fi); %lower cutoff in MHz
                            [b,a]=butter(order,lowF/sampling_rate*2,'low');
                            
                        case 'highpass'
                            highF = highs(Fi); %upper cutoff in MHz
                            [b,a]=butter(order,highF/sampling_rate*2,'high');
                    end
                    switch filterType
                        case 'filter'
                            d(:,xcol) = filter(b,a,dd);
                        case 'filfilt'
                            d(:,xcol) = filtfilt(b,a,dd);
                        case 'doubleFiltering'
                            d1 = filter(b,a,dd);
                            d1r = wrev(d1);
                            d2 = filter(b,a,d1r);
                            d2r = wrev(d2);
                            d(:,xcol) = d2r;
                    end
                    switch filterStructure
                        case 'bandpass'
                            switch filterType
                                case 'filter'
                                    titleF = {['MATLAB filtered ',num2str(lowF),'-',num2str(highF),' MHz, Butterworth order ',num2str(order)]}; %for fig title
                                    fileF = ['data_filtered_butter',num2str(lowF),'-',num2str(highF),'MHzOrder',num2str(order)]; %for filename
                                case 'doubleFiltering'
                                    titleF = {['MATLAB zero-phase filtered (double filtering)'],[num2str(lowF),'-',num2str(highF),' MHz, Butterworth order ',num2str(order)]}; %for fig title
                                    fileF = ['data_zero-phase-filtered_butter',num2str(lowF),'-',num2str(highF),'MHzOrder',num2str(order)]; %for filename
                            end
                        case 'lowpass'
                            switch filterType
                                case 'filter'
                                    titleF = {['MATLAB filtered -',num2str(lowF),' MHz, Butterworth order ',num2str(order)]}; %for fig title
                                    fileF = ['data_filtered_butter-',num2str(lowF),'MHzOrder',num2str(order)]; %for filename
                                case 'doubleFiltering'
                                    titleF = {['MATLAB zero-phase filtered (double filtering)'],[num2str(lowF),' MHz, lowpass Butterworth order ',num2str(order)]}; %for fig title
                                    fileF = ['data_zero-phase-filtered_butter-',num2str(lowF),'MHzOrder',num2str(order)]; %for filename
                            end
                        case 'highpass'
                            switch filterType
                                case 'doubleFiltering'
                                    titleF = {['MATLAB zero-phase filtered (double filtering)'],[num2str(highF),' MHz, highpass Butterworth order ',num2str(order)]}; %for fig title
                                    fileF = ['data_zero-phase-filtered_butter+',num2str(highF),'MHzOrder',num2str(order)]; %for filename
                            end
                    end
                else d(:,xcol) = dd;
                    titleF = {'No MATLAB filtering'};
                    fileF = 'data_unfiltered';
                end
                
                %%%%%%%%%% Calculate cross-correlation %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if rem(frame,2)==0 %second frame in a pair
                    xc = xcorr(d(:,1),d(:,2),'unbiased');
                    c(:,col/2) = xc./max(xc); %normalise
                    [C,I]=max(c(:,col/2));
                    if plotXcorr
                        figure(figXcorr); subplot(5,4,sub); sub=sub+1;
                        plot(t,c(:,col/2)); title([num2str((figNO-1)*20+sub-1),'; ts = ',num2str(t(I)*1000),' ns'])
                        if sub==21
                            if savePlots; saveas(figXcorr,[datestr(now,'ddmmmyy_HH-MM-SS-AM'),'_',fileF,'_meanXcorr(',num2str(figNO),').fig']); end
                            figXcorr = figure('Position',[0 0 scrsz(3) scrsz(4)]); sub=1; figNO=figNO+1;
                        end
                    end
                end
                D = fft(d(:,xcol),nopoints);
                PD = D.* conj(D) / nopoints; %normalise for number of points
                PD = PD./max(PD(1:nopoints/2+1)); %normalise for amplitude
                P(:,col) = PD; %power spectra
                col=col+1;
            end
            c = c(:,((col-1)/2-evenFramesPerFile/2+1:(col-1)/2)); %select relevant columns for file just analysed
            meanXcorr = mean(c,2);
            interpolate = true;
            if interpolate
                tI = dt*(nopoints-1):-dt/10:-dt*(nopoints-1); %vector of time shifts for xcorr (microseconds)
                meanXcorrI = interp1(t,meanXcorr,tI,'spline'); %10 point interpolation
            else
                meanXcorrI = meanXcorr; tI = t;
            end
            searchLims = [1,0.8]; %percentage of xcorr function to search
            centre = ceil(length(meanXcorrI)/2);
            limits = [centre-floor(searchLims/2*length(meanXcorrI)); centre+floor(searchLims/2*length(meanXcorrI))];
            for sL = 1:2
                [C,I] = max(meanXcorrI(limits(1,sL):limits(2,sL)));
                ts(fileNo,sL) = tI(I+limits(1,sL)-1)*1000; %time shift in ns
            end
            if plotXcorr
                figure(figXcorr); subplot(5,4,sub); sub=sub+1;
                plot(t,meanXcorr);
                xlabel('Time shift (\mus)')
                title({filename{1},['Time shift = ',num2str(ts(fileNo,1)),' ns (search ',num2str(searchLims(1)),'), ',num2str(ts(fileNo,2)),' ns (search ',num2str(searchLims(2)),')']},'Interpreter','none','FontSize',8)
                
            end
        end
        if and(plotXcorr,savePlots)
            saveas(figXcorr,[datestr(now,'ddmmmyy_HH-MM-SS-AM'),'_',fileF,'_meanXcorr(',num2str(figNO),').fig'])
        end
        
        fileErrors = ERRno-1
        
        s1=size(P,2);
        nans = find(isnan(sum(P)));
        P(:,nans)=[]; %delete any NaN columns in P
        NaNfiles = (s1-size(P,2))./evenFramesPerFile
        meanP = mean(P,2);
        stdP = std(P,0,2);
        nanFrames=(nans-1)/evenFramesPerFile; %NaN frames
        nanFileNos = nanFrames(1:25:length(nanFrames)); %numbers of NaN files
        
        %%%%%%%%%% Plots the frequency spectra %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure('Position',[scrsz(3)/3 0 scrsz(3)/2 scrsz(4)*0.8])
        subplot(211)
        errorbar(f(1:nopoints/plotprop)./plotHz,meanP(1:nopoints/plotprop),stdP(1:nopoints/plotprop),'.')
        xlabel(['Frequency (',unit,')'])
        ylabel('Power')
        title({pwd,datestr(now),['Calculations assume ',num2str(framesPerFile),...
            ' frames per file, giving each frame ',num2str(nopoints),...
            ' points, and dt = ',num2str(dt*1000),' ns'],titleF{:}},'Interpreter','none','FontSize',8)
        
        subplot(212)
        plot(f(1:nopoints/plotprop)./plotHz,meanP(1:nopoints/plotprop),'Marker','.')
        xlabel(['Frequency (',unit,')'])
        ylabel('Power')
        title('No error bars')
        
        if savePlots
            saveas(gcf,[datestr(now,'ddmmmyy_HH-MM-SS-AM'),'_meanFFT.fig'])
        end
        
        xlswrite([datestr(now,'ddmmmyy_HH-MM-SS-AM'),'_files_meanFFT'],filenames)
        % save the names of the files that have been analysed
        
        
        %% Analyse velocities?
        if analyse_velocities
            
            tank_data = dlmread('tank_data.txt');
            radius = tank_data(3)*10^-4/2;
            soundSpeed = tank_data(6); %sound speed in m/s
            T = tank_data(9); %pulse separation in ms
            theta = tank_data(5)/180*pi; %theta in radians
            %theta = 58/180*pi;
            
            manualShifting = false;
            if manualShifting
                rateV_factor = 1000./T; %conversion from shifts in mm to V in mm/s
                xlabelS = 'Shifts (mm)';
            else
                rateV_factor = 10/(3600*pi*radius^2); %conversion from rates in ml/hr to V in mm/s
                xlabelS = 'Rates (ml/hr)';
            end
            factor = soundSpeed.*10^-9./(T.*cos(theta))*10^6; %conversion from ts in ns to V in mm/s
            
            filenames(nanFileNos)=[];
            if plotVspectra
                fig1=figure('Position',[0 0 scrsz(3) scrsz(4)]);
                fig2=figure('Position',[0 0 scrsz(3) scrsz(4)]);
            end
            v=1; sv=0; m=1; sub=1; figNO=1; rep=zeros(length(filenames),1);
            rates = zeros(size(ts,1),1);
            
            while v<=length(filenames)
                fileS = filenames(sv+1); fileS=fileS{1}; %start file
                if v~=length(filenames)
                    file1 = filenames(v); file1=file1{1}; underScores1=strfind(file1,'_'); dot1=strfind(file1,'.');
                    r1 = str2double(file1(underScores1(1)+1:underScores1(2)-1)); %find rate (ml/hr)
                    rep(v) = str2double(file1(underScores1(length(underScores1))+1:dot1(length(dot1))-1)); %find record number
                    file2 = filenames(v+1); file2=file2{1}; underScores2=strfind(file2,'_'); dot2=strfind(file2,'.');
                    r2 = str2double(file2(underScores2(1)+1:underScores2(2)-1)); %find rate (ml/hr)
                    rep(v+1) = str2double(file2(underScores2(length(underScores2))+1:dot2(length(dot2))-1)); %find record number
                else
                    if length(filenames)==1
                        file1 = filenames(v); file1=file1{1}; underScores1=strfind(file1,'_'); dot1=strfind(file1,'.');
                        r1 = str2double(file1(underScores1(1)+1:underScores1(2)-1)); %find rate (ml/hr)
                        rep(v) = str2double(file1(underScores1(length(underScores1))+1:dot1(length(dot1))-1)); %find record number
                        r2 = r1;
                        file2 = file1;
                    else
                        sv=v; m=m+1;
                    end
                    
                end
                if r1~=r2 || v==length(filenames)-1 || length(filenames)==1
                    if and(v==length(filenames)-1,r1==r2); v=v+1; file1=file2; end
                    r(m) = r1;
                    rates(sv+1:v,1) = r1;
                    relevant_ts = ts(sv+1:v,:);%mean time shifts, two searching limits
                    ex=1;
                    for extremes = [extreme(1),extreme(2)]./factor
                        meants(m,1,ex) = mean(relevant_ts(abs(relevant_ts(:,1))<extremes,1));
                        stdts(m,1,ex) = std(relevant_ts(abs(relevant_ts(:,1))<extremes,1));
                        meants(m,2,ex) = mean(relevant_ts(abs(relevant_ts(:,2))<extremes,2));
                        stdts(m,2,ex) = std(relevant_ts(abs(relevant_ts(:,2))<extremes,2));
                        ex=ex+1;
                    end
                    
                    meanPv(:,m) = mean(P(:,sv*evenFramesPerFile+1:v*evenFramesPerFile),2); %mean power spectrum
                    stdPv(:,m) = std(P(:,sv*evenFramesPerFile+1:v*evenFramesPerFile),0,2);
                    if plotVspectra
                        if sub==21
                            if savePlots
                                saveas(fig1,[datestr(now,'ddmmmyy_HH-MM-SS-AM'),'_meanFFT_velocities_errorBars(',num2str(figNO),').fig'])
                                saveas(fig2,[datestr(now,'ddmmmyy_HH-MM-SS-AM'),'_meanFFT_velocities_NOerrorBars(',num2str(figNO),').fig'])
                            end
                            fig1=figure('Position',[0 0 scrsz(3) scrsz(4)]);
                            fig2=figure('Position',[0 0 scrsz(3) scrsz(4)]);
                            sub=1; figNO=figNO+1;
                        else
                            figure(fig1); subplot(5,4,sub);
                            errorbar(f(1:nopoints/plotprop)./plotHz,meanPv(1:nopoints/plotprop,m),stdPv(1:nopoints/plotprop,m),'.')
                            xlabel(['Frequency (',unit,')'])
                            ylabel('Power')
                            title({[fileS,' to'],file1},'FontSize',8,'Interpreter','none')
                            %xlim([0 60])
                            
                            figure(fig2); subplot(5,4,sub); sub=sub+1;
                            plot(f(1:nopoints/plotprop)./plotHz,meanPv(1:nopoints/plotprop,m))
                            xlabel(['Frequency (',unit,')'])
                            ylabel('Power')
                            title({[fileS,' to'],file1},'FontSize',8,'Interpreter','none')
                            %xlim([0 60])
                        end
                    end
                    sv=v; m=m+1;
                end
                v=v+1;
            end
            if and(plotVspectra,savePlots)
                saveas(fig1,[datestr(now,'ddmmmyy_HH-MM-SS-AM'),'_meanFFT_velocities_errorBars(',num2str(figNO),').fig'])
                saveas(fig2,[datestr(now,'ddmmmyy_HH-MM-SS-AM'),'_meanFFT_velocities_NOerrorBars(',num2str(figNO),').fig'])
            end
            
            plotAbsolute = false;
            
            V = rates.*rateV_factor.*sign(rates).^(not(plotAbsolute)+1); %known velocities in mm/s
            mean_V = r.*rateV_factor.*sign(r).^(not(plotAbsolute)+1); %known velocities in mm/s
            
            ts = ts.*sign(ts).^(not(plotAbsolute)+1); %choose sign according to whether plotting absolute values
            meants = meants.*sign(meants).^(not(plotAbsolute)+1);
            
            measV = factor.*ts; %measured velocities in mm/s
            mean_measV = factor.*meants; %measured velocities in mm/s
            std_measV = factor.*stdts; %measured standard deviations in mm/s
            
            ex = 1;
            for extremes = [extreme(1),extreme(2)]
                figure('Position',[0 0 scrsz(3) scrsz(4)]); sub=-1;
                for sL=1:2 %plot for two different searching limits
                    sub=sub+2;
                    subplot(2,4,sub); plot(V./factor,V./factor,'k'); hold on; plot(V(abs(ts(:,sL))<(extremes./factor))./factor,ts(abs(ts(:,sL))<(extremes./factor),sL),'.'); xlabel('Known time shift (ns)'); ylabel('Measured time shift (ns)');
                    pos=get(gca,'Position'); set(gca,'Position',[pos(1) pos(2) pos(3) pos(4)*0.8])
                    axes('Position',get(gca,'Position'),'YTick',[],'XAxisLocation','top','Color','none','box','on','XLim',(get(gca,'XLim').*factor)./rateV_factor); xlabel(xlabelS); title({'All time shifts',['Searching ',num2str(searchLims(sL)*100),'% xcorr; Excluding |V| > ',num2str(extremes),' mm/s'],titleF{:}},'FontSize',8)
                    subplot(2,4,sub+1); plot(V,V,'k'); hold on; plot(V(abs(measV(:,sL))<extremes),measV(abs(measV(:,sL))<extremes,sL),'.'); xlabel('Known velocity (mm/s)'); ylabel('Measured velocity (mm/s)');
                    pos=get(gca,'Position'); set(gca,'Position',[pos(1) pos(2) pos(3) pos(4)*0.8])
                    axes('Position',get(gca,'Position'),'YTick',[],'XAxisLocation','top','Color','none','box','on','XLim',get(gca,'XLim')./rateV_factor); xlabel(xlabelS); title({'All velocities',['Searching ',num2str(searchLims(sL)*100),'% xcorr; Excluding |V| > ',num2str(extremes),' mm/s'],titleF{:}},'FontSize',8)
                    subplot(2,4,sub+4); plot(V./factor,V./factor,'k'); hold on; errorbar(mean_V'./factor,meants(:,sL,ex),stdts(:,sL,ex),'.'); xlabel('Known time shift (ns)'); ylabel('Measured time shift (ns)');
                    pos=get(gca,'Position'); set(gca,'Position',[pos(1) pos(2) pos(3) pos(4)*0.8])
                    axes('Position',get(gca,'Position'),'YTick',[],'XAxisLocation','top','Color','none','box','on','XLim',(get(gca,'XLim').*factor)./rateV_factor); xlabel(xlabelS); title({'Mean time shifts',['Searching ',num2str(searchLims(sL)*100),'% xcorr; Excluding |V| > ',num2str(extremes),' mm/s'],titleF{:}},'FontSize',8)
                    subplot(2,4,sub+5); plot(V,V,'k'); hold on; errorbar(mean_V',mean_measV(:,sL,ex),std_measV(:,sL,ex),'.'); xlabel('Known velocity (mm/s)'); ylabel('Measured velocity (mm/s)');
                    pos=get(gca,'Position'); set(gca,'Position',[pos(1) pos(2) pos(3) pos(4)*0.8])
                    axes('Position',get(gca,'Position'),'YTick',[],'XAxisLocation','top','Color','none','box','on','XLim',get(gca,'XLim')./rateV_factor); xlabel(xlabelS); title({'Mean velocities',['Searching ',num2str(searchLims(sL)*100),'% xcorr; Excluding |V| > ',num2str(extremes),' mm/s'],titleF{:}},'FontSize',8)
                    
                    
                    meanF_all(sL,ex) = mean(abs((V(abs(V)>0 & abs(V)<vlim & abs(measV(:,sL))<extremes)-measV(abs(V)>0 & abs(V)<vlim & abs(measV(:,sL))<extremes,sL))./V(abs(V)>0 & abs(V)<vlim & abs(measV(:,sL))<extremes)));
                    stdF_all(sL,ex) = std(abs((V(abs(V)>0 & abs(V)<vlim & abs(measV(:,sL))<extremes)-measV(abs(V)>0 & abs(V)<vlim & abs(measV(:,sL))<extremes,sL))./V(abs(V)>0 & abs(V)<vlim & abs(measV(:,sL))<extremes)));
                    try
                        meanF_means(sL,ex) = mean(abs((mean_V(abs(mean_V)>0 & abs(mean_V)<vlim & abs(mean_measV(:,sL,ex))'<extremes)'-mean_measV(abs(mean_V)>0 & abs(mean_V)<vlim & abs(mean_measV(:,sL,ex))'<extremes,sL,ex))./mean_V(abs(mean_V)>0 & abs(mean_V)<vlim & abs(mean_measV(:,sL,ex))'<extremes)'));
                        stdF_means(sL,ex) = std(abs((mean_V(abs(mean_V)>0 & abs(mean_V)<vlim & abs(mean_measV(:,sL,ex))'<extremes)'-mean_measV(abs(mean_V)>0 & abs(mean_V)<vlim & abs(mean_measV(:,sL,ex))'<extremes,sL,ex))./mean_V(abs(mean_V)>0 & abs(mean_V)<vlim & abs(mean_measV(:,sL,ex))'<extremes)'));
                    catch %if only zero velocities, will get an error
                        warning = 'unable to calculate fractional error of mean measurements'
                        meanF_means(sL,ex)=NaN;
                        stdF_means(sL,ex)=NaN;
                    end
                    meanR_means(sL,ex) = mean(std_measV(abs(mean_V)>0 & abs(mean_V)<vlim  & abs(mean_measV(:,sL,ex))'<extremes,sL,ex));
                    stdR_means(sL,ex) = std(std_measV(abs(mean_V)>0 & abs(mean_V)<vlim  & abs(mean_measV(:,sL,ex))'<extremes,sL,ex));
                    
                end
                if savePlots
                    saveas(gcf,[datestr(now,'ddmmmyy_HH-MM-SS-AM'),'_',fileF,',excludingV-',num2str(extremes),'.fig'])
                end
                ex = ex+1;
            end
            
            
            if plotAnalysis
                
                figure('Position',[0 0 scrsz(3) scrsz(4)/2]);
                subplot(131); bar(meanF_all(:,:)); hold on
                errorbar(0.85:1:length(meanF_all),meanF_all(:,1),stdF_all(:,1),'.');
                errorbar(1.15:1:length(meanF_all)+1,meanF_all(:,2),stdF_all(:,2),'r.');
                title({titleF{:},'Mean fractional error of all data',['Searching ',num2str(searchLims(1)*100),'% xcorr (1) and ',num2str(searchLims(2)*100),'% xcorr (2)']})
                legend(['Excluding |V| > ',num2str(extreme(1)),' mm/s'],['Excluding |V| > ',num2str(extreme(2)),' mm/s'])
                ylabel('Fractional error')
                pos=get(gca,'Position'); set(gca,'Position',[pos(1) pos(2) pos(3) pos(4)*0.8])
                
                subplot(132); bar(meanF_means(:,:)); hold on
                errorbar(0.85:1:length(meanF_means),meanF_means(:,1),stdF_means(:,1),'.');
                errorbar(1.15:1:length(meanF_means)+1,meanF_means(:,2),stdF_means(:,2),'r.');
                title({titleF{:},'Mean fractional error of data means',['Searching ',num2str(searchLims(1)*100),'% xcorr (1) and ',num2str(searchLims(2)*100),'% xcorr (2)']})
                legend(['Excluding |V| > ',num2str(extreme(1)),' mm/s'],['Excluding |V| > ',num2str(extreme(2)),' mm/s'])
                ylabel('Fractional error')
                pos=get(gca,'Position'); set(gca,'Position',[pos(1) pos(2) pos(3) pos(4)*0.8])
                
                subplot(133); bar(meanR_means(:,:)); hold on
                errorbar(0.85:1:length(meanR_means),meanR_means(:,1),stdR_means(:,1),'.');
                errorbar(1.15:1:length(meanR_means)+1,meanR_means(:,2),stdR_means(:,2),'r.');
                title({titleF{:},'Mean data resolution (std of means)',['Searching ',num2str(searchLims(1)*100),'% xcorr (1) and ',num2str(searchLims(2)*100),'% xcorr (2)']})
                legend(['Excluding |V| > ',num2str(extreme(1)),' mm/s'],['Excluding |V| > ',num2str(extreme(2)),' mm/s'])
                ylabel('Resolution (mm/s)')
                pos=get(gca,'Position'); set(gca,'Position',[pos(1) pos(2) pos(3) pos(4)*0.8])
                
                if savePlots
                    saveas(gcf,[datestr(now,'ddmmmyy_HH-MM-SS-AM'),'_',fileF,'_analysis.fig'])
                end
                
                dlmwrite([datestr(now,'ddmmmyy_HH-MM-SS-AM'),'_',fileF,'.csv'],[rates rep V./factor ts(:,1) ts(:,2) V measV(:,1) measV(:,2)],'delimiter','\t','newline','pc','precision',8)
                for ex=1:2
                    dlmwrite([datestr(now,'ddmmmyy_HH-MM-SS-AM'),'_mean-',fileF,',excludingV-',num2str(extreme(ex)),'.csv'],[r' mean_V'./factor ...
                        meants(:,1,ex) stdts(:,1,ex) meants(:,2,ex) stdts(:,2,ex) mean_V' mean_measV(:,1,ex) std_measV(:,1,ex) mean_measV(:,2,ex) std_measV(:,2,ex)],'delimiter','\t','newline','pc','precision',8)
                end
                
            end
            
        end
        
        
    end
    
    
end

