clear;
% path for PA data
basepath=...
'/Users/Thore/Documents/PIVLabData/oct10_10muSpheres_400muTube/PA/';  
% basepath=...
% '/Users/Thore/Documents/PIVLabData/07.19_10muSpheresTween_400muTube_3.62kpa/';  
cd(basepath)
cases=struct2cell(dir);  % get a list of the file names
cases=cases(1,3:end);

% names: it is important that the same number of symbols are in every
% filename. add " "(space) after filename to make sure.
names=cellstr(['rate_-0.1_4.0sep_full_16.300trig_10files_3.csv  ';
    'rate_-0.1_8.0sep_full_16.300trig_10files_0.csv  ';
    'rate_-0.1_12.0sep_full_16.300trig_10files_1.csv ';
    'rate_-0.1_16.0sep_full_16.300trig_10files_3.csv ';
    'rate_-0.1_20.0sep_full_16.300trig_10files_0.csv ';
    'rate_-0.1_40.0sep_full_16.300trig_10files_0.csv ';   
    'rate_-0.1_60.0sep_full_16.300trig_10files_0.csv ';
    'rate_-0.1_80.0sep_full_16.300trig_10files_0.csv ']);
% 
% names=cellstr(['rate_-0.20_2.0sep_full_16.300trig_1files_5.csv  ';
%     'rate_-0.20_4.0sep_full_16.300trig_1files_4.csv  ';
%     'rate_-0.20_6.0sep_full_16.300trig_1files_3.csv  ';
%     'rate_-0.20_8.0sep_full_16.300trig_1files_2.csv  ';
%     'rate_-0.20_10.0sep_full_16.300trig_1files_1.csv ';
%     'rate_-0.20_15.0sep_full_16.300trig_1files_6.csv ';
%     'rate_-0.20_20.0sep_full_16.300trig_1files_7.csv ';
%     'rate_-0.20_30.0sep_full_16.300trig_1files_9.csv ';
%     'rate_-0.20_40.0sep_full_16.300trig_1files_10.csv']);

% %Read files that were accepted
% fileID = fopen('velocitymeasurements_full_peak0_Interpolant.csv','r');
% formatSpec = '%f %f %f %f %f %f %f %f %f %f %f %f Good result';
% % 1   Dial/rate (ml/h)
% % 2   Record number
% % 3   Number of files
% % 4   Known velocity (mm/s)
% % 5   Known resolution (mm/s)
% % 6   Measured velocity from fit to xcorr peak (mm/s)
% % 7   Measured resolution from fit to xcorr peak (mm/s)
% % 8   Measured pulse separation (ms)
% % 9   Theoretical resolution based on oscilloscope sampling interval (mm/s)
% % 10  Amplitude of selected peak
% % 11  Amplitude of selected peak relative to the next largest peak
% % 12  Amplitude of selected peak relative to the RMS of the xcorr amplitude
% % 13  Quality of result (options: "Good result", Default; "Bad result 
% % (incorrect pulse sep - laser misfiring)"; "Bad 
% size = [12 Inf];
% files_raw=fscanf(fileID,formatSpec, size);
% fclose(fileID);
%Convert files_raw into filenames

% names=cellstr('rate_+15_0.5sep_250MHz_14.012trig_5files_')
maxprofile=[];
maxprofile_interp=[];
shift_all=[];
draw=0; %only works with time_gating=1
remove_outliers=1;
time_gating=1;
dt=0.25; %sampling interval in ns
number_of_files=length(names);

for u=1:number_of_files;
    display(u);
    filename=names{u};
    
    %load all files into "all_points" (just some variable)
    all_points = csvread([filename,'']); 
    %whatever it says here needs to add up to the full name of the file
%     all_points=cat(1,all_points,csvread([filename,'2.csv']));
    %
    i=0;
    %take 5000 point, split them off, filter attach to "pressure"
    %step by 5000, repeat, until end of file
    %split all data inte separate acquistions (indext in dimension 2)
    for x=5000:5000:size(all_points,1);
        i=i+1;
        
        %Split file
        dd=all_points(x-4999:x,2);
        
        %Filter
        order=2;
        lowF=250;
        sampling_rate = 4000;
        [b,a]=butter(order,lowF/sampling_rate*2,'low');
        d1 = filter(b,a,dd);
        d1r = wrev(d1);
        d2 = filter(b,a,d1r);
        dd = wrev(d2);
        
        pressure(1:5000,i)= dd;
    end 
    clear('all_points','dd','d1','d1r','d1','b','a');
       
    %% XCORR
    meansignal=squeeze(mean(pressure,2));
    %cross correlate all signal and normalise by maximum
    for i=1:2:size(pressure,2)-1;
        xcorrs(:,(i+1)/2)=xcorr(pressure(:,i),pressure(:,i+1),'unbiased');
        xcorrs_norm(:,(i+1)/2)=xcorrs(:,(i+1)/2)./max(xcorrs(:,(i+1)/2));
        %xcorrs_norm: normalised xcorr of all pairs
        %xcorrs_norm(xcorr_value, Pair_index)        
    end
    
    
    %% find maximum of x-correlation of entire waveform
    %take average of all normalises cross correlations and find position of
    %maximum
    %"shift" is position of maximum
    [~, shift_i]=max(xcorrs_norm);
    shift_std=std(shift_i); clear('shift_i');
    meanxcorr=squeeze(mean(xcorrs_norm,2));
    [~, shift]=max(meanxcorr);
    
    %Interpolate location of maximum
    t=-(5000-1)*dt:dt:(5000-1)*dt; %map of time points
    %create 100+100 points around maximum data point
    %if I=1 then peak is at centre. i.e. need to shift by one
    %to get rid of offset (peak at I-1 in absolute time)
    tI = dt*(shift-5002):dt/100:dt*(shift-5000);
    xcorr_interp = ... %100 point interpolation
        interp1(t,meanxcorr,tI,'spline');
    [~,max_pos]=max(xcorr_interp');
    xcorr_peak_pos=tI(max_pos);

        
    shift=(shift-5000)*dt;
    shift_all=cat(1,shift_all,xcorr_peak_pos);
    
    clear('max_shift_i','meanxcorr','shift', 'xcorrs','xcorrs_norm');
    %% Time gating
    if time_gating    
        %starting point for xcorr
        N=1000;%starting poing
        q=50;%stepsize
        w=250; %interrogation windows
        W=2000; %walking length
        jmax=ceil(W/q); %number of plots
        corr_bounds=0.6; %between 0and1.defines size of search for corr peak 
        %assume first pair correlates
        for i=1:2:size(pressure,2)-1;
            for j=1:jmax
                xcorrwindows(1:2*w-1,(i+1)/2,j)=xcorr(...
                    pressure((N+j*q):(N+j*q+w-1),i),...
                    pressure((N+j*q):(N+j*q+w-1),i+1),'unbiased');
            end
        end
        
        for i=1:length(xcorrwindows(1,:,1))
            for j=1:length(xcorrwindows(1,1,:))
                xcorrwindows_norm(:,i,j)=...
                    xcorrwindows(:,i,j)./max(xcorrwindows(:,i,j));
            end
        end
        
        xcorrwindows_ensemble=squeeze(mean(xcorrwindows,2));
        [M, I]=max(xcorrwindows_ensemble(w:w+ceil(w*0.7),:),[],1);
        %     [M, I]=max(xcorrwindows_ensemble(w-ceil(w*0.7):w,:),[],1);
        
        if remove_outliers
            for i=1:length(I);
                if I(i)<1 || I(i)>w*corr_bounds %peak out of bounds
                    %play with 1 and 0.6 to adjust bound for outliers
                    if i~=1
                        %replace outlier by previous data point
                        I_0(i)=I_0(i-1);                        
                    else
                        I_0(i)=0;
                    end
                    display(['changing ',num2str(i)])
                else
                    I_0(i)=I(i);
                end
            end
            I=I_0;
        end
        maxoffset=w;
        I=I';
        
        %Interpolate location of maximum
        t=-(w-1)*dt:dt:(w-1)*dt; %map of time points        
        for i=1:size(I,1)
            %create 100+100 points around maximum data point
            %if I=1 then peak is at centre. i.e. need to shift by one
            %to get rid of offset (peak at I-1 in absolute time)
            tI(i,:) = dt*(I(i)-2):dt/100:dt*(I(i));
            xcorr_interp(i,:) = ... %100 point interpolation
                interp1(t,xcorrwindows_ensemble(:,i),tI(i,:),'spline');            
        end
        [~,max_pos]=max(xcorr_interp');
        for i=1:size(max_pos,2) 
            xcorr_peak_pos(i)=tI(i,max_pos(i));
        end
        
        maxprofile_interp=cat(2,maxprofile_interp,xcorr_peak_pos');
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
                plot(I(i:i+1)+maxoffset,M(i:i+1),'b-o',...
                    'LineWidth',2,'Color',[i/jmax,.5-i/(jmax*2),1-i/jmax]);
            end
            hold off
        end
        
       
        %Mean signal
        subplot(4,2,1:2)
        plot(meansignal);
        hold on;
        for i=1:jmax-1;
            plot(N+i*q:N+i*q+49,meansignal(N+i*q:N+i*q+49),...
                'LineWidth',2,'Color',[i/jmax,.5-i/(jmax*2),1-i/jmax]);
        end ;
        
        %all xcorrs
        subplot(4,2,3:4);hold on
        for i=1:jmax-1;
            plot(squeeze(xcorrwindows_ensemble(:,i)./...
                max(abs(xcorrwindows_ensemble(:,i)))),...
                'Color',[i/jmax,.5-i/(jmax*2),1-i/jmax]);
        end
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
end
   %adjust the number in for loop to equal number of files anaylsed (n)
   %make sure subplot has enough space: subplot(X,Y,i)-> X*Y>n
   %make sure that you increase the subplot no as the no of files increase
  
   
  
 %shift_all contains the measured shift obtained by correlating the entire
 %waveform
 shift_all=padarray(shift_all, [0 size(maxprofile_interp,1)-1], ...
     'symmetric', 'post');
 %% create plot
 figure;
 y_pos=0;
 x_pos=0;
 ncol=4;
 nrow=3;
 plotsizeX=0.6/ncol;
 plotsizeY=0.5/nrow;
 for i=1:number_of_files;     
     if ~rem(i-1,ncol);y_pos=y_pos+1;x_pos=0; end;
     ax1=axes('Position',...
         [0.1+x_pos*(plotsizeX+0.045) 1-(plotsizeY+0.1)*y_pos plotsizeX plotsizeY],...
         'XTickLabel','',...
         'YTickLabel','shift (ns)');
     ax2=axes('Position',ax1.Position);
     plot(ax1,maxprofile_interp(:,i), 'Color', [0, 0, 0]);hold on;
     plot(ax2,shift_all(i,:)', 'Color', [1, 0, 0]);
     set(gca, 'Color', 'None','XColor','r','YColor','r',...
         'XAxisLocation','top',...
         'YAxisLocation','right',...
         'Box','off',...
         'XLim', ax1.XLim,'YLim',ax1.YLim,... 'XTick',[],'YTick',[],...
         'XColor', [0 0 0], 'YColor', [0 0 0]);
     S=strrep(names(i),'_','\_');
     S=S{1};
     S=[S(1:20),'\ xcorr', num2str(shift_all(i,1))];
     title(S);
     x_pos=x_pos+1;
 end
 i=i+1;
 if ~rem(i-1,ncol);y_pos=y_pos+1;x_pos=0; end;
 ax1=axes('Position',...
         [0.1+x_pos*(plotsizeX+0.1) 1-(plotsizeY+0.1)*y_pos plotsizeX*2 plotsizeY],...
         'XTickLabel','',...
         'YTickLabel','shift (ns)');
 plot(pressure(:,1),'Color',[0 0 0]);hold on;
 for i=1:jmax-1;
     plot(ax1,N+i*q:N+i*q+49,pressure(N+i*q:N+i*q+49,1),...
         'LineWidth',2,'Color',[i/jmax,.5-i/(jmax*2),1-i/jmax]);
 end;hold off;