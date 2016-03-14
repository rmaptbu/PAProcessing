classdef PAF < handle %PAF is a handle class
    properties
        %Oscilloscope options
        dt
        sampling_rate
        %extrinsic data
        flow_rate
        %Xcorr of entire waveform options
        corrmin
        corrmax
        %Data
        pressure %pressure(data, #waveform)       
        shift %shift calculated from xcorrs: shift(flow_rate, shift(ns))
        profile %time gating profile: profile(shift(ns))
        imf %intrinsic mode functions
        %format: imf{#waveform}(mode, pressure)
    end
    methods
        function obj = PAF(dt)
        obj.dt=dt;        
        obj.sampling_rate=1.0/dt*1E3;
        end
        function ReadFlowrate(obj, filename)
            expr='-?[\d.]+_0\.5sep';
            [start, fin]=regexp(filename,expr);
            obj.flow_rate = str2double(filename(start:fin-7));
        end
        function ReadData(obj,filename, other_files)
            obj.pressure=[];
            %shift(flowrate, xcorr_peak_pos in ns)
            %profile is timegating xcorr peak pos plotted against windows position
            expr='-?[\d.]+_0\.5sep';
            [start, fin]=regexp(filename,expr);
            obj.flow_rate = str2double(filename(start:fin-7));
            
            %load all files into "all_points" (just some variable)
            all_points = csvread(filename);
            filename=filename(1:end-5);
            %add all other files with same rates to this dataset
            for i=1:length(other_files);
                all_points=cat(1,all_points,...
                    csvread([filename,num2str(other_files(i)),'.csv']));
            end
            %
            i=0;
            %take 5000 point, split them off, filter attach to "pressure"
            %step by 5000, repeat, until end of file
            %split all data inte separate acquistions (indext in dimension 2)
            for x=5000:5000:size(all_points,1);
                i=i+1;
                %Split file
                obj.pressure(1:5000,i)=all_points(x-4999:x,2);
            end
            obj.corrmin=1;
            obj.corrmax=5000;
            
        end
        function lowpass(obj,frequency)
            %Filter
            Fp=frequency; %enter in MHz
            Fp=Fp*2*pi/obj.sampling_rate; %sampling in samples per microscond
            %stop band frequency in rad/sample
            Fst=frequency*1.1; %in Mhz
            Fst=Fst*2*pi/obj.sampling_rate;
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
            for i=1:size(obj.pressure,2)
                obj.pressure(:,i) = filter(hd,obj.pressure(:,i));
                %             d1r = wrev(d1);
                %             d2 = filter(hd,d1r);
                %             pressure_lowp(:,i) = wrev(d2);
            end
        end
        function highpass(obj,frequency)
            %pass band frequency in rad/sample
            Fp=frequency; %enter in MHz
            Fp=Fp*2*pi/obj.sampling_rate; %sampling in samples per microscond
            %stop band frequency in rad/sample
            Fst=frequency/2; %in Mhz
            Fst=Fst*2*pi/obj.sampling_rate;
            d = fdesign.highpass('Fst,Fp,Ast,Ap',Fst,Fp,60,1);
            hd = design(d,'butter');
            for i=1:size(obj.pressure,2)
                obj.pressure(:,i) = filter(hd,obj.pressure(:,i));
                %             d1r = wrev(d1);
                %             d2 = filter(hd,d1r);
                %             pressure_lowp(:,i) = wrev(d2);
            end
        end
        function wallfilter(obj)
            meansignal=squeeze(mean(obj.pressure,2));
            obj.pressure=obj.pressure-repmat(meansignal,...
                [1 size(obj.pressure,2)]);
        end
        function EMD2Pressure(obj,IMFs,low,high)
            N=size(IMFs,2);
            obj.pressure=[];
            for i=1:N
                obj.pressure(:,i)=sum(IMFs{i}(1+high:end-low,:),1);
            end
        end
        function emd(obj,low,high) %empirical mode decomposition
            obj.imf={};
            h = waitbar(0, 'Initialising Waitbar');
            for i=1:size(obj.pressure,2)
                msg=['Emperical Mode Decomposition: ',...
                    num2str(i/size(obj.pressure,2)*100),'%'];
                waitbar(i/size(obj.pressure,2),h,msg);
                imf = emd(obj.pressure(:,i));
                p = cell2mat(imf');
                obj.pressure(:,i)=sum(p(1+high:end-low,:),1);
                obj.imf{i} = p;
            end
            close(h);
        end
        function shift = xcorr(obj, corrmin, corrmax)
            obj.corrmin = corrmin;
            obj.corrmax = corrmax;

            %cross correlate all signal and normalise by maximum
            width=obj.corrmax-obj.corrmin;
            for i=1:2:size(obj.pressure,2)-1;
                xcorrs(:,(i+1)/2)=xcorr(obj.pressure(obj.corrmin:obj.corrmax,i+1),...
                    obj.pressure(obj.corrmin:obj.corrmax,i),'biased');
                xcorrs_norm(:,(i+1)/2)=xcorrs(:,(i+1)/2)./max(xcorrs(:,(i+1)/2));
                %xcorrs_norm: normalised xcorr of all pairs
                %xcorrs_norm(xcorr_value, Pair_index)
            end
            
            % find maximum of x-correlation of entire waveform
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
            t = t * obj.dt;%map of time points
            %create 100+100 points around maximum data point
            %if I=1 then peak is at centre. i.e. need to shift by one
            %to get rid of offset (peak at I-1 in absolute time)
            tI = (shift-width-2):1.0/100:(shift-width);
            tI = tI * obj.dt;
            xcorr_interp = ... %100 point interpolation
                interp1(t,meanxcorr,tI,'spline');
            [~,max_pos]=max(xcorr_interp');
            xcorr_peak_pos=tI(max_pos);
            
            
            shift=(shift-width)*obj.dt;
            shift=[obj.flow_rate, shift];
            obj.shift=shift;
            %SHIFT_ALL=cat(1,SHIFT_ALL,xcorr_peak_pos);
        end
        function profile = TimeGating(obj, N, q, w, W, normalise, remove_outliers)
            jmax=W/q;
            assert(~rem(W,q),'W/q not integer');
            %starting point for xcorr
            corr_bounds=ceil(0.2*w);
            %between 0and1.defines size of search for corr peak
            %assume first pair correlates
            h = waitbar(0, 'Initialising Waitbar'); 
            xcorr_it_max=size(obj.pressure,2)-1;
            xcorrwindows=zeros(2*w+1,size(obj.pressure,2)/2,jmax);
            for i=1:2:size(obj.pressure,2)-1;
                msg=['Time Gating: ',num2str(i/xcorr_it_max*100),'%'];
                waitbar(i/xcorr_it_max,h,msg);
                for j=1:jmax                    
                    xcorrwindows(1:2*w+1,(i+1)/2,j)=xcorr(...
                        obj.pressure((N+(j-1)*q):(N+(j-1)*q+w),i+1),...
                        obj.pressure((N+(j-1)*q):(N+(j-1)*q+w),i),'biased');
                end
            end
            close(h);
            %Normalise
%             if normalise
%                 for i=1:length(xcorrwindows(1,:,1))
%                     for j=1:length(xcorrwindows(1,1,:))
%                         xcorrwindows(:,i,j)=...
%                             xcorrwindows(:,i,j)./max(xcorrwindows(:,i,j));
%                     end
%                 end
%             end
            
            xcorrwindows_ensemble=squeeze(mean(xcorrwindows,2));
            %w+1 is centre
            [M, I]=max(xcorrwindows_ensemble(w+1-corr_bounds:w+1+corr_bounds,:),[],1);
            %adjust I to work with indexing of xcorrwindows_ensemble
            %coor_bounds+1 is at centre. needs to be shifted to w+1
            I=I+w-corr_bounds;
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
            t=t*obj.dt;
            %map of time points
            for i=1:size(I,1)
                %create 100+100 points around maximum data point
                %if I=1 then peak is at centre. i.e. need to shift by one
                %to get rid of offset (peak at I-1 in absolute time)
                tI(i,:) = (I(i)-w-2):1.0/100:(I(i)-w);
                tI(i,:) = obj.dt*tI(i,:);
                xcorr_interp(i,:) = ... %100 point interpolation
                    interp1(t,xcorrwindows_ensemble(:,i),tI(i,:),'spline');
            end
            [~,max_pos]=max(xcorr_interp');
            for i=1:size(max_pos,2)
                xcorr_peak_pos(i)=tI(i,max_pos(i));
            end
            
            profile=xcorr_peak_pos';
            obj.profile = profile;
        end
        function profile = IMFTimeGating(obj, N, q, w, W, order)
            assert(~isempty(obj.imf), ...
                'Need to generate intrinsic mode functions first. Run obj.emd.');
            warning('Pressure Data is being overwritten');
%             pressure = obj.pressure; %store pressure as backup
            for i=1:size(obj.imf,2) %iterate of intrinsic mode functions 
            obj.pressure(:,i) = obj.imf{i}(order,:)';
            end
            profile = obj.TimeGating(N,q,w,W,0,0);
%             obj.pressure = pressure;
        end
        function draw(obj)
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
            meansignal=squeeze(mean(obj.pressure,2));
            plot(meansignal);
            hold on;
            for i=1:jmax-1;
                plot(N+i*q:N+i*q+49,meansignal(N+i*q:N+i*q+49),...
                    'LineWidth',2,'Color',[i/jmax,.5-i/(jmax*2),1-i/jmax]);
            end ;
            
            %all xcorrs
            subplot(4,2,3:4);
            figure;hold on
            for i=1:jmax-1;
                fig1=plot(squeeze(xcorrwindows_ensemble(:,i)./...
                    max(abs(xcorrwindows_ensemble(:,i)))),...
                    'Color',[.2,.2,.2]);
                fig1.Color(4)=0.1;
            end
            hold off
            %         figure;plot(maxprofile);
            %     savefig([filename,'.fig']);
            
            figure;
            subplot(1,2,1);hold on
            plot(obj.pressure(:,1));
            for i=1:jmax-1;
                plot(N+i*q:N+i*q+49,obj.pressure(N+i*q:N+i*q+49,1),...
                    'LineWidth',2,'Color',[i/jmax,.5-i/(jmax*2),1-i/jmax]);
            end ;hold off
            subplot(1,2,2);hold on
            plot(obj.pressure(:,2));
            for i=1:jmax-1;
                plot(N+i*q:N+i*q+49,obj.pressure(N+i*q:N+i*q+49,2),...
                    'LineWidth',2,'Color',[i/jmax,.5-i/(jmax*2),1-i/jmax]);
            end ;
        end
    end
end

