classdef PAF < handle %PAF is a handle class
    properties
        %Oscilloscope options
        dt
        points
        sampling_rate
        %extrinsic data
        flow_rate
        %Xcorr of entire waveform options
        corrmin
        corrmax
        %Data
        pressure %pressure(data, #waveform)
        xcorrwindows_ensemble
        shift %shift calculated from xcorrs: shift(flow_rate, shift(ns))
        profile %time gating profile: profile(shift(ns))
        imf %intrinsic mode functions
        %format: imf{#waveform}(mode, pressure)
        %Timegating option
        N
        q
        w
        W
        jmax
    end
    methods
        function obj = PAF(dt)
            obj.dt=dt;
            obj.sampling_rate=1.0/dt*1E3;
        end
        function ReadFlowrate(obj, input)
            if ischar(input)
                filename=input;
            expr='-?[\d.]+_0\.5sep';
            [start, fin]=regexp(filename,expr);
            obj.flow_rate = str2double(filename(start:fin-7));
            elseif isnumeric(input)
                obj.flow_rate = input
            else
                error('filetype(input) is is not accepted as input');
            end  
        end
        function ReadData(obj,file, other_files)
            if ischar(file)
                obj.pressure=[];
                %shift(flowrate, xcorr_peak_pos in ns)
                %profile is timegating xcorr peak pos plotted against windows position
                expr='-?[\d.]+_0\.5sep';
                [start, fin]=regexp(file,expr);
                obj.flow_rate = str2double(file(start:fin-7));
                
                %load all files into "all_points" (just some variable)
                all_points = csvread(file);
                file=file(1:end-5);
                %add all other files with same rates to this dataset
                for i=1:length(other_files);
                    all_points=cat(1,all_points,...
                        csvread([file,num2str(other_files(i)),'.csv']));
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
            elseif isreal(file)
                %remember to set flowrate explicitly (obj.flow_rate=..)
                obj.pressure=file;
            else
                error('filetype(file) is is not accepted as input');
            end            
            
            obj.corrmin=1;
            obj.corrmax=5000;
            obj.points=5000;
            
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
        function ReloadPressure(obj)
            obj.EMD2Pressure(obj.imf,0,0);
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
        function shift = xcorr41(obj, corrmin, corrmax)
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
            obj.N=N;
            obj.q=q;
            obj.w=w;
            obj.W=W;
            obj.jmax=jmax;
            assert(~rem(W,q),'W/q not integer');
            %starting point for xcorr
            corr_bounds=ceil(0.9*w);
            %between 0and1.defines size of search for corr peak
            %assume first pair correlates
            h = waitbar(0, 'Initialising Waitbar');
            xcorr_it_max=size(obj.pressure,2)-1;
            xcorrwindows=zeros(6*w+1,size(obj.pressure,2)/2,jmax);
            for i=1:2:size(obj.pressure,2)-1;
                msg=['Time Gating: ',num2str(i/xcorr_it_max*100),'%'];
                waitbar(i/xcorr_it_max,h,msg);
                for j=1:jmax                    
                    f1=obj.pressure((N+(j-1)*q)-w:(N+(j-1)*q+w)+w,i);
                    f2=[zeros(w,1);...
                        obj.pressure((N+(j-1)*q):(N+(j-1)*q+w),i+1);...
                        zeros(w,1)]; 
                    %Linear fit of f1
                    X=[ones(length(f1),1),[1:length(f1)]'];
                    b=f1'/X'; %b(1) intercept, b(2) gradient
                    f_lin=b(2)*X(:,2);
                    f1=f1 - f_lin;
                    f2(w+1:2*w+1)=f2(w+1:2*w+1) - f_lin(w+1:2*w+1);
                    xcorrwindows(1:6*w+1,(i+1)/2,j)=xcorr(f1,f2,'none');
                    if i==1 && (j==1 || j==6)
                        figure;
                        subplot(2,2,1);hold on;
                        plot(f1,'LineWidth',2,'Color',[0 0 0]);                        
                        plot(f2,'LineWidth',2,'Color',[0.5 0.5 0.5]);
                        plot(f_lin,'LineWidth',2,'Color',[1 0.5 0.0]);
                        xlim([0 size(f2,1)]);
                        hold off
                        subplot(2,2,3:4);
                        plot(xcorr(f1,f2,'none'),...
                            'LineWidth',2,'Color',...
                            [j/obj.jmax,.5-j/(obj.jmax*2),1-j/obj.jmax]);
                        xlim([0 size(f2,1)]*2);
                    end
                end
            end
            close(h);
            xcorrwindows=xcorrwindows(2*w+1:4*w+1,:,:);
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
            obj.xcorrwindows_ensemble = xcorrwindows_ensemble;
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
        function draw(obj,varargin)
            if isempty(varargin)
                figname=['flowrate = ',num2str(obj.flow_rate)];
                savefig = 0;
            else 
                figname = varargin{1};
                savefig = 1;
            end
            
            fig=figure('Position', [500, 500, 640, 790]);
            set(0,'DefaultAxesFontSize', 11)
            
            %profile
            subplot(3,1,3)   
            hold on; box on
            dt_w=obj.dt*obj.q;
            for i=1:obj.jmax-1;
                plot(dt_w*(i-1):dt_w:dt_w*i,obj.profile(i:i+1),'LineWidth',2,'Color',...
                    [i/obj.jmax,.5-i/(obj.jmax*2),1-i/obj.jmax]);
            end ;
            xlim([1 inf]);
            xlabel('Time (ns)');
            ylabel('Xcorr Shift (ns)')
            
            %overview plot           
            subplot(3,1,1)
            t=obj.dt:obj.dt:(obj.points*obj.dt);
            sfig{2}=plot(t,obj.pressure(:,2),'LineWidth',2,'Color',[0.5 0.5 0.5]);
            hold on;
            sfig{1}=plot(t,obj.pressure(:,1),'LineWidth',2,'Color','Black');
            for i=1:obj.jmax;
                plot(obj.dt*(obj.N+(i-1)*obj.q):obj.dt:obj.dt*(obj.N+(i-1)*obj.q+obj.w),...
                    obj.pressure(obj.N+(i-1)*obj.q:obj.N+(i-1)*obj.q+obj.w,1),...
                    'LineWidth',2,'Color',[i/obj.jmax,.5-i/(obj.jmax*2),1-i/obj.jmax]);
            end ;
            legend([sfig{:}],{'Signal 1', 'Signal 2'});
            xlim([t(1) t(end)]);
            xlabel('Time (ns)');
            ylabel('Pressure (a.u.)');
            title(figname,'Interpreter', 'none');
            %Individual xcorrs
            subplot(3,1,2)
            hold on; box on
            for i=1:obj.jmax-1;
                plot(-obj.w*obj.dt:obj.dt:obj.w*obj.dt,...
                    obj.xcorrwindows_ensemble(:,i)./max(obj.xcorrwindows_ensemble(:,i)),...
                    'LineWidth',2,'Color',...
                    [i/obj.jmax,.5-i/(obj.jmax*2),1-i/obj.jmax]);
            end ; 
            xlabel('Xcorr Shift (ns)');
            ylabel('Corraltion (a.u.)');
            
            %Save figure
            if savefig
            set(gcf,'PaperPositionMode','auto')
            print(fig,figname,'-dpng','-r0')
            close(fig);
            end
        end
    end
end

