classdef SignalGenerator < handle
    properties
        %oscilloscope settings
        dt
        points
        tmax
        %experimental setup
        num_sig %number of succesive signals in one acquisition
        num_acq %total number fo signals (e.g. pairs if num_sig=2)
        delay %known time delay
        %Data
        time
        pressure
        voltage
    end
    methods
        function obj = SignalGenerator(points,dt,num_sig,num_acq)
            obj.points = points; %Number of points per acquisition
            obj.dt = dt; %Sampling intervall in ns
            obj.tmax = dt*points;
            obj.num_sig = num_sig; %number of succesive signals
            obj.num_acq = num_acq; %total number of signal(pairs,tripples,)
            obj.time = 0:obj.dt:obj.tmax-obj.dt;
        end
        function GeneratePressure(obj,delay)
            obj.delay=delay;
            obj.pressure = [];
            obj.voltage = [];            
            ret_time = obj.time'; %Retarded time
            for i=1:obj.num_sig-1;
                ret_time=cat(2,ret_time,obj.time'+i*delay);
            end
            h = waitbar(0, 'Initialising Waitbar');
            for j=1:obj.num_acq 
                msg='Generating Pressure...';
                waitbar(j/obj.num_acq,h,msg);
                p=zeros(size(ret_time));
                for i=0:5
                    % gaussian envelope
                    sig = 150;
                    c = obj.tmax/2;
                    % ultrasound
                    f = 10*(1+randn) ; %frequency in Mhz
                    omega = 2*pi*f*1E-3; %angular frequency
                    %envelope
                    env = gaussmf(ret_time, [sig c]);
                    % generate shape
                    p = p + 1/(100)*env.*exp(1i*ret_time*omega+2*pi*rand);
                end
                
                % measure pressure
                env = gaussmf(repmat(obj.time,obj.num_sig,1), [50 obj.tmax/2]);
                snr = 10.0; %in dB
                v = awgn(p.*env',snr,'measured');
                obj.pressure = cat(2,obj.pressure,p);
                obj.voltage = cat(2,obj.voltage,v);                
            end
            close(h);
        end
        function draw(obj,idx)
            hold on; box on;
            %add content
            axis([0, obj.tmax, -inf, inf]);
            for i=1:obj.num_sig               
                plot(obj.time,real(obj.voltage(:,(idx-1)*obj.num_sig+i)),...
                    '-','LineWidth', 2, 'Color',...
                    [i/obj.num_sig, .5-i/(obj.num_sig*2),1-i/obj.num_sig]);
            end
            %figure setup
            title(['Detected Pressures, delay = ', num2str(obj.delay)]);
            %legend('Absolute', 'Real');
            set(gca, 'FontName', 'Times New Roman');
            set(gca, 'FontSize', 12);
            xlabel('Time (ns)', 'FontSize', 16);
            ylabel('Pressure (a.u.)', 'FontSize', 16);
            hold off;
        end
    end
end