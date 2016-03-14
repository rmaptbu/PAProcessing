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
            for j=1:obj.num_acq                
                p=zeros(size(ret_time));
                for i=0:2
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
        end
    end
end