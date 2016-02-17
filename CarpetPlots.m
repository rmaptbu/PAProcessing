clear;
% path for PA data
basepath=...
'/Users/Thore/Documents/MATLAB/PAProcessing/Speckle/0205_Slab02OnGlass_NearPerpendicular/';
cd(basepath)

names={};
for i=1:10
names=[names,['rate_0_66.6sep_250MHz_16.300trig_1files_',num2str(i),'.csv']];
end


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
    %take 1000 point, split them off, filter attach to "pressure"
    %step by 1000, repeat, until end of file
    %split all data inte separate acquistions (indext in dimension 2)
    for x=1000:1000:size(all_points,1);
        i=i+1;
        
        %Split file
        dd=all_points(x-999:x,2);
        
        %Filter
%         order=2;
%         lowF=250;
%         sampling_rate = 4000;
%         [b,a]=butter(order,lowF/sampling_rate*2,'low');
%         d1 = filter(b,a,dd);
%         d1r = wrev(d1);
%         d2 = filter(b,a,d1r);
%         dd = wrev(d2);
        
        pressure(u,1:1000,i)= dd;
    end 
    clear('all_points','dd','d1','d1r','d1','b','a');
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