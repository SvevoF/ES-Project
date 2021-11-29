
%% 0001
% 100 µF 5.5V da carico solo corrente con 1mV = 1mA. Misuriamo la corrente del consensatore. Ping ricevuto.
clc 
clearvars
close all
FigID=0;
% importing clean data measired with power supply 
filePattern = fullfile('..\..\ping\boot_ping_total','20210601-0001.csv');% name pattern for files
theFiles = dir(filePattern);% properties of the given name pattern
allFileNames = {theFiles.name};% list of all the files with the given name pattern
folder=theFiles.folder;
T1=readmatrix(fullfile(folder,allFileNames{1}));

% loading data
filePattern = fullfile('20211119-0001','*20211119-0001_*.csv');% name pattern for files
theFiles = dir(filePattern);% properties of the given name pattern
allFileNames = {theFiles.name};% list of all the files with the given name pattern
folder=theFiles.folder;
T=[];
for i=1:length(allFileNames)
    T=[T;readmatrix(fullfile(folder,allFileNames{i}))]; 
end
T1=T1(~isnan(T1(:,2)),:);
t1=T1(:,1)/1000;% [s]
i1=T1(:,2);% [mA]
% sampling frequency
Ts = T(2,1) - T(1,1);
Fs = 1/Ts;
T=T(~isnan(T(:,2)),:);% eliminating NaN values(problem with filtering)
% multiple measures and each measure start from zero,
% when we concatenate the we have to recreate the time scale:
T(:,1)=(0:1:length(T(:,1))-1)*Ts;
% extraction of variables
t=T(:,1);%[s]
i=-T(:,2).*1000;%[mA]
v=-T(:,3);% [v]
% filtering
i_f=lowpass(i,500,Fs,'Steepness',0.95);
i_f=smooth(i,100);
% aligning signals
test=i_f>5;
test1=i1>5;
or=find(test);
or1=find(test1,1);
t=t-t(or(100));
t1=t1-t1(or1);
% cutting to match the length of the other signal
i=i(1:find(t>t1(end),1));
i_f=i_f(1:find(t>t1(end),1));
v=v(1:find(t>t1(end),1));
t=t(1:find(t>t1(end),1));

% psd
[psd,f]=pwelch(i,[],[],[],Fs);
[psd_f,f_f]=pwelch(i_f,[],[],[],Fs);
% plot psd
FigID=FigID+1;
figure(FigID)
    hold on
    plot(f,10*log10(psd),'LineWidth',1)
    plot(f_f,10*log10(psd_f),'LineWidth',1)
    xlabel('Frequency (Hz)')
    ylabel('current PSD (dB/Hz)')
    legend('Unfiltered','Filtered')

% plot signal current
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,i,'LineWidth',1)
    xlabel('time')
    ylabel('current  [mA]')
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,i_f,'LineWidth',1)    
    plot(t1,i1,'LineWidth',1)  
    xlabel('time [s]')
    ylabel('current(filtered) [mA]')
% plot signal tension
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,v,'LineWidth',1)
    xlabel('time [s]')
    ylabel('tension [V]')
% energy computation
P= v.*i_f;% [mW]
E=trapz(t,P);% [mJ]
% plot power
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,P,'LineWidth',1)
    xlabel('time [s]')
    ylabel('power [mW]')

%% 0002
% uguale a 100 µF
clc 
clearvars
close all
FigID=0;
% importing clean data measired with power supply 
filePattern = fullfile('..\..\ping\boot_ping_total','20210601-0001.csv');% name pattern for files
theFiles = dir(filePattern);% properties of the given name pattern
allFileNames = {theFiles.name};% list of all the files with the given name pattern
folder=theFiles.folder;
T1=readmatrix(fullfile(folder,allFileNames{1}));

% loading data
filePattern = fullfile('20211119-0002.csv');% name pattern for files
theFiles = dir(filePattern);% properties of the given name pattern
allFileNames = {theFiles.name};% list of all the files with the given name pattern
folder=theFiles.folder;
T=[];
for i=1:length(allFileNames)
    T=[T;readmatrix(fullfile(folder,allFileNames{i}))]; 
end
T1=T1(~isnan(T1(:,2)),:);
t1=T1(:,1)/1000;% [s]
i1=T1(:,2);% [mA]
% sampling frequency
Ts = T(2,1) - T(1,1);
Fs = 1/Ts;
T=T(~isnan(T(:,2)),:);% eliminating NaN values(problem with filtering)
% multiple measures and each measure start from zero,
% when we concatenate the we have to recreate the time scale:
T(:,1)=(0:1:length(T(:,1))-1)*Ts;
% extraction of variables
t=T(:,1);%[s]
i=-T(:,2).*1000;%[mA]
v=-T(:,3);% [v]
% filtering
i_f=lowpass(i,500,Fs,'Steepness',0.95);
i_f=smooth(i,100);
% aligning signals
test=i_f>5;
test1=i1>5;
or=find(test);
or1=find(test1,1);
t=t-t(or(100));
t1=t1-t1(or1);
% cutting to match the length of the other signal
i=i(1:find(t>t1(end),1));
i_f=i_f(1:find(t>t1(end),1));
v=v(1:find(t>t1(end),1));
t=t(1:find(t>t1(end),1));

% psd
[psd,f]=pwelch(i,[],[],[],Fs);
[psd_f,f_f]=pwelch(i_f,[],[],[],Fs);
% plot psd
FigID=FigID+1;
figure(FigID)
    hold on
    plot(f,10*log10(psd),'LineWidth',1)
    plot(f_f,10*log10(psd_f),'LineWidth',1)
    xlabel('Frequency (Hz)')
    ylabel('current PSD (dB/Hz)')
    legend('Unfiltered','Filtered')

% plot signal current
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,i,'LineWidth',1)
    xlabel('time')
    ylabel('current  [mA]')
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,i_f,'LineWidth',1)    
    plot(t1,i1,'LineWidth',1)  
    xlabel('time [s]')
    ylabel('current(filtered) [mA]')
% plot signal tension
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,v,'LineWidth',1)
    xlabel('time [s]')
    ylabel('tension [V]')
% energy computation
P= v.*i_f;% [mW]
E=trapz(t,P);% [mJ]
% plot power
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,P,'LineWidth',1)
    xlabel('time [s]')
    ylabel('power [mW]')
%% 0003
% corrente come 0001 + Tensione SuperCap
clc 
clearvars
close all
FigID=0;
% importing clean data measired with power supply 
filePattern = fullfile('..\..\ping\boot_ping_total','20210601-0001.csv');% name pattern for files
theFiles = dir(filePattern);% properties of the given name pattern
allFileNames = {theFiles.name};% list of all the files with the given name pattern
folder=theFiles.folder;
T1=readmatrix(fullfile(folder,allFileNames{1}));

% loading data
filePattern = fullfile('20211119-0003','*20211119-0003_*.csv');% name pattern for files
theFiles = dir(filePattern);% properties of the given name pattern
allFileNames = {theFiles.name};% list of all the files with the given name pattern
folder=theFiles.folder;
T=[];
for i=1:length(allFileNames)
    T=[T;readmatrix(fullfile(folder,allFileNames{i}))]; 
end
T1=T1(~isnan(T1(:,2)),:);
t1=T1(:,1)/1000;% [s]
i1=T1(:,2);% [mA]
% sampling frequency
Ts = T(2,1) - T(1,1);
Fs = 1/Ts;
T=T(~isnan(T(:,2)),:);% eliminating NaN values(problem with filtering)
% multiple measures and each measure start from zero,
% when we concatenate the we have to recreate the time scale:
T(:,1)=(0:1:length(T(:,1))-1)*Ts;
% extraction of variables
t=T(:,1);%[s]
i=-T(:,2).*1000;%[mA]
v=-T(:,3);% [v]
% filtering
i_f=lowpass(i,500,Fs,'Steepness',0.95);
i_f=smooth(i,100);
% aligning signals
test=i_f>5;
test1=i1>5;
or=find(test,1);
or1=find(test1,1);
t=t-t(or);
t1=t1-t1(or1);
% cutting to match the length of the other signal
i=i(1:find(t>t1(end),1));
i_f=i_f(1:find(t>t1(end),1));
v=v(1:find(t>t1(end),1));
t=t(1:find(t>t1(end),1));

% psd
[psd,f]=pwelch(i,[],[],[],Fs);
[psd_f,f_f]=pwelch(i_f,[],[],[],Fs);
% plot psd
FigID=FigID+1;
figure(FigID)
    hold on
    plot(f,10*log10(psd),'LineWidth',1)
    plot(f_f,10*log10(psd_f),'LineWidth',1)
    xlabel('Frequency (Hz)')
    ylabel('current PSD (dB/Hz)')
    legend('Unfiltered','Filtered')

% plot signal current
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,i,'LineWidth',1)
    xlabel('time')
    ylabel('current  [mA]')
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,i_f,'LineWidth',1)    
    plot(t1,i1,'LineWidth',1)  
    xlabel('time [s]')
    ylabel('current(filtered) [mA]')
% plot signal tension
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,v,'LineWidth',1)
    xlabel('time [s]')
    ylabel('tension [V]')
% energy computation
P= v.*i_f;% [mW]
E=trapz(t,P);% [mJ]
% plot power
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,P,'LineWidth',1)
    xlabel('time [s]')
    ylabel('power [mW]')
%% 0004
% stesso file di 0003 salvato due volte.

%% 0005
% stesso setup 0003 senza led (SB18 traccia interrotta)
clc 
clearvars
close all
FigID=0;
% importing clean data measired with power supply 
filePattern = fullfile('..\..\ping\boot_ping_total','20210601-0001.csv');% name pattern for files
theFiles = dir(filePattern);% properties of the given name pattern
allFileNames = {theFiles.name};% list of all the files with the given name pattern
folder=theFiles.folder;
T1=readmatrix(fullfile(folder,allFileNames{1}));

% loading data
filePattern = fullfile('20211119-0005','*20211119-0005_*.csv');% name pattern for files
theFiles = dir(filePattern);% properties of the given name pattern
allFileNames = {theFiles.name};% list of all the files with the given name pattern
folder=theFiles.folder;
T=[];
for i=1:length(allFileNames)
    T=[T;readmatrix(fullfile(folder,allFileNames{i}))]; 
end
T1=T1(~isnan(T1(:,2)),:);
t1=T1(:,1)/1000;% [s]
i1=T1(:,2);% [mA]
% sampling frequency
Ts = T(2,1) - T(1,1);
Fs = 1/Ts;
T=T(~isnan(T(:,2)),:);% eliminating NaN values(problem with filtering)
% multiple measures and each measure start from zero,
% when we concatenate the we have to recreate the time scale:
T(:,1)=(0:1:length(T(:,1))-1)*Ts;
% extraction of variables
t=T(:,1);%[s]
i=-T(:,2).*1000;%[mA]
v=-T(:,3);% [v]
% filtering
i_f=lowpass(i,500,Fs,'Steepness',0.95);
i_f=smooth(i,100);
% aligning signals
test=i_f>5;
test1=i1>5;
or=find(test,1);
or1=find(test1,1);
t=t-t(or);
t1=t1-t1(or1);
% cutting to match the length of the other signal
i=i(1:find(t>t1(end),1));
i_f=i_f(1:find(t>t1(end),1));
v=v(1:find(t>t1(end),1));
t=t(1:find(t>t1(end),1));

% psd
[psd,f]=pwelch(i,[],[],[],Fs);
[psd_f,f_f]=pwelch(i_f,[],[],[],Fs);
% plot psd
FigID=FigID+1;
figure(FigID)
    hold on
    plot(f,10*log10(psd),'LineWidth',1)
    plot(f_f,10*log10(psd_f),'LineWidth',1)
    xlabel('Frequency (Hz)')
    ylabel('current PSD (dB/Hz)')
    legend('Unfiltered','Filtered')

% plot signal current
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,i,'LineWidth',1)
    xlabel('time')
    ylabel('current  [mA]')
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,i_f,'LineWidth',1)    
    plot(t1,i1,'LineWidth',1)  
    xlabel('time [s]')
    ylabel('current(filtered) [mA]')
% plot signal tension
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,v,'LineWidth',1)
    xlabel('time [s]')
    ylabel('tension [V]')
% energy computation
P= v.*i_f;% [mW]
E=trapz(t,P);% [mJ]
% plot power
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,P,'LineWidth',1)
    xlabel('time [s]')
    ylabel('power [mW]')
%% 0006
% 30 mF 5.5V carico, corrente + tensione parte superC, senza led. Ping ricevuto
clc 
clearvars
close all
FigID=0;
% importing clean data measired with power supply 
filePattern = fullfile('..\..\ping\boot_ping_total','20210601-0001.csv');% name pattern for files
theFiles = dir(filePattern);% properties of the given name pattern
allFileNames = {theFiles.name};% list of all the files with the given name pattern
folder=theFiles.folder;
T1=readmatrix(fullfile(folder,allFileNames{1}));

% loading data
filePattern = fullfile('20211119-0006','*20211119-0006_*.csv');% name pattern for files
theFiles = dir(filePattern);% properties of the given name pattern
allFileNames = {theFiles.name};% list of all the files with the given name pattern
folder=theFiles.folder;
T=[];
for i=1:length(allFileNames)
    T=[T;readmatrix(fullfile(folder,allFileNames{i}))]; 
end
T1=T1(~isnan(T1(:,2)),:);
t1=T1(:,1)/1000;% [s]
i1=T1(:,2);% [mA]
% sampling frequency
Ts = T(2,1) - T(1,1);
Fs = 1/Ts;
T=T(~isnan(T(:,2)),:);% eliminating NaN values(problem with filtering)
% multiple measures and each measure start from zero,
% when we concatenate the we have to recreate the time scale:
T(:,1)=(0:1:length(T(:,1))-1)*Ts;
% extraction of variables
t=T(:,1);%[s]
i=-T(:,2).*1000;%[mA]
v=-T(:,3);% [v]
% filtering
i_f=lowpass(i,500,Fs,'Steepness',0.95);
i_f=smooth(i,100);
% aligning signals
test=i_f>5;
test1=i1>5;
or=find(test,1);
or1=find(test1,1);
t=t-t(or);
t1=t1-t1(or1);
% cutting to match the length of the other signal
i=i(1:find(t>t1(end),1));
i_f=i_f(1:find(t>t1(end),1));
v=v(1:find(t>t1(end),1));
t=t(1:find(t>t1(end),1));

% psd
[psd,f]=pwelch(i,[],[],[],Fs);
[psd_f,f_f]=pwelch(i_f,[],[],[],Fs);
% plot psd
FigID=FigID+1;
figure(FigID)
    hold on
    plot(f,10*log10(psd),'LineWidth',1)
    plot(f_f,10*log10(psd_f),'LineWidth',1)
    xlabel('Frequency (Hz)')
    ylabel('current PSD (dB/Hz)')
    legend('Unfiltered','Filtered')

% plot signal current
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,i,'LineWidth',1)
    xlabel('time')
    ylabel('current  [mA]')
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,i_f,'LineWidth',1)    
    plot(t1,i1,'LineWidth',1)  
    xlabel('time [s]')
    ylabel('current(filtered) [mA]')
% plot signal tension
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,v,'LineWidth',1)
    xlabel('time [s]')
    ylabel('tension [V]')
% energy computation
P= v.*i_f;% [mW]
E=trapz(t,P);% [mJ]
% plot power
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,P,'LineWidth',1)
    xlabel('time [s]')
    ylabel('power [mW]')
%% 0007 
% 33 mF 5.5V carico, corrente + tensione parte Cap, senza led. Ping ricevuto.
clc 
clearvars
close all
FigID=0;
% importing clean data measired with power supply 
filePattern = fullfile('..\..\ping\boot_ping_total','20210601-0001.csv');% name pattern for files
theFiles = dir(filePattern);% properties of the given name pattern
allFileNames = {theFiles.name};% list of all the files with the given name pattern
folder=theFiles.folder;
T1=readmatrix(fullfile(folder,allFileNames{1}));

% loading data
filePattern = fullfile('20211119-0007.csv');% name pattern for files
theFiles = dir(filePattern);% properties of the given name pattern
allFileNames = {theFiles.name};% list of all the files with the given name pattern
folder=theFiles.folder;
T=[];
for i=1:length(allFileNames)
    T=[T;readmatrix(fullfile(folder,allFileNames{i}))]; 
end
T1=T1(~isnan(T1(:,2)),:);
t1=T1(:,1)/1000;% [s]
i1=T1(:,2);% [mA]
% sampling frequency
Ts = T(2,1) - T(1,1);
Fs = 1/Ts;
T=T(~isnan(T(:,2)),:);% eliminating NaN values(problem with filtering)
% multiple measures and each measure start from zero,
% when we concatenate the we have to recreate the time scale:
T(:,1)=(0:1:length(T(:,1))-1)*Ts;
% extraction of variables
t=T(:,1);%[s]
i=-T(:,2).*1000;%[mA]
v=-T(:,3);% [v]
% filtering
i_f=lowpass(i,500,Fs,'Steepness',0.95);
i_f=smooth(i,100);
% aligning signals
test=i_f>5;
test1=i1>5;
or=find(test,1);
or1=find(test1,1);
t=t-t(or);
t1=t1-t1(or1);
% cutting to match the length of the other signal
i=i(1:find(t>t1(end),1));
i_f=i_f(1:find(t>t1(end),1));
v=v(1:find(t>t1(end),1));
t=t(1:find(t>t1(end),1));

% psd
[psd,f]=pwelch(i,[],[],[],Fs);
[psd_f,f_f]=pwelch(i_f,[],[],[],Fs);
% plot psd
FigID=FigID+1;
figure(FigID)
    hold on
    plot(f,10*log10(psd),'LineWidth',1)
    plot(f_f,10*log10(psd_f),'LineWidth',1)
    xlabel('Frequency (Hz)')
    ylabel('current PSD (dB/Hz)')
    legend('Unfiltered','Filtered')

% plot signal current
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,i,'LineWidth',1)
    xlabel('time')
    ylabel('current  [mA]')
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,i_f,'LineWidth',1)    
    plot(t1,i1,'LineWidth',1)  
    xlabel('time [s]')
    ylabel('current(filtered) [mA]')
% plot signal tension
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,v,'LineWidth',1)
    xlabel('time [s]')
    ylabel('tension [V]')
% energy computation
P= v.*i_f;% [mW]
E=trapz(t,P);% [mJ]
% plot power
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,P,'LineWidth',1)
    xlabel('time [s]')
    ylabel('power [mW]')
%% 0008
% 30 mF 5.5V, corrente + tensione parte Cap, con led. 
clc 
clearvars
close all
FigID=0;
% importing clean data measired with power supply 
filePattern = fullfile('..\..\ping\boot_ping_total','20210601-0001.csv');% name pattern for files
theFiles = dir(filePattern);% properties of the given name pattern
allFileNames = {theFiles.name};% list of all the files with the given name pattern
folder=theFiles.folder;
T1=readmatrix(fullfile(folder,allFileNames{1}));

% loading data
filePattern = fullfile('20211119-0008.csv');% name pattern for files
theFiles = dir(filePattern);% properties of the given name pattern
allFileNames = {theFiles.name};% list of all the files with the given name pattern
folder=theFiles.folder;
T=[];
for i=1:length(allFileNames)
    T=[T;readmatrix(fullfile(folder,allFileNames{i}))]; 
end
T1=T1(~isnan(T1(:,2)),:);
t1=T1(:,1)/1000;% [s]
i1=T1(:,2);% [mA]
% sampling frequency
Ts = T(2,1) - T(1,1);
Fs = 1/Ts;
T=T(~isnan(T(:,2)),:);% eliminating NaN values(problem with filtering)
% multiple measures and each measure start from zero,
% when we concatenate the we have to recreate the time scale:
T(:,1)=(0:1:length(T(:,1))-1)*Ts;
% extraction of variables
t=T(:,1);%[s]
i=-T(:,2).*1000;%[mA]
v=-T(:,3);% [v]
% filtering
i_f=lowpass(i,500,Fs,'Steepness',0.95);
i_f=smooth(i,100);
% psd
[psd,f]=pwelch(i,[],[],[],Fs);
[psd_f,f_f]=pwelch(i_f,[],[],[],Fs);
% plot psd
FigID=FigID+1;
figure(FigID)
    hold on
    plot(f,10*log10(psd),'LineWidth',1)
    plot(f_f,10*log10(psd_f),'LineWidth',1)
    xlabel('Frequency (Hz)')
    ylabel('current PSD (dB/Hz)')
    legend('Unfiltered','Filtered')
% aligning signals
test=i_f>5;
test1=i1>5;
or=find(test,1);
or1=find(test1,1);
t=t-t(or);
t1=t1-t1(or1);
% plot signal current
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,i,'LineWidth',1)
    xlabel('time [s]')
    ylabel('current  [mA]')
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,i_f,'LineWidth',1)    
    plot(t1,i1,'LineWidth',1)  
    xlabel('time [s]')
    ylabel('current (filtered) [mA]')
% plot signal tension
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,v,'LineWidth',1)
    xlabel('time [s]')
    ylabel('tension [V]')

% cutting to match the length of the other signal
i=i(1:find(t>t1(end),1));
i_f=i_f(1:find(t>t1(end),1));
v=v(1:find(t>t1(end),1));
t=t(1:find(t>t1(end),1));
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,i_f,'LineWidth',1)    
    plot(t1,i1,'LineWidth',1)  
    xlabel('time [s]')
    ylabel('current (filtered) [mA]')
% energy computation
P= v.*i_f;% [mW]
E=trapz(t,P);% [mJ]
% plot power
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,P,'LineWidth',1)
    xlabel('time [s]')
    ylabel('power [mW]')

%% 0009
% 50 mF, 5.5V, corr + tens parte superC, con led.
clc 
clearvars
close all
FigID=0;
% importing clean data measired with power supply 
filePattern = fullfile('..\..\ping\boot_ping_total','20210601-0001.csv');% name pattern for files
theFiles = dir(filePattern);% properties of the given name pattern
allFileNames = {theFiles.name};% list of all the files with the given name pattern
folder=theFiles.folder;
T1=readmatrix(fullfile(folder,allFileNames{1}));

% loading data
filePattern = fullfile('20211119-0009.csv');% name pattern for files
theFiles = dir(filePattern);% properties of the given name pattern
allFileNames = {theFiles.name};% list of all the files with the given name pattern
folder=theFiles.folder;
T=[];
for i=1:length(allFileNames)
    T=[T;readmatrix(fullfile(folder,allFileNames{i}))]; 
end
T1=T1(~isnan(T1(:,2)),:);
t1=T1(:,1)/1000;% [s]
i1=T1(:,2);% [mA]
% sampling frequency
Ts = T(2,1) - T(1,1);
Fs = 1/Ts;
T=T(~isnan(T(:,2)),:);% eliminating NaN values(problem with filtering)
% multiple measures and each measure start from zero,
% when we concatenate the we have to recreate the time scale:
T(:,1)=(0:1:length(T(:,1))-1)*Ts;
% extraction of variables
t=T(:,1);%[s]
i=-T(:,2).*1000;%[mA]
v=-T(:,3);% [v]
% filtering
i_f=lowpass(i,500,Fs,'Steepness',0.95);
i_f=smooth(i,100);
% aligning signals
test=i_f>5;
test1=i1>5;
or=find(test,1);
or1=find(test1,1);
t=t-t(or);
t1=t1-t1(or1);
% cutting to match the length of the other signal
i=i(1:find(t>t1(end),1));
i_f=i_f(1:find(t>t1(end),1));
v=v(1:find(t>t1(end),1));
t=t(1:find(t>t1(end),1));

% psd
[psd,f]=pwelch(i,[],[],[],Fs);
[psd_f,f_f]=pwelch(i_f,[],[],[],Fs);
% plot psd
FigID=FigID+1;
figure(FigID)
    hold on
    plot(f,10*log10(psd),'LineWidth',1)
    plot(f_f,10*log10(psd_f),'LineWidth',1)
    xlabel('Frequency (Hz)')
    ylabel('current PSD (dB/Hz)')
    legend('Unfiltered','Filtered')

% plot signal current
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,i,'LineWidth',1)
    xlabel('time')
    ylabel('current  [mA]')
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,i_f,'LineWidth',1)    
    plot(t1,i1,'LineWidth',1)  
    xlabel('time [s]')
    ylabel('current(filtered) [mA]')
% plot signal tension
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,v,'LineWidth',1)
    xlabel('time [s]')
    ylabel('tension [V]')
% energy computation
P= v.*i_f;% [mW]
E=trapz(t,P);% [mJ]
% plot power
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,P,'LineWidth',1)
    xlabel('time [s]')
    ylabel('power [mW]')
%% 0010
% 50 mF, 5.5V, corr + tens parte cap, senza led.( aimentazione staccata e riattaccata per
% errore , doppio ping)
clc 
clearvars
close all
FigID=0;
% importing clean data measired with power supply 
filePattern = fullfile('..\..\ping\boot_ping_total','20210601-0001.csv');% name pattern for files
theFiles = dir(filePattern);% properties of the given name pattern
allFileNames = {theFiles.name};% list of all the files with the given name pattern
folder=theFiles.folder;
T1=readmatrix(fullfile(folder,allFileNames{1}));

% loading data
filePattern = fullfile('20211119-0010.csv');% name pattern for files
theFiles = dir(filePattern);% properties of the given name pattern
allFileNames = {theFiles.name};% list of all the files with the given name pattern
folder=theFiles.folder;
T=[];
for i=1:length(allFileNames)
    T=[T;readmatrix(fullfile(folder,allFileNames{i}))]; 
end
T1=T1(~isnan(T1(:,2)),:);
t1=T1(:,1)/1000;% [s]
i1=T1(:,2);% [mA]
% sampling frequency
Ts = T(2,1) - T(1,1);
Fs = 1/Ts;
T=T(~isnan(T(:,2)),:);% eliminating NaN values(problem with filtering)
% multiple measures and each measure start from zero,
% when we concatenate the we have to recreate the time scale:
T(:,1)=(0:1:length(T(:,1))-1)*Ts;
% extraction of variables
t=T(:,1);%[s]
i=-T(:,2).*1000;%[mA]
v=-T(:,3);% [v]
% filtering
i_f=lowpass(i,500,Fs,'Steepness',0.95);
i_f=smooth(i,100);
% aligning signals
test=i_f>5;
test1=i1>5;
or=find(test,1);
or1=find(test1,1);
t=t-t(or);
t1=t1-t1(or1);
% cutting to match the length of the other signal
i=i(1:find(t>t1(end),1));
i_f=i_f(1:find(t>t1(end),1));
v=v(1:find(t>t1(end),1));
t=t(1:find(t>t1(end),1));

% psd
[psd,f]=pwelch(i,[],[],[],Fs);
[psd_f,f_f]=pwelch(i_f,[],[],[],Fs);
% plot psd
FigID=FigID+1;
figure(FigID)
    hold on
    plot(f,10*log10(psd),'LineWidth',1)
    plot(f_f,10*log10(psd_f),'LineWidth',1)
    xlabel('Frequency (Hz)')
    ylabel('current PSD (dB/Hz)')
    legend('Unfiltered','Filtered')

% plot signal current
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,i,'LineWidth',1)
    xlabel('time')
    ylabel('current  [mA]')
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,i_f,'LineWidth',1)    
    plot(t1,i1,'LineWidth',1)  
    xlabel('time [s]')
    ylabel('current(filtered) [mA]')
% plot signal tension
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,v,'LineWidth',1)
    xlabel('time [s]')
    ylabel('tension [V]')
% energy computation
P= v.*i_f;% [mW]
E=trapz(t,P);% [mJ]
% plot power
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,P,'LineWidth',1)
    xlabel('time [s]')
    ylabel('power [mW]')
%% 0011
% 30 mF, 5.5.V, corr + tens parte cap, senza led, carica con ALPS. da 0.3 V a - V. Non registrato tutto. Controllare.
clc 
clearvars
close all
FigID=0;
% loading data
filePattern = fullfile('20211119-0011','*20211119-0011_*.csv');% name pattern for files
theFiles = dir(filePattern);% properties of the given name pattern
allFileNames = {theFiles.name};% list of all the files with the given name pattern
folder=theFiles.folder;
T=[];
for i=1:length(allFileNames)
    T=[T;readmatrix(fullfile(folder,allFileNames{i}))]; 
end
% sampling frequency
Ts = T(2,1) - T(1,1);
Fs = 1/Ts;
% time scale adjusing
T(:,1)=(0:1:length(T(:,1))-1)*Ts;
T=T(~isnan(T(:,2)),:);% eliminating NaN values(problem with filtering)
% extraction of variables
t=T(:,1);%[s]
i=T(:,2).*1000;%[mA]
v=-T(:,3);% [v]
puls=T(:,4); % [v]
% plot signal current
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,i,'LineWidth',1)
    xlabel('time')
    ylabel('current  [mA]')
% plot signal tension
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,v,'LineWidth',1)
    xlabel('time [s]')
    ylabel('tension [V]') 
% pressures count
puls=puls(~isnan(puls));%removing nan values due to channel ovverange
t=t(~isnan(puls));
pres = (puls>3);% trasform il segnale in binario e dove è maggiore di 3 lo metto a 1 e dove è minore a 0 in modo da eliminare micropicchi
pres = double(pres);
% ci sono ancora le oscillazioni attorno al 3 da eliminare, per farlo
% prendo picchi che abbiano un minimo di distanza fra loro(in qunato le
% oscillazioni le si trova in genere fra una campionatura e l'altra e quindi sono mpolo vicine):
[pk,lk] = findpeaks(pres,'MinPeakDistance',1000);
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,puls,t(lk),puls(lk),'o','LineWidth',1)
pressures=length(pk);% valore approssimativo ma moto vicino
%% 0012 
% 30 mF, 5.5V, corr, tens parte cap + tasto, senza led, carica da 3.9 V a V accensione con ALPS. Ping ricevuto.
clc 
clearvars
close all
FigID=0;
% loading data
filePattern = fullfile('20211119-0012','*20211119-0012_*.csv');% name pattern for files
theFiles = dir(filePattern);% properties of the given name pattern
allFileNames = {theFiles.name};%l ist of all the files with the given name pattern
folder=theFiles.folder;
T=[];
for i=1:length(allFileNames)
    T=[T;readmatrix(fullfile(folder,allFileNames{i}))]; 
end
% sampling frequency
Ts = T(2,1) - T(1,1);
Fs = 1/Ts;
% time scale adjusing
T(:,1)=(0:1:length(T(:,1))-1)*Ts;
T=T(~isnan(T(:,2)),:);% eliminating NaN values(problem with filtering)
% extraction of variables
t=T(:,1);%[s]
i=T(:,2).*1000;%[mA]
v=-T(:,3);% [v]
puls=T(:,4); % [v]
% plot signal current
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,i,'LineWidth',1)
    xlabel('time')
    ylabel('current  [mA]')
% plot signal tension
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,v,'LineWidth',1)
    xlabel('time [s]')
    ylabel('tension [V]') 
% pressures count
puls=puls(~isnan(puls));%removing nan values due to channel ovverange
t=t(~isnan(puls));
pres = (puls>3);% trasform il segnale in binario e dove è maggiore di 3 lo metto a 1 e dove è minore a 0 in modo da eliminare micropicchi
pres = double(pres);
% ci sono ancora le oscillazioni attorno al 3 da eliminare, per farlo
% prendo picchi che abbiano un minimo di distanza fra loro(in qunato le
% oscillazioni le si trova in genere fra una campionatura e l'altra e quindi sono mpolo vicine):
[pk,lk] = findpeaks(pres,'MinPeakDistance',1000);
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,puls,t(lk),puls(lk),'o','LineWidth',1)
pressures=length(pk);% valore approssimativo ma moto vicino
%% 0013
% bmp

%% 0014
% 30 mF, 5.5V, corr tens parte cap + tasto, senza led, carica da 3.9 V con AFIG. Ping ricevuto.
clc 
clearvars
close all
FigID=0;
% loading data
filePattern = fullfile('20211119-0014','*20211119-0014_*.csv');% name pattern for files
theFiles = dir(filePattern);% properties of the given name pattern
allFileNames = {theFiles.name};% list of all the files with the given name pattern
folder=theFiles.folder;
T=[];
for i=1:length(allFileNames)
    T=[T;readmatrix(fullfile(folder,allFileNames{i}))]; 
end
% sampling frequency
Ts = T(2,1) - T(1,1);
Fs = 1/Ts;
% time scale adjusing
T(:,1)=(0:1:length(T(:,1))-1)*Ts;
T=T(~isnan(T(:,2)),:);% eliminating NaN values(problem with filtering)
% extraction of variables
t=T(:,1);%[s]
i=T(:,2).*1000;%[mA]
v=-T(:,3);% [v]
puls=T(:,4); % [v]
% plot signal current
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,i,'LineWidth',1)
    xlabel('time')
    ylabel('current  [mA]')
% plot signal tension
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,v,'LineWidth',1)
    xlabel('time [s]')
    ylabel('tension [V]') 
% pressures count
puls=puls(~isnan(puls));%removing nan values due to channel ovverange
t=t(~isnan(puls));
pres = (puls>3);% trasform il segnale in binario e dove è maggiore di 3 lo metto a 1 e dove è minore a 0 in modo da eliminare micropicchi
pres = double(pres);
% ci sono ancora le oscillazioni attorno al 3 da eliminare, per farlo
% prendo picchi che abbiano un minimo di distanza fra loro(in qunato le
% oscillazioni le si trova in genere fra una campionatura e l'altra e quindi sono mpolo vicine):
[pk,lk] = findpeaks(pres,'MinPeakDistance',1000);
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,puls,t(lk),puls(lk),'o','LineWidth',1)
pressures=length(pk);% valore approssimativo ma moto vicino
%% 0015
% 30 mF, 5.5V, corr tens parte cap + tasto, senza led, carica da 0 a V accensione con AFIG. Ping ricevuto.
clc 
clearvars
close all
FigID=0;
% loading data
filePattern = fullfile('20211119-0015','*20211119-0015_*.csv');% name pattern for files
theFiles = dir(filePattern);% properties of the given name pattern
allFileNames = {theFiles.name};% list of all the files with the given name pattern
folder=theFiles.folder;
T=[];
for i=1:length(allFileNames)
    T=[T;readmatrix(fullfile(folder,allFileNames{i}))]; 
end
% sampling frequency
Ts = T(2,1) - T(1,1);
Fs = 1/Ts;
% time scale adjusing
T(:,1)=(0:1:length(T(:,1))-1)*Ts;
T=T(~isnan(T(:,2)),:);% eliminating NaN values(problem with filtering)
% extraction of variables
t=T(:,1);%[s]
i=T(:,2).*1000;%[mA]
v=-T(:,3);% [v]
puls=T(:,4); % [v]
% plot signal current
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,i,'LineWidth',1)
    xlabel('time')
    ylabel('current  [mA]')
% plot signal tension
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,v,'LineWidth',1)
    xlabel('time [s]')
    ylabel('tension [V]') 
% pressures count
puls=puls(~isnan(puls));%removing nan values due to channel ovverange
t=t(~isnan(puls));
pres = (puls>2);% trasform il segnale in binario e dove è maggiore di 3 lo metto a 1 e dove è minore a 0 in modo da eliminare micropicchi
pres = double(pres);
% ci sono ancora le oscillazioni attorno al 3 da eliminare, per farlo
% prendo picchi che abbiano un minimo di distanza fra loro(in qunato le
% oscillazioni le si trova in genere fra una campionatura e l'altra e quindi sono mpolo vicine):
[pk,lk] = findpeaks(pres,'MinPeakDistance',1000);
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,puls,t(lk),puls(lk),'o','LineWidth',1)
    xlabel('time [s]')
    ylabel('voltage [v]')
pressures=length(pk);% valore approssimativo ma moto vicino
%% 0016
% bmp

%% 0017
% ping pong 100mF con led di segnalazione e led alimentazione, misura corrente + tensione supercap
clc 
clearvars
close all
FigID=0;
% loading data
filePattern = fullfile('20211119-0017','*20211119-0017_*.csv');% name pattern for files
theFiles = dir(filePattern);% properties of the given name pattern
allFileNames = {theFiles.name};% list of all the files with the given name pattern
folder=theFiles.folder;
T=[];
for i=1:length(allFileNames)
    T=[T;readmatrix(fullfile(folder,allFileNames{i}))]; 
end
% sampling frequency
Ts = T(2,1) - T(1,1);
Fs = 1/Ts;
% time scale adjusing
T(:,1)=(0:1:length(T(:,1))-1)*Ts;
T=T(~isnan(T(:,2)),:);% eliminating NaN values(problem with filtering)
% extraction of variables
t=T(:,1);%[s]
i=-T(:,2).*1000;%[mA]
v=T(:,3);% [v]
% filtering
%i_f=lowpass(i,500,Fs,'Steepness',0.95);
i_f=smooth(i,100);
% psd
[psd,f]=pwelch(i,[],[],[],Fs);
[psd_f,f_f]=pwelch(i_f,[],[],[],Fs);
% plot psd
FigID=FigID+1;
figure(FigID)
    hold on
    plot(f,10*log10(psd),'LineWidth',1)
    plot(f_f,10*log10(psd_f),'LineWidth',1)
    xlabel('Frequency (Hz)')
    ylabel('current PSD (dB/Hz)')
    legend('Unfiltered','Filtered')
% plot signal current
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,i,'LineWidth',1)
    xlabel('time')
    ylabel('current  [mA]')
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,i_f,'LineWidth',1)    
    xlabel('time [s]')
    ylabel('current(filtered) [mA]')
% plot signal tension
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,v,'LineWidth',1)
    xlabel('time [s]')
    ylabel('tension [V]')
%% 0018
% ping pong 100mF diverso sampling time
clc 
clearvars
close all
FigID=0;
% loading data
filePattern = fullfile('20211119-0018.csv');% name pattern for files
theFiles = dir(filePattern);% properties of the given name pattern
allFileNames = {theFiles.name};% list of all the files with the given name pattern
folder=theFiles.folder;
T=[];
for i=1:length(allFileNames)
    T=[T;readmatrix(fullfile(folder,allFileNames{i}))]; 
end
% sampling frequency
Ts = T(2,1) - T(1,1);
Fs = 1/Ts;
% time scale adjusing
T(:,1)=(0:1:length(T(:,1))-1)*Ts;
T=T(~isnan(T(:,2)),:);% eliminating NaN values(problem with filtering)
% extraction of variables
t=T(:,1);%[s]
i=-T(:,2).*1000;%[mA]
v=T(:,3);% [v]
% filtering
%i_f=lowpass(i,500,Fs,'Steepness',0.95);
i_f=smooth(i,100);
% psd
[psd,f]=pwelch(i,[],[],[],Fs);
[psd_f,f_f]=pwelch(i_f,[],[],[],Fs);
% plot psd
FigID=FigID+1;
figure(FigID)
    hold on
    plot(f,10*log10(psd),'LineWidth',1)
    plot(f_f,10*log10(psd_f),'LineWidth',1)
    xlabel('Frequency (Hz)')
    ylabel('current PSD (dB/Hz)')
    legend('Unfiltered','Filtered')
% plot signal current
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,i,'LineWidth',1)
    xlabel('time')
    ylabel('current  [mA]')
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,i_f,'LineWidth',1)    
    xlabel('time [s]')
    ylabel('current(filtered) [mA]')
% plot signal tension
FigID=FigID+1;
figure(FigID)
    hold on    
    plot(t,v,'LineWidth',1)
    xlabel('time [s]')
    ylabel('tension [V]')