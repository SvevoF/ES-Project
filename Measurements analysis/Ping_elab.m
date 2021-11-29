clc 
clearvars
close all
filter_en=1;
%% Loading data
filePattern1 = fullfile({'boot_ping_lra';'boot_ping_mcu';'boot_ping_total'},'20210601-0001.csv');% name pattern for files
%% LRA

    % data extraction
    theFiles = dir(filePattern1{1});% properties of the given name pattern
    allFileNames = {theFiles.name};%list of all the files with the given name pattern
    folder=theFiles.folder;
    T=readmatrix(fullfile(folder,allFileNames{1}));%reading file
    t_lra=T(:,1)/1000; % [S]
    I_lra=T(:,2); % [mA]
    V_lra=T(:,3); % [V]
    % energy
    E_lra=-trapz(t_lra,I_lra.*V_lra);%  [mJ]



    % frequency extraction
    Ts = t_lra(2) - t_lra(1);
    Fs = 1/Ts;

    % filter
    d = designfilt('bandstopiir','FilterOrder',2, ...
        'HalfPowerFrequency1',49,'HalfPowerFrequency2',51, ...
        'DesignMethod','butter','SampleRate',Fs);
    %plot
    FigID=0;
    FigID=FigID+1;
    fvtool(d,'Fs',Fs)
    title('filter lra')

    % filtering
    I_lra_filt=filtfilt(d,I_lra);
    % plot
    FigID=FigID+1;
    figure(FigID)
    plot(t_lra,I_lra,t_lra,I_lra_filt,'LineWidth',1)
    ylabel('Current (A)')
    xlabel('Time (s)')
    title('Signal filtering lra')
    legend('Unfiltered','Filtered')
    grid

    % psd
    [psd,f]=pwelch(I_lra,[],[],[],Fs);
    [psd_f,f_f]=pwelch(I_lra_filt,[],[],[],Fs);
    % plot
    FigID=FigID+1;
    figure(FigID)
    hold on
    plot(f(f<1000),10*log10(psd(f<1000)),'LineWidth',1)
    plot(f_f(f_f<1000),10*log10(psd_f(f_f<1000)),'LineWidth',1)
    xlabel('Frequency (Hz)')
    ylabel('PSD (dB/Hz)')
    legend('Unfiltered','Filtered')
    title('Power spectral density lra')

    %energy
    E_lra_filt=-trapz(t_lra,I_lra_filt.*V_lra);% [mJ]
   
%% MCU
    % data extraction
    theFiles = dir(filePattern1{2});% properties of the given name pattern
    allFileNames = {theFiles.name};%list of all the files with the given name pattern
    folder=theFiles.folder;
    T=readmatrix(fullfile(folder,allFileNames{1}));%reading file
    t_mcu=T(:,1); % [S]
    I_mcu=T(:,2); % [mA]
    V_mcu=T(:,3); % [V]
    % energy
    E_mcu=-trapz(t_mcu,I_mcu.*V_mcu);%  [mJ]



    % frequency extraction
    Ts = t_mcu(2) - t_mcu(1);
    Fs = 1/Ts;

    % filter
    d = designfilt('bandstopiir','FilterOrder',2, ...
        'HalfPowerFrequency1',49,'HalfPowerFrequency2',51, ...
        'DesignMethod','butter','SampleRate',Fs);
    %plot
    FigID=FigID+1;
    fvtool(d,'Fs',Fs)
    title(' filter mcu')

    % filtering
    I_mcu_filt=filtfilt(d,I_mcu);
    % plot
    FigID=FigID+1;
    figure(FigID)
    plot(t_mcu,I_mcu,t_mcu,I_mcu_filt,'LineWidth',1)
    ylabel('Current (A)')
    xlabel('Time (s)')
    title('Signal filtering mcu')
    legend('Unfiltered','Filtered')
    grid

    % psd
    [psd,f]=pwelch(I_mcu,[],[],[],Fs);
    [psd_f,f_f]=pwelch(I_mcu_filt,[],[],[],Fs);
    % plot
    FigID=FigID+1;
    figure(FigID)
    hold on
    plot(f(f<1000),10*log10(psd(f<1000)),'LineWidth',1)
    plot(f_f(f_f<1000),10*log10(psd_f(f_f<1000)),'LineWidth',1)
    xlabel('Frequency (Hz)')
    ylabel('PSD (dB/Hz)')
    legend('Unfiltered','Filtered')
    title('Power spectral density mcu')

    %energy
    E_mcu_filt=-trapz(t_mcu,I_mcu_filt.*V_mcu);% [mJ]
    
%% TOTAL
    % data extraction
    theFiles = dir(filePattern1{3});% properties of the given name pattern
    allFileNames = {theFiles.name};%list of all the files with the given name pattern
    folder=theFiles.folder;
    T=readmatrix(fullfile(folder,allFileNames{1}));%reading file
    t_tot=T(:,1)/1000; % [S]
    I_tot=T(:,2); % [mA]
    % energy
    E_tot=trapz(t_tot,I_tot.*3.3);%  [mJ]



    % frequency extraction
    Ts = t_tot(2) - t_tot(1);
    Fs = 1/Ts;

    % filter
    d = designfilt('bandstopiir','FilterOrder',2, ...
        'HalfPowerFrequency1',47,'HalfPowerFrequency2',51, ...
        'DesignMethod','butter','SampleRate',Fs);
    %plot
    FigID=FigID+1;
    fvtool(d,'Fs',Fs)
    title('filter total')
    % filtering
    I_tot_filt=filtfilt(d,I_tot);
    % plot
    FigID=FigID+1;
    figure(FigID)
    plot(t_tot,I_tot,t_tot,I_tot_filt,'LineWidth',1)
    ylabel('Current (A)')
    xlabel('Time (s)')
    title('Signal filtering total')
    legend('Unfiltered','Filtered')
    grid

    % psd
    [psd,f]=pwelch(I_tot,[],[],[],Fs);
    [psd_f,f_f]=pwelch(I_tot_filt,[],[],[],Fs);
    % plot
    FigID=FigID+1;
    figure(FigID)
    hold on
    plot(f(f<1000),10*log10(psd(f<1000)),'LineWidth',1)
    plot(f_f(f_f<1000),10*log10(psd_f(f_f<1000)),'LineWidth',1)
    xlabel('Frequency (Hz)')
    ylabel('PSD (dB/Hz)')
    legend('Unfiltered','Filtered')
    title('Power spectral density total')

    %energy
    E_tot_filt=trapz(t_tot,I_tot_filt.*3.3);% [mJ]
%% signal allining
test=I_lra_filt>3;
test1=I_mcu_filt>3;
test2=I_tot_filt>3;
or=find(test,1);
or1=find(test1,1);
or2=find(test2,1);
t_lra=t_lra-t_lra(or);
t_mcu=t_mcu-t_mcu(or1);
t_tot=t_tot-t_tot(or2);
I_mcu_filt=I_mcu_filt(t_mcu>t_lra(1) & t_mcu<t_lra(end));
V_mcu=V_mcu(t_mcu>t_lra(1) & t_mcu<t_lra(end));
t_mcu=t_mcu(t_mcu>t_lra(1) & t_mcu<t_lra(end));
I_tot_filt=I_tot_filt(t_tot>t_lra(1) & t_tot<t_lra(end));
t_tot=t_tot(t_tot>t_lra(1)& t_tot<t_lra(end));

%% recoputing energy with alligned signals
E_tot_al=trapz(t_tot,I_tot_filt.*3.3);% [mJ]
E_mcu_al=-trapz(t_mcu,I_mcu_filt.*V_mcu);% [mJ]
E_lra_al=-trapz(t_lra,I_lra_filt.*V_lra);% [mJ]
%% interpolationg signal 
Ts=t_lra(2)-t_lra(1);
I_mcu_int=interp1(t_mcu,I_mcu_filt,t_lra);
%% Final Plots
FigID=FigID+1;
figure(FigID)
    hold on
    plot(t_tot,I_tot_filt,'LineWidth',1)
    plot(t_lra,I_lra_filt,'LineWidth',1)
    plot(t_mcu,I_mcu_filt,'LineWidth',1)
    plot(t_lra,I_mcu_int+I_lra_filt,'LineWidth',1)
    ylabel('Current (mA)')
    xlabel('Time (s)')
    title('Mesaured signals')
    legend('total','lra','mcu','lra+mcu')
    
FigID=FigID+1;
figure(FigID)
    X = categorical({'lra','mcu','sum','total mes'});
    Y = [E_lra_filt;E_mcu_filt; E_lra_filt+E_mcu_filt;E_tot_filt];
    bar(X,Y,'stacked')
    title('energy overestimated')
    ylabel('E (mJ)')

FigID=FigID+1;
    figure(FigID)
    X = categorical({'lra','mcu','mcu+lra','total mes'});
    Y = [E_lra_al;E_mcu_al; E_lra_al+E_mcu_al;E_tot_al];
    bar(X,Y,'stacked')
    title('energy consume')
    ylabel('E (mJ)')

