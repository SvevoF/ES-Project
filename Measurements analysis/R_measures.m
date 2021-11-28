%% Power consumption analysis
clc 
clearvars
close all
%% Loading data

myFolder = 'misure_harvester';% starting folder
filePattern1 = fullfile(myFolder,'AFIG',{'10_ohm';'100_ohm';'1000_ohm'},'20210413-0002', '*20210413-0002_*.csv');% name pattern for files AFIG 
filePattern2 = fullfile(myFolder,'ALPS',{'10_ohm';'100_ohm';'1000_ohm'},'20210413-0002', '*20210413-0002_*.csv');% name pattern for files ALPS

%% energy
AFIG=zeros(3,1);% we have 3 values of resistance 10 100 1000 ohm at wich we measure voltage
ALPS=zeros(3,1);
%AFIG 
for j=1:3
    theFiles = dir(filePattern1{j,1});% properties of the given name pattern
    allFileNames = {theFiles.name};%list of all the files with the given name pattern
    folder=theFiles.folder;
    E = zeros(10,1);
    
    for i=1:10 % for each resistance we take 10 measurements of the generated voltage 
        T=readmatrix(fullfile(folder,allFileNames{1,i}));%reading file
        % considering only voltage >0.6 to account rectifier diode
        V=T(T(:,2)>0.7,:);
        E(i)=trapz(V(:,1)/1000,(V(:,2).^2)./(10^j));%computing energy integrating the power
    end
    
    AFIG(j)=mean(E);%mean energy
end
%ALPS 
for j=1:3
    theFiles = dir(filePattern2{j,1});
    allFileNames = {theFiles.name};
    folder=theFiles.folder;
    E = zeros(10,1);
    for i=1:10
        T=readmatrix(fullfile(folder,allFileNames{1,i}));
        if j==1 % for the firstmeasure we use millivolts since we have small voltages so we have to convert
            div=1000;
        else
            div=1;
        end
        V=T(T(:,2)>0.7*div,:);
        E(i)=trapz(V(:,1)/1000,((V(:,2)/div).^2)./(10^j));
    end
    ALPS(j)=mean(E);
end

%% Plots
figure(1), clf
hold on
plot([10 100 1000], AFIG*1e6,'o','linewidth',1.5);
plot([10 100 1000], ALPS*1e6,'o','linewidth',1.5);
xlabel('resistance [ohm]')
ylabel('energy [uJ]')
legend('AFIG','ALPS')
set(gca,'XScale','log')