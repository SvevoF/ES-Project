%% Power consumption analysis
clc 
clearvars
close all
%% Loading data
filePattern= fullfile({'07_0';'07_1';'07_2';'08_0';'08_1';'08_2';'09_0';'09_1';'09_2';'10_0';'10_1';'10_2';'11_0';'11_1';'11_2';'12_0';'12_1';'12_2'}...
, '20210525-0001.csv');% name pattern

%% energy
mean_power=zeros(6,3);
energy=zeros(6,3);
i=0;
for j=1:18
    if i==3
        i=1;
    else
        i=i+1;
    end
    if j==13||j==16||j==17
        continue
    end
    theFiles = dir(filePattern{j});% properties of the given name pattern
    allFileNames = {theFiles.name};%list of all the files with the given name pattern
    folder=theFiles.folder;   
    T=readmatrix(fullfile(folder,allFileNames{1}));%reading file    
    V=T(T(:,2)>0,:);
    mean_power(ceil(j/3),i)=mean((V(:,2))*3.3);
    energy(ceil(j/3),i)=trapz(V(:,1)/1000,(V(:,2))*3.3);
end

%% plots
%power
figure(1)
plot(7:12, mean_power,'d','LineWidth',1)
legend('bandwidth: 125','bandwidth: 250','bandwidth: 500')
xlabel('Spreading Factor')
ylabel('Power [mW]')
title('Mean power consumption in 1 second measure')
curtick = get(gca, 'xTick');
xticks(unique(round(curtick)));
%energy
figure(2)
hold on
sz=100;
scatter(7:12,energy(:,1),sz,'o','filled')
scatter(7:12,energy(:,2),sz,'o','filled')
scatter(7:12,energy(:,3),sz,'o','filled')
legend('bandwidth: 125 kHz','bandwidth: 250 kHz','bandwidth: 500 kHz','Location','southwest')
xlabel('Spreading Factor')
ylabel('Energy [mJ]')
title('Energy consumption in 1 second measure')
curtick = get(gca, 'xTick');
xticks(unique(round(curtick)));

% we get the same values for energy and mean power because we measure for 1
% second