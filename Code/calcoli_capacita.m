%% iniialization
clc;
clearvars;
close all;
fprintf('_____________________________________________________________________\n')
fprintf('All the following computations are done for the following devices:\n')
fprintf('Alps LTC3588-1 energy harvester configured for a 3.3 [V] output \n')
fprintf('AFIG-0007 Energy Generator \n')
fprintf('The computations are done negleting the power losses through the wire\n')
fprintf('The enegy values comes from datasheets and our measurements\n')

%% Capacity computation
fprintf('_____________________________________________________________________\n')

Vr_min=4.73;
Vf_max=3.91;
Vr_max=5.37;
Vf_min=3.42;
Dv=Vr_min-Vf_max;
E=0.03705;% sbagliata (energia totale ma sovrastimata)
E=0.00648;% sbagliata(energia solo della raadio ma sovrastimata)
E=0.002157;% energia consumata da un ping , presa dal altro file di elaborazione
C=E*2/Dv^2;
Cs=C;
Vs=5.5;

fprintf('Computing the capacity value:\n')
fprintf('The minimum rising voltage treshold is %4.2f [V]\n',Vr_min)
fprintf('The maximum falling voltage treshold is %4.2f [V]\n',Vf_max)
fprintf('The maximum rising voltage treshold is %4.2f [V]\n',Vr_max)
fprintf('The minimum falling voltage treshold is %4.2f [V]\n',Vf_min)
fprintf('In the worst case scenario the voltage difference is %4.2f [V]\n',Dv)
fprintf('The energy required to send a simple ping is %8.6f [J]\n',E)
fprintf('The maximum capacity values is %5.4f [F] \n',C)
fprintf('The super capacitor selected from datasheet has C=%7.6f [F] V=%3.2f [V] \n',Cs,Vs)
%% Estimation of the number of button presses
% calcolo motlo approssimativo vito che l'energi assorbita dipende dalla
% tensione a cui è il condensatore(come visto anche l'harvester da diverse
% tensioni ed energie in base al carico)
fprintf('_____________________________________________________________________\n')

E_press=0.000122;
E_full=Cs*Vr_max^2/2;
E_norm=Cs*(Vr_max-Vf_min)^2/2;
press_full=ceil(E_full/E_press);
press_norm=ceil(E_norm/E_press);

fprintf('Compunting the necessary number of button presses to charge the Capacitor:\n')
fprintf('The minimum energy produced by a single press of the button is %8.7f [J] \n',E_press)
fprintf('In the first cycle for the worst case scenario we have to reach a tension of %5.3f [V]\n',Vr_max)
fprintf('In a normal cycle for the worst case scenario we have to increase the tensios  by %5.3f [V]\n',Vr_max-Vf_min)
fprintf('In the first cycle for the worst case scenario we have to press the button at least %i times\n',press_full)
fprintf('In the normal cycle for the worst case scenario we have to press the button at least %i times\n',press_norm)





