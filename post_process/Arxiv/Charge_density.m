function [carrier_dens]=Charge_density

waitfor(msgbox(['Select extracted E field']));
[file,path] = uigetfile('*.mat');
a=load([path '\' file]');

%L=who('E*');
L=fieldnames(a);
[indx,tf] = listdlg('PromptString',{'Select a file.',...
     'Only one file can be selected at a time.',''},...
     'SelectionMode','multiple','ListString',L);


eps_CZT=10.9;
eps_0=8.8542e-12; %F/m
eps=eps_CZT*eps_0;
e_charge=1.6e-19;

%E=E_900v_25mA_7336_98_9;
E=a.(char(L(indx)));
diff_E=diff(E);
dEdx=diff_E./(10^-5); %Calculates dE/dx
charge_dens=dEdx.*eps./10^6; %C/cm^3
carrier_dens=charge_dens./e_charge;

figure
plot((1:length(smooth(carrier_dens,8)))./100,smooth(carrier_dens,8))
%xlim([1 3.2])
hold on
plot((1:length(carrier_dens))./100,carrier_dens)
title(['Carrier density across the sensor' char(L(indx))])
ylabel('carrier density (/cm^3)')
xlabel('cm')

figure
plot((1:length(smooth(E,8)))./100,smooth(E./100,8))
%xlim([1 3.2])
hold on
plot((1:length(E))./100,E./100)
title(['Electic field' char(L(indx))])
ylabel('E field (V/cm)')
xlabel('cm')


end