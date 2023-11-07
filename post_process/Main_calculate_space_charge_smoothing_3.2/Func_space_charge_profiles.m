% calculate net space charge profile based on E-field profile from Pockels
% for (1)dark and (2) under X-ray at the same bias, calculate the space
% charge profiles and then calculate the net space charge profile
% 2021-10-20 Created by Yuxin
%------------------!!!----------------------------------------------------
% the difference of this code from Func_space_charge_profiles_mu_tau_h.m
% is: a 5th order of polynominal fitting is implemented for both dark and
% Xray data
%-------------------------------------------------------------------------

function [rho_net,rho_net_mean,rho_dark, rho_Xray,mu_tau_h,xx]=Func_space_charge_profiles(bias,current)
%load data with selected bias and X-ray tube current into workspace
Bias=num2str(bias);
Current=num2str(current);
string1=['*' Bias 'V_' Current 'mA' '*output.mat'];
S = dir(fullfile(string1));
load(S.name);
Output_Xray=Output;

string2=['*' Bias 'V_0mA' '*output.mat'];
S = dir(fullfile(string2));
load(S.name);
Output_dark=Output;

string3=['*calib*.mat'];
S = dir(fullfile(string3));
load(S.name);

%% calculate space charge profile
[rho_dark,xx_dark]=Func_space_charge_line_profile(Output_dark,Calib);%[e/cm^3]
[rho_Xray,xx_Xray]=Func_space_charge_line_profile(Output_Xray,Calib);%[e/cm^3]

% if part of the E-field profile was selected for curve fitting, recsulting
% in different length of rho_dark and rho_Xray

if length(xx_dark)~=length(xx_Xray)
    [xx,position_dark,position_Xray]=intersect(xx_dark,xx_Xray);
    rho_dark=rho_dark(position_dark(1:end-1));
    rho_Xray=rho_Xray(position_Xray(1:end-1));
else
    xx=xx_dark;
end
    
rho_net=rho_Xray-rho_dark;%[e/cm^3]% net space charge profile caused by X-ray
rho_net_mean=mean(rho_net,'omitnan');%[e/cm^3]% average net charge density caused by X-ray
%% calculate mu.tau-h
% [mu_tau_h,errorcode]=Func_mu_tau_h(rho_net,xx);%[cm^2/V] 
[mu_tau_h]=Func_mu_tau_h_method2(rho_net_mean,current);%[cm^2/V] 
%% plot figure
% figure,hold on
f93=figure(93);
clf(f93)
movegui(f93,[1184 575]);%movegui(f1,[40 903]);
% xx=Calib.x_sensor(1:end-1);
plot(xx(1:end-1),rho_dark,'b','displayname','dark')
hold on
plot(xx(1:end-1),rho_Xray,'g','displayname','X-ray')
plot(xx(1:end-1),rho_net,'r','LineWidth',2,'displayname','net space charge')
ylim([-1.5E11 1.5E11])
zeromark=zeros(1,length(xx)-1);
plot(xx(1:end-1),zeromark,'--k','displayname','zero charge')

xlabel('position [mm]')
ylabel('Space charge [e/cm^3]')
legend
title(['Space charge profiles-' Output_dark.sensor_name '-' Bias 'V,' Current 'mA'])
str1=num2str(mu_tau_h);
% str2=num2str(mu_tau_h2);
% annotation('textbox',[0.15260434458171 0.13015873015873 0.346205179227814 0.0801739926739944],'String',{['mu.tau-h=' str1 '[cm^2/V]']},'FitBoxToText','off');
annotation('textbox',[0.141890058867424 0.13968253968254 0.346205179227814 0.116681929181931],'String',{['mu.tau-h=' str1 '[cm^2/V]']},'FitBoxToText','off');
grid
box
hold off
savefig([Output_dark.sensor_name '_' Bias 'V_' Current 'mA_Space_Charge_profiles'])

%% save workspace
filename=[Output_dark.sensor_name '_' Bias 'V_' Current 'mA' '_Space_charge.mat'];
save(filename,'rho_dark','rho_net','rho_net_mean','rho_Xray','mu_tau_h')



end