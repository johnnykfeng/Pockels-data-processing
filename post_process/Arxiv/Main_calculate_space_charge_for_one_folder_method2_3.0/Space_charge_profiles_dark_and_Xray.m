%calculate net space charge profile based on E-field profile from Pockels
%2021-10-20 Created by Yuxin

% required variables in the Workspace:
% Output_dark
% Output_Xray
% Calib
[rho_dark]=Func_space_charge_line_profile(Output_dark,Calib);%[e/cm^3]
[rho_Xray]=Func_space_charge_line_profile(Output_Xray,Calib);%[e/cm^3]
rho_net=rho_Xray-rho_dark;%[e/cm^3]% net space charge profile caused by X-ray
rho_net_mean=mean(rho_net);%[e/cm^3]% average net charge density caused by X-ray
%% plot figure
figure,hold on
xx=Calib.x_sensor(1:end-1);
plot(xx,rho_dark,'b','displayname','dark')
plot(xx,rho_Xray,'g','displayname','X-ray')
plot(xx,rho_net,'r','LineWidth',2,'displayname','net space charge')

zeromark=zeros(1,length(xx));
plot(xx,zeromark,'--k','displayname','zero charge')

xlabel('position [mm]')
ylabel('Space charge [e/cm^3]')
legend
% title('')
grid
box
