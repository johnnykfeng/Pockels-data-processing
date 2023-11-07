%-------------------------------------------------------------------------
% For calculating the space charge profile and mu.tau-h for all
% measurements within the same folder
%------------------!!!----------------------------------------------------
% the difference of method 2 compared to method 1
% is: a 5th order of polynominal fitting is implemented for both dark and
% Xray data
%-------------------------------------------------------------------------
% Step1: Check whether the P.txt file exists with all biases and currents
% Step2: set the working directory to the folder with Pockels raw data
% Step3: run this code, choose "add to path" if directory selection window
% pops up

% Version 3.0
% modification hisory see the end
%-------------------------------------------------------------------------
% addpath 'C:\Users\10404\OneDrive - Redlen\Documents\MATLAB'
% addpath('C:\Users\10404\OneDrive - Redlen\Documents\MATLAB\Pockels\post_process\Space_charge')
% addpath('C:\Users\10404\OneDrive - Redlen\Documents\MATLAB\Pockels\post_process')
% addpath('C:\Users\10404\OneDrive - Redlen\Documents\MATLAB\Pockels')
%%
clear
load 'P.txt'
% P=load('P_test.txt')
Version= '3.0';
Psize=size(P);
N=Psize(1);
rho_dark_all=NaN(N,1);
rho_Xray_all=NaN(N,1);
rho_net_all=NaN(N,1);
mu_tau_h_all=NaN(N,1);
% mu_tau_h_all2=NaN(N,1);
Bias=num2str(P(1,1));
Current=num2str(P(1,2));
string1=['*' Bias 'V_' Current 'mA' '*output.mat'];
S = dir(fullfile(string1));
if size(S)>1
    fprintf(2,'Warning:there are more than 1 data set with the same bias and current');
end

load(S.name);
sensor_name=Output.sensor_name;

for i=1:N
    if P(i,2)>0
        [rho_net,rho_net_mean,rho_dark, rho_Xray,mu_tau_h,xx]=Func_space_charge_profiles(P(i,1),P(i,2));
        rho_dark_all(i)=mean(rho_dark(20:end-20));
        if i-1>0
            if P(i-1,2)==0
                rho_dark_all(i-1)=mean(rho_dark(20:end-20));
            end
        end
        rho_Xray_all(i)=mean(rho_Xray(20:end-20));
        rho_net_all(i)=mean(rho_net(20:end-20));
        mu_tau_h_all(i)=mu_tau_h;
%         mu_tau_h_all2(i)=mu_tau_h2;
    else
        bias_str=num2str(P(i,1));
        string2=['*' bias_str 'V_0mA' '*.mat'];
        S = dir(fullfile(string2));
        load(S.name);
        Output_dark=Output;

        string3=['*calib*.mat'];
        S = dir(fullfile(string3));
        load(S.name);

        [rho_dark,xx_dark]=Func_space_charge_line_profile(Output_dark,Calib);%[e/cm^3]       
        rho_dark_all(i)=mean(rho_dark(20:end-20));
    end
end

%% plot figures
% mu.tau-h vs bias
N_currents=unique(P(:,2));
N_bias=unique(P(:,1));
fig = figure;
hold on
for i=1:length(N_bias)
    for j=1:length(N_currents)
        bias_axis_add(j)=(j-2)*10;
%         biasStr=num2str(N_bias(i));
%         currentStr=num2str(N_currents(j));
    end
    bias_axis=ones(length(N_currents),1)*N_bias(i)+bias_axis_add';
    plot(bias_axis,mu_tau_h_all((i-1)*length(N_currents)+1:(i-1)*length(N_currents)+length(N_currents)),'-o');%,'DisplayName',currentStr)
end
xlabel('bias [V] (+ tube current [5,15,25 mA])')
ylabel('mu.tau h [cm^2/V]')
box
grid
title(['mu.tau h vs bias:' sensor_name])
set(gca, 'YScale', 'log')
savefig([sensor_name '_mu.tau h vs bias.fig'])
saveas(gcf,[Output.sensor_name '_mu.tau h vs bias'], 'png')

% space charge vs bias
N_currents=unique(P(:,2));
N_bias=unique(P(:,1));
fig = figure;
hold on
for i=1:length(N_bias)
    for j=1:length(N_currents)
        bias_axis_add(j)=(j-2)*10;
%         biasStr=num2str(N_bias(i));
%         currentStr=num2str(N_currents(j));
    end
    bias_axis=ones(length(N_currents),1)*N_bias(i)+bias_axis_add';
    plot(bias_axis,rho_dark_all((i-1)*length(N_currents)+1:(i-1)*length(N_currents)+length(N_currents)),'-bo','DisplayName','dark');
    plot(bias_axis,rho_Xray_all((i-1)*length(N_currents)+1:(i-1)*length(N_currents)+length(N_currents)),'-rs','DisplayName','Xray');
    plot(bias_axis,rho_net_all((i-1)*length(N_currents)+1:(i-1)*length(N_currents)+length(N_currents)),'-g*','DisplayName','net');
end
xlabel('bias [V] (jittered with increasing tube current)')
ylabel('space charge density [e/cm^3]')
box
grid
title(['space charge vs bias:' sensor_name])
legend('dark','Xray','net')
% set(gca, 'YScale', 'log')
savefig([sensor_name '_space charge vs bias.fig'])

%% abs(space charge) vs bias
% N_currents=unique(P(:,2));
% N_bias=unique(P(:,1));
% fig = figure;
% hold on
% for i=1:length(N_bias)
%     for j=1:length(N_currents)
%         bias_axis_add(j)=(j-2)*10;
% %         biasStr=num2str(N_bias(i));
% %         currentStr=num2str(N_currents(j));
%     end
%     bias_axis=ones(length(N_currents),1)*N_bias(i)+bias_axis_add';
%     plot(bias_axis,abs(rho_dark_all((i-1)*4+1:(i-1)*4+length(N_currents))),'-bo','DisplayName','dark');
%     plot(bias_axis,abs(rho_Xray_all((i-1)*4+1:(i-1)*4+length(N_currents))),'-rs','DisplayName','Xray');
%     plot(bias_axis,abs(rho_net_all((i-1)*4+1:(i-1)*4+length(N_currents))),'-g*','DisplayName','net');
% end
% xlabel('bias [V] (+ tube current [5,15,25 mA])')
% ylabel('mean abs space charge density [e/cm^3]')
% box
% grid
% title(['mean-abs-space charge vs bias:' sensor_name])
% legend('dark','Xray','net')
% % set(gca, 'YScale', 'log')
% savefig([sensor_name '_abs_space charge vs bias.fig'])
%% save workspace
filename=[sensor_name '_all_space_charge_mu_tau_h.mat'];
save(filename,'rho_dark_all','rho_net_all','rho_Xray_all','mu_tau_h_all','sensor_name','N_bias','N_currents','Version')

%--save to excel
Voltage=P(:,1);
TubeCurrent=P(:,2);
outputtable=table(Voltage,TubeCurrent,rho_dark_all,rho_Xray_all,rho_net_all, mu_tau_h_all);
filename=[sensor_name '_summary.xlsx'];
writetable(outputtable, filename);
%--------modification history---------------------------------------------------
% 2021-11-02 created by Yuxin
% 2021-12-07 modified by Yuxin, version 1.1 save data summary in Excel file
% 2022-09-09 modified by Yuxin, version 1.2 modify the sectoin about
% checking fitting quality in 'Func_space_charge_line_profile.m'
% 2022-09-12 modified by Yuxin, version 1.3 Enabled to process data with
% only dark measurement.
% 2022-09-26 modified by Yuxin, version 2.0 (1) removed the dependence of mu.tau-h on X-ray tube current in Func_mu_tau_h_method2.m(2)not stop at broken data
% point 
% 2022-10-25 modified by Yuxin, version 3.0 
% (1)for calculating the dderivative of E-field, the new version uses smoothing by movingaverage instead of 5
% order polinominal fitting. most change happen to "Space_charge\Func_space_charge_line_profile.m"
% (2) remove 20 data points at both cathode and anode for calculating average space charge density
% (3) remove the mu_tau_h calculations and results based on the first cross-section method
% (4) if calculated mu_tau_h is wrongly high, it is given NaN (0.01 before)
% (5) figure location adjusted