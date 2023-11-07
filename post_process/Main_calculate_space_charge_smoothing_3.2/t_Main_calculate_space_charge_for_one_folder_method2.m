%-------------------------------------------------------------------------
% for time dependent Pockels measurement at the same X-ray flux
% For calculating the space charge profile and mu.tau-h for all
% measurements within the same folder 
%-------------------------------------------------------------------------
% Step0: prepare t.txt file with all the measurement time

% Step1: Check whether the P.txt file exists with all biases and currents
% Step2: set the working directory to the folder with Pockels raw data
% Step3: run this code, choose "add to path" if directory selection window
% pops up

% 2022-09-06 created by Yuxin based on
% Main_calculate_space_charge_for_one_folder_method2.m 
% 2021-12-07 modified by Yuxin, version 1.1 save data summary in Excel file

%-------------------------------------------------------------------------
addpath 'C:\Users\10404\OneDrive - Redlen\Documents\MATLAB'
addpath('C:\Users\10404\OneDrive - Redlen\Documents\MATLAB\Pockels\post_process\Space_charge')
addpath('C:\Users\10404\OneDrive - Redlen\Documents\MATLAB\Pockels\post_process')
addpath('C:\Users\10404\OneDrive - Redlen\Documents\MATLAB\Pockels')
%%
clear
%---------input-------
tube_current=25;%[mA] tube current is defined here instead of the P.txt file
timeUnit='run';% unit of time
%---------------------
load 'P.txt'
load 't.txt'
N=length(t);
rho_dark_all=NaN(N,1);
rho_Xray_all=NaN(N,1);
rho_net_all=NaN(N,1);
mu_tau_h_all=NaN(N,1);
mu_tau_h_all2=NaN(N,1);
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
        [rho_net,rho_net_mean,rho_dark, rho_Xray,mu_tau_h,mu_tau_h2,errorcode,coef_dark,coef_Xray,xx]=Func_space_charge_profiles_time(P(i,1),P(i,2),tube_current,t(i),timeUnit);
        rho_dark_all(i)=mean(rho_dark);
        if P(i-1,2)==0
            rho_dark_all(i-1)=mean(rho_dark);
        end
        rho_Xray_all(i)=mean(rho_Xray);
        rho_net_all(i)=mean(rho_net);
        mu_tau_h_all(i)=mu_tau_h;
        mu_tau_h_all2(i)=mu_tau_h2;
    end
end

%% plot figures
% % mu.tau-h vs time
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
%     plot(bias_axis,mu_tau_h_all((i-1)*length(N_currents)+1:(i-1)*length(N_currents)+length(N_currents)),'-o');%,'DisplayName',currentStr)
% end
% xlabel('bias [V] (+ tube current [5,15,25 mA])')
% ylabel('mu.tau h [cm^2/V]')
% box
% grid
% title(['mu.tau h vs bias:' sensor_name])
% set(gca, 'YScale', 'log')
% savefig([sensor_name '_mu.tau h vs bias.fig'])
% saveas(gcf,[Output.sensor_name '_mu.tau h vs bias'], 'png')
% mu.tau-h2 vs bias
fig = figure;
plot(t,mu_tau_h_all2,'-o')
xlabel(['time [' timeUnit ']'])
ylabel('mu.tau h [cm^2/V]')
box
grid
title(['mu.tau h (method2) vs time:' sensor_name])
% set(gca, 'YScale', 'log')
savefig([sensor_name '_mu.tau h vs time.fig'])

% space charge vs time
% N_currents=unique(P(:,2));
% N_bias=unique(P(:,1));

fig = figure;
hold on
plot(t,rho_dark_all,'-bo','DisplayName','dark');
plot(t,rho_Xray_all,'-rs','DisplayName','Xray');
plot(t,rho_net_all,'-g*','DisplayName','net');
xlabel(['time [' timeUnit ']'])
ylabel('space charge density [e/cm^3]')
box
grid
title(['space charge vs time:' sensor_name])
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
save(filename,'rho_dark_all','rho_net_all','rho_Xray_all','mu_tau_h_all','mu_tau_h_all2','sensor_name','t','timeUnit')

%--save to excel
Voltage=P(:,1);
TubeCurrent=ones(length(t),1)*tube_current;
% TubeCurrent=P(:,2);
TimeColumn=t;
outputtable=table(TimeColumn,Voltage,TubeCurrent,rho_dark_all,rho_Xray_all,rho_net_all,mu_tau_h_all2);
filename=[sensor_name '_summary.xlsx'];
writetable(outputtable, filename);