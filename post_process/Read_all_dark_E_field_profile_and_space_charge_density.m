%-------------------------------------------------------------------------
% For extracting all the E-field profiles under dark within one folder and
% manually select a linear region for calculating rho_dark

% 2023-6-23 created by Yuxin
%-------------------------------------------------------------------------
load 'P.txt'
close all
Bias=num2str(P(1,1));
Current=num2str(P(1,2));
string1=['*' Bias 'V_' Current 'mA' '*output.mat'];
S = dir(fullfile(string1));
if size(S)>1
    fprintf(2,'Warning:there are more than 1 data set with the same bias and current');
end
load(S.name);
sensor_name=Output.sensor_name;
string3=['*calib*.mat'];
S = dir(fullfile(string3));
load(S.name);

%%
bias=unique(P(:,1));
N=length(bias);
f51=figure(51);
hold on
movegui(f51,[40 575]);%movegui(f1,[40 903]);
f52=figure(52);
movegui(f52,[612 575]);%movegui(f1,[40 903]);
rho_dark=zeros(length(bias),1);
constants
eps_CZT=C.eps_0*C.epsilon;
for i=1:N
    Bias=num2str(bias(i));
    Current='0';
    string1=['*' Bias 'V_' Current 'mA' '*output.mat'];
    S = dir(fullfile(string1));
    if size(S)>1
        fprintf(2,'Warning:there are more than 1 data set with the same bias and current');
    end
    load(S.name);
    figure(52)
    plot(Output.E_cross_section_average_corrected*1e-5,'displayname',[Bias 'V']);
    title(['Dark E-field profiles, sensor-' Output.sensor_name])
    xlabel('element number')
    ylabel('E field (V/m)')
    box
    if i==1
        waitfor(msgbox(['Select the region for linear fitting']));
    end
    [points,y]=ginput(2);
    points=round(points);
    E_selected=Output.E_cross_section_average_corrected(points(1):points(2));%[V/m]
    x_selected=Calib.x_all(points(1):points(2))-Calib.x_all(Calib.cathode);%[mm]
    p = polyfit(x_selected, E_selected, 1);
    E_fit = polyval(p,x_selected);
    rho_dark(i)=-p(1).*eps_CZT./10^3/C.q; %[e/cm^3] space charge density 

    figure(51)
    plot(Calib.x_all-Calib.x_all(Calib.cathode),Output.E_cross_section_average_corrected,'displayname',[Bias 'V']);
    plot(x_selected,E_fit,'--k','displayname',[Bias 'V fit'])
end
figure(51)
zeromark=zeros(1,length(Calib.x_all));
plot(Calib.x_all-Calib.x_all(Calib.cathode),zeromark,'k--','displayname','zero');
title(['Dark E-field profiles, sensor-' Output.sensor_name])
xlabel('thickness (mm)')
ylabel('E field (V/m)')
box
xlim([0,2])
savefig([Output.sensor_name '_all_E_field_under_dark_profiles'])
saveas(gcf,[Output.sensor_name '_all_E_field_under_dark_profiles'], 'png')
%%
f53=figure(53);
movegui(f53,[40 66]);
plot(bias,rho_dark,'-o','DisplayName','sensor_name')
title(['\rho_d_a_r_k, sensor-' Output.sensor_name])
xlabel('Bias (V)')
ylabel('\rho_d_a_r_k (1/cm^3)')
grid
box
savefig([Output.sensor_name '_rho_dark_vs_bias'])
saveas(gcf,[Output.sensor_name '_rho_dark_vs_bias'], 'png')
%% save workspace
filename=[sensor_name '_all_rho_dark.mat'];
save(filename,'bias','rho_dark')

%--save to excel
outputtable=table(bias,rho_dark);
filename=[sensor_name '_rho_dark_vs_bias.xlsx'];
writetable(outputtable, filename);