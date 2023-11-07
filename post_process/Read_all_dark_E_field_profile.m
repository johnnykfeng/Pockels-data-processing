%-------------------------------------------------------------------------
% For extracting all the E-field profiles under dark within one folder

% 2022-10-31 created by Yuxin
%-------------------------------------------------------------------------

load 'P.txt'
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
figure, hold on
for i=1:N
    Bias=num2str(bias(i));
    Current='0';
    string1=['*' Bias 'V_' Current 'mA' '*output.mat'];
    S = dir(fullfile(string1));
    if size(S)>1
        fprintf(2,'Warning:there are more than 1 data set with the same bias and current');
    end
    load(S.name);
    plot(Calib.x_all-Calib.x_all(Calib.cathode),Output.E_cross_section_average_corrected*1e-5,'displayname',[Bias 'V']);

end
zeromark=zeros(1,length(Calib.x_all));
plot(Calib.x_all-Calib.x_all(Calib.cathode),zeromark,'k--','displayname','zero');
title(['Dark E-field profiles, sensor-' Output.sensor_name])
xlabel('thickness (mm)')
ylabel('E field (100V/mm)')
box
xlim([0,2])

%% save workspace
% filename=[sensor_name '_missing_voltage.mat'];
% save(filename,'Int_E_dark','Int_E_Xray','Missing_voltage')
% 
% %--save to excel
% Voltage=P(:,1);
% TubeCurrent=P(:,2);
% outputtable=table(Voltage,TubeCurrent,Int_E_dark,Int_E_Xray,Missing_voltage);
% filename=[sensor_name '_missing_voltage.xlsx'];
% writetable(outputtable, filename);