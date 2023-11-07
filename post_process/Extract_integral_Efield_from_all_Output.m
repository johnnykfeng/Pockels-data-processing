%-------------------------------------------------------------------------
% For extracting "integral_Efield" from all Output files within one folder
% 2021-11-08 created by Yuxin
%-------------------------------------------------------------------------

clear
load 'P.txt'
N=length(P);


%% load calibration file and sensor name
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

%% read out integral_Efield from Output files

quantity='integral_Efield';
for i=1:N
    bias_str=num2str(P(i,1));
    current_str=num2str(P(i,2));
    [integral_Efield(i)]=Func_extract_from_output(P(i,1),P(i,2),quantity);

end

% figure
plot(P(:,2),integral_Efield)
xlabel('X-ray tube current [mA])')
ylabel('Integration of E-field [V]')
box
grid
title(['Integration of E-field vs tube current:' sensor_name])