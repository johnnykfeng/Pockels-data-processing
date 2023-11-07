%-------------------------------------------------------------------------
% For extracting quantities from the Output file of Pockels data, for group
% extraction of one quantity for many measurement conditions in a loop
% Example: Output.integral_Efield
% 2021-11-08 created by Yuxin

% List of quantities:
% 'sensor_name'
% 'I_crossed_biased_field'
% 'I_over_I0'
% 'E_field_biased_corrected'
% 'E_cross_section_average_corrected'
% 'integral_Efield'
% 'version'
% 'thickness'
%-------------------------------------------------------------------------

function [output_quantity]=Func_extract_from_output(bias,current,quantity)
%load data with selected bias and X-ray tube current into workspace
Bias=num2str(bias);
Current=num2str(current);
string1=['*' Bias 'V_' Current 'mA' '*output.mat'];
S = dir(fullfile(string1));
load(S.name);
% Output_Xray=Output;
% 
% string2=['*' Bias 'V_0mA' '*.mat'];
% S = dir(fullfile(string2));
% load(S.name);
% Output_dark=Output;

% string3=['*calib*.mat'];
% S = dir(fullfile(string3));
% load(S.name);

quantity_name=['Output.' quantity];
output_quantity=eval(quantity_name);

end