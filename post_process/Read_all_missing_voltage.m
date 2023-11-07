%-------------------------------------------------------------------------
% For calculating the "missing voltage" for all
% measurements within the same folder
%-------------------------------------------------------------------------
% 2022-09-14 created by Yuxin
% 2023-06-12 improved figure plot and result output by Yuxin
%-------------------------------------------------------------------------
addpath 'C:\Users\10404\OneDrive - Redlen\Documents\MATLAB'
% addpath('C:\Users\10404\OneDrive - Redlen\Documents\MATLAB\Pockels\post_process\Space_charge')
addpath('C:\Users\10404\OneDrive - Redlen\Documents\MATLAB\Pockels\post_process')
addpath('C:\Users\10404\OneDrive - Redlen\Documents\MATLAB\Pockels')
%%
clear
load 'P.txt'
N=length(P);
Int_E_dark=NaN(N,1);
Int_E_Xray=NaN(N,1);

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
    Bias=num2str(P(i,1));
    Current=num2str(P(i,2));
    string1=['*' Bias 'V_' Current 'mA' '*output.mat'];
    S = dir(fullfile(string1));
    if size(S)>1
        fprintf(2,'Warning:there are more than 1 data set with the same bias and current');
    end
    load(S.name);
    if P(i,2)==0
        Int_E_dark(i)=Output.integral_Efield;
        Int_E_dark_temp=Int_E_dark(i);
        bias_dark_temp=P(i,1);
    else
        Int_E_Xray(i)=Output.integral_Efield;
        if P(i,1)==bias_dark_temp
            Int_E_dark(i)=Int_E_dark_temp;
        else
            fprintf(2,'Warning:cannot find data under dark');
        end

    end
end
Missing_voltage=Int_E_dark-Int_E_Xray;

%% plot figures
Voltage=P(:,1);
TubeCurrent=P(:,2);
I_uni=unique(TubeCurrent);
M=length(I_uni);
% colormap = jet(numel(I_uni)); % Using 'jet' colormap
colormap = hsv(numel(I_uni)); % Using 'jet' colormap
figure; hold on; 
for i = 1:N
    j=find(I_uni==TubeCurrent(i));
    plot(Voltage(i), Missing_voltage(i), 'o', 'Color', colormap(j, :));
end
xlabel('bias [V]')
ylabel('missing voltage [V]')
box
grid
title(['missing voltage:' sensor_name])
hold off;


%% save workspace
filename=[sensor_name '_missing_voltage.mat'];
save(filename,'Int_E_dark','Int_E_Xray','Missing_voltage')

%--save to excel
Voltage=P(:,1);
TubeCurrent=P(:,2);
outputtable=table(Voltage,TubeCurrent,Int_E_dark,Int_E_Xray,Missing_voltage);
filename=[sensor_name '_missing_voltage.xlsx'];
writetable(outputtable, filename);