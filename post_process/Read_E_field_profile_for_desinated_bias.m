%-------------------------------------------------------------------------
% For extracting all the E-field profiles under dark within one folder

% 2023-05-03 created by Yuxin
%-------------------------------------------------------------------------
bias_aim=700;
% current_aim=[0 5 25];
%-------------------------------------------------------------------------
load 'P.txt'
bias=P(:,1);
current=P(:,2);
Bias=num2str(P(1,1));
Current=num2str(P(1,2));
index=find(bias==bias_aim);
%%
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
% bias=unique(P(:,1));
N=length(index);
% figure, hold on
for i=1:N
    Bias=num2str(bias(index(i)));
    Current=num2str(current(index(i)));
    string1=['*' Bias 'V_' Current 'mA' '*output.mat'];
    S = dir(fullfile(string1));
    if size(S)>1
        fprintf(2,'Warning:there are more than 1 data set with the same bias and current');
    end
    load(S.name);
    if current(index(i))==0
        plot(Calib.x_all-Calib.x_all(Calib.cathode),Output.E_cross_section_average_corrected*1e-5,'--','displayname',[Bias 'V ' Current 'mA']);
    else
        plot(Calib.x_all-Calib.x_all(Calib.cathode),Output.E_cross_section_average_corrected*1e-5,'displayname',[Bias 'V ' Current 'mA']);
    end
    hold on
end
% zeromark=zeros(1,length(Calib.x_all));
% plot(Calib.x_all-Calib.x_all(Calib.cathode),zeromark,'k--','displayname','zero');
title(['E-field profiles, sensor-' Output.sensor_name])
xlabel('thickness (mm)')
ylabel('E field (100V/mm)')
legend
box
xlim([0,2])