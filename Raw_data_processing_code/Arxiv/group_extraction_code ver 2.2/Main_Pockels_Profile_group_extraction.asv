%-------------------------------------------------------------------------
% For group extraction of all Pockels profile for the same sensor with
% different bias and current conditons
%-------------------------------------------
%Version:2.2
%-------------------------------------------
% modification history is appended in the end 
%-------------------------------------------
%Procedure:
%(1)Prepare a file called "P.txt". Each row is for one measurement
%condition. the first column is for bias and the second for X-ray tube
%current. the code will ask you if you have prepared this file and wait.
%Example"
% 500 0
% 500 25
% 700 0
% 700 25
% 900 0
% 900 25

%(2) directories
%use the data folder as the "current folder" for ease of later process
%add the folder of the code to directory
%(3)run this code and follow the pop-up instructions

%The output data are saved in the same folder of the raw data, not in
%Workspace
%-------------------------------------------------------------------------
%%
clear
close all
prompt = {'Enter sensor name','No. of conditions to be extracted'};
dlgtitle = 'Test information ';
%---other inputs----
Z_sensor=2;%[mm] thickness of sensor. default is 2mm
version='2.2';
%-------------
dims = [1 35];
answer = inputdlg(prompt,dlgtitle,dims);
sensor_name=char(answer{1,1});
N_char=answer{2,1};
N=str2num(N_char);

waitfor(msgbox(['please prepare a table with the first column of "Bias" and the second of "Xray tube current" (enter 0 if no Xray) saved as "P.txt"']))
waitfor(msgbox(['Select the P.txt file']));
[file_biased_crossed,path_biased_crossed] = uigetfile('*.txt');
filename=([path_biased_crossed  file_biased_crossed]);
parameter=importdata(filename);

%%
f1=figure(1);
movegui(f1,[40 575]);%movegui(f1,[40 903]);
f2=figure(2);
movegui(f2,[612 575]);%movegui(f2,[612 903]);
hold on
f3=figure(3);
movegui(f3,[40 66]);%movegui(f3,[40 392]);
plot([0 900],[0 900],'k--')
hold on
f4=figure(4);
movegui(f4,[612  66]);%movegui(f4,[612  392]);
% f5=figure(5);
% movegui(f5,[1184 575]);%movegui(f2,[612 903]);
% hold on
%% Calibration
[Output,Calib,calib_status]=calibration_ini(parameter,sensor_name,Z_sensor);
%% process raw data of each measurement condition
for i=1:N
    bias_string=num2str(parameter(i,1));
    flux_string=num2str(parameter(i,2));
    
    %%%%% Read the raw measurement data
    [Output]=Read_Pockels_data_group_extraction(sensor_name,parameter(i,1),parameter(i,2));
    %%%%% Calculates the E-field
    [Output]=Calculate_Efield(Output, Calib);
    
    figure(1)
    plot(Calib.x_all,Output.E_cross_section_average_corrected,'displayname',[bias_string 'V,' flux_string 'mA']);
    title(['Pockels average profile sensor-' Output.sensor_name '@' bias_string 'V,' flux_string 'mA'])
    xlabel('thickness (mm)')
    ylabel('E field (V/m)')
    %---------added-----------------
    hold on
    E_max=(1/Calib.alpha) * (asin(1));% max possible E-field calculated from light intensity
    L=length(Output.E_cross_section_average_corrected);
    E_max_vec=ones(1,L)*E_max;
    plot(Calib.x_all,E_max_vec,'k--','displayname','E-max')
    hold off
    %-----------------------------
if Output.bias>600
    sin_max = questdlg('Does the graph hit the maximum of the sine wave?');
    switch sin_max
        case 'Yes'
            [x_max,y_max]=ginput(1);
            %             close(f1)
            %%%%%%Correcting for sin
            Nom3(:,:,1)=(Output.I_crossed_biased_field(:,:,1)-Calib.I_crossed_nobias_field);
            A3=Nom3(:,:,1)<0;
            Nom3(A3)=0;
            I_over_I0=(Nom3(:,:,1)./Calib.I_parallel_nobias_field);
            I_over_I0(I_over_I0>1)=1;
            E_biased_corrected_recons(:,round(x_max):520)=(1/Calib.alpha) *(max(mean(asin(sqrt(I_over_I0))))*2-asin(sqrt(I_over_I0(:,round(x_max):end,1))));
            
            Output.E_cross_section_average_corrected(:,round(x_max):end,1)=mean(E_biased_corrected_recons(:,round(x_max):end,1),'omitnan');
            
            figure
            plot([1:520]*9.3e-3,Output.E_cross_section_average_corrected(:,:,1));
            title(['CORRECTED-Pockels average profile sensor-' Output.sensor_name ', bias=' Output.bias, ', Xray tube current=' Output.flux])
            xlabel('thickness (mm)')
            ylabel('E field (V/m)')
            
            %Integrate E-field over thickness between chosen cathode and anode
            y=Output.E_cross_section_average_corrected(round(Calib.cathode):round(Calib.anode));
            Output.integral_Efield=trapz(Calib.x_sensor*1e-3,y);
        case 'No'
            %Integrate E-field over thickness between chosen cathode and anode
            y=Output.E_cross_section_average_corrected(round(Calib.cathode):round(Calib.anode));
            Output.integral_Efield=trapz(Calib.x_sensor*1e-3,y);
            
    end
else
    y=Output.E_cross_section_average_corrected(round(Calib.cathode):round(Calib.anode));
    Output.integral_Efield=trapz(Calib.x_sensor*1e-3,y);
end 
    figure(4)
    imagesc(Output.E_field_biased_corrected)
    axes1=gca;
    set(axes1,'DataAspectRatio',[1 1 1],'Layer','top');
    set(axes1,'View',[-90 90]);
    colorbar
    box
    if Output.flux>0
        figure(2)
        plot(Calib.x_all,Output.E_cross_section_average_corrected,'displayname',[bias_string 'V_' flux_string 'mA']);
        title(['Pockels average profile sensor-' Output.sensor_name])
        xlabel('thickness (mm)')
        ylabel('E field (V/m)')
        box
        
        figure(3)
        plot(Output.bias,Output.integral_Efield,'o','displayname',[bias_string 'V,' flux_string 'mA']);
        title(['Integral E-field sensor-' Output.sensor_name])
        xlabel('Bias (V)')
        ylabel('Integral E-field (V)')
        box
        figure(4)
        imagesc(Output.E_field_biased_corrected)
        title(['E-field sensor-' Output.sensor_name '@' bias_string 'V,' flux_string 'mA'])
        axes1=gca;
        set(axes1,'DataAspectRatio',[1 1 1],'Layer','top');
        set(axes1,'View',[-90 90]);
        colorbar
        box    
    else
        figure(2)
        plot(Calib.x_all,Output.E_cross_section_average_corrected,'--','displayname',[bias_string 'V_' flux_string 'mA']);
        title(['Pockels average profile sensor-' Output.sensor_name])
        xlabel('thickness (mm)')
        ylabel('E field (V/m)')
        box
        
        figure(3)
        plot(Output.bias,Output.integral_Efield,'s','displayname',[bias_string 'V,' flux_string 'mA']);
        title(['Integral E-field sensor-' Output.sensor_name])
        xlabel('Bias (V)')
        ylabel('Integral E-field (V)')
        box
        figure(4)
        imagesc(Output.E_field_biased_corrected)
        title(['E-field sensor-' Output.sensor_name '@' bias_string 'V,' flux_string 'mA'])
        axes1=gca;
        set(axes1,'DataAspectRatio',[1 1 1],'Layer','top');
        set(axes1,'View',[-90 90]);
        colorbar
        box
    end
    
    %Saves the output files into a directory specified by user
    %--------------
    Output.version=version;
    %--------------
    Output.thickness=Z_sensor;
    Dir=path_biased_crossed;
    save([Dir '\' Output.sensor_name '_' bias_string 'V_' flux_string 'mA_Pockels_output.mat'], 'Output');
    figure(4)
    savefig([Output.sensor_name '_' bias_string 'V_' flux_string 'mA_E-field_uncorrected'])
    if strcmp(calib_status,'No')
        save([Dir '\' Output.sensor_name '_Pockels_calib_file.mat'], 'Calib');
    end
    calib_status='yes';
end
figure(2)
savefig([Output.sensor_name '_all_E_field_profiles'])
figure(3)
savefig([Output.sensor_name '_Integral_E_field_vs_bias'])
figure(1)
hold off
figure(2)
hold off
figure(3)
hold off
waitfor(msgbox(['process is complete']));
clearvars -except Output f1 f2 f3 f4 Calib

%% Revision history
%2021-09-08 Version 2.2 by Yuxin. Improvement of juding where the signal is
%in "calibration_ini.m" (line 78-87)
%2021-08-09 Version 2.1 by Yuxin. Minor improvements for calibration data
%selection and simplified P.txt file preparation. Corrected some bugs
%2021-07-27 Version 2.0 by Yuxin
    % *redefined calibration process in the subprogram of
    % *calibration_ini.m,including
    % (1) automatic determining sensor edges
    % (2) the thickness of sensor is fixed at 2mm (can be modified in the
    % main code) instead of being measured by nuber of camera pixels
    % (3) the 'alpha' (and subsequently, the 'r_41') parameters are
    % calibrated by forcing the integration of E-field equal to the applied
    % bias, in the seleted region for calibration
    % *improved some part of the codes to handel more general situations
%2021-07-05 Version 1.2 Yuxin modifed (1)fix bugs (2) adjust the E-field image
%2021-06-29 Version 1.1 Yuxin modified (1)fix bugs (2) include the 2D E-field in
%2021-06-21 Version 1.0 Created by Yuxin Song based on Niloofar's single profile code