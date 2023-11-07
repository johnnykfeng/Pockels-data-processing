%% -------------------------------------------------------------------------
% For group extraction of all Pockels profile for the same sensor with
% different bias and current conditions
%-------------------------------------------
% Version: 5.2
%-------------------------------------------
% Modification history is appended at the end
%-------------------------------------------
% Procedure:
% (1) Prepare a file called "P.txt". Each row is for one measurement
% condition. The first column is for bias, and the second is for X-ray tube
% current. The code will ask you if you have prepared this file and wait.
% Example:
% 500 0
% 500 25
% 700 0
% 700 25
% 900 0
% 900 25
%
% (2) Directories:
% Use the data folder as the "current folder" for ease of later processing.
% Add the folder of the code to the directory.
%
% (3) Run this code and follow the pop-up instructions.
%
% The output data is saved in the same folder as the raw data, not in
% the workspace.
%-------------------------------------------------------------------------
%%

clear
close all
prompt = {'Enter sensor name', 'Number of conditions to be extracted', 'Sensor thickness [mm]'};
dlgtitle = 'Test information ';
% Other inputs:
% Z_sensor = 2; %[mm] thickness of sensor. Default is 2mm
version = '5.2';
% Dimensions for the input dialog box
dims = [1, 35];
answer = inputdlg(prompt, dlgtitle, dims);
sensor_name = char(answer{1});
N_char = answer{2};
N = str2num(N_char);
N_char = answer{3};
Z_sensor = str2num(N_char);

waitfor(msgbox('Please prepare a table with the first column of "Bias" and the second of "X-ray tube current" (enter 0 if no X-ray) saved as "P.txt"'));
waitfor(msgbox('Select the "P.txt" file'));
[file_biased_crossed, path_biased_crossed] = uigetfile('*.txt'); % Reads the "P.txt" file
filename = [path_biased_crossed, file_biased_crossed];
parameter = importdata(filename);

%%
% These coordinates only make sense on a large 1920 x 1080 screen
f1 = figure('Name', 'Pockels average E-field profile per data');
movegui(f1, [40, 575]);
f2 = figure('Name', 'Pockels average E-field profile all data');
movegui(f2, [612, 575]);
hold on
f3 = figure('Name', 'Integrated E-field vs Bias');
movegui(f3, [40, 66]);
plot([0, 1000], [0, 1000], 'k--')
hold on
f4 = figure('Name', 'Pockels Image Data');
movegui(f4, [612, 66]);
%% Calibration
% read calibration_ini, it's does a lot of things
[Output, Calib, calib_status]=calibration_ini(parameter,sensor_name, Z_sensor);
%% process raw data of each measurement condition

Integral_Efield_all=zeros(N,1);
for i=1:N
    bias_string=num2str(parameter(i,1));
    flux_string=num2str(parameter(i,2));

    %%%%% Read the raw measurement data
    [Output]=Read_Pockels_data_group_extraction(sensor_name,parameter(i,1),parameter(i,2));
    %%%%% Calculates the E-field
    [Output, Calib]=Calculate_Efield(Output, Calib);

    figure(1)
    plot(Calib.x_all-Calib.x_all(Calib.cathode),Output.E_cross_section_average_corrected,'displayname',[bias_string 'V,' flux_string 'mA']);
    title(['Pockels average profile sensor-' Output.sensor_name '@' bias_string 'V,' flux_string 'mA'])
    xlabel('Cathode to anode distance (mm)')
    ylabel('E field (V/m)')

    y = Output.E_cross_section_average_corrected(round(Calib.cathode): round(Calib.anode));
    y(isnan(y)) = 0;
    Output.integral_Efield = trapz(Calib.x_sensor*1e-3, y);
    Integral_Efield_all(i) = Output.integral_Efield;
    %     end
    figure(4)
    imagesc(Output.E_field_biased_corrected)
    axes1=gca;
    set(axes1,'DataAspectRatio',[1 1 1],'Layer','top');
    set(axes1,'View',[-90 90]);
    colorbar
    box
    
    if Output.flux > 0 % for non dark current
        figure(2)
        plot(Calib.x_all-Calib.x_all(Calib.cathode),Output.E_cross_section_average_corrected*1e-5,'displayname',[bias_string 'V,' flux_string 'mA']);
        title(['Pockels average profile sensor-' Output.sensor_name])
        xlabel('Cathode to anode distance (mm)')
        ylabel('E field (100V/mm)')
        box

        figure(3)
        plot(Output.bias,Output.integral_Efield,'o','displayname',[bias_string 'V,' flux_string 'mA']);
        title(['Integral E-field sensor-' Output.sensor_name])
        xlabel('Bias (V)')
        ylabel('Integral E-field (V)')
        box

        figure(4)
        imagesc(Output.E_field_biased_corrected)
        E_ave=parameter(i,1)/Z_sensor*1E3;
        E_max=E_ave*1.7;
        clim([0,E_max])
        title(['Corrected E-field for sensor-' Output.sensor_name ' @ ' bias_string 'V, ' flux_string 'mA'])
        xlabel('x-pixel')
        ylabel('y-pixel')
        axes1=gca;
        set(axes1,'DataAspectRatio',[1 1 1],'Layer','top');
        set(axes1,'View',[-90 90]);
        colorbar
        box
    else  % Dark current plots
        figure(2)
        plot(Calib.x_all-Calib.x_all(Calib.cathode), Output.E_cross_section_average_corrected*1e-5,'--', 'displayname', [bias_string 'V,' flux_string 'mA']);
        title(['Pockels average profile sensor-' Output.sensor_name])
        xlabel('thickness (mm)')
        ylabel('E field (100V/mm)')
        box

        figure(3)
        plot(Output.bias,Output.integral_Efield,'s','displayname',[bias_string 'V,' flux_string 'mA']);
        title(['Integral E-field sensor-' Output.sensor_name])
        xlabel('Bias (V)')
        ylabel('Integral E-field (V)')
        box

        figure(4)
        imagesc(Output.E_field_biased_corrected)
        E_ave=parameter(i,1)/Z_sensor*1E3;
        E_max=E_ave*1.7;
        clim([0,E_max])
        title(['Corrected E-field for sensor-' Output.sensor_name ' @ ' bias_string 'V, ' flux_string 'mA'])
        xlabel('x-pixel')
        ylabel('y-pixel')
        axes1=gca;
        set(axes1,'DataAspectRatio',[1 1 1],'Layer','top');
        set(axes1,'View',[-90 90]);
        colorbar
        box
    end

    %Saves the output files into a directory specified by user
    %--------------
    Output.version = version;
    %--------------
    Output.thickness = Z_sensor;
    Dir = path_biased_crossed;
    save([Dir '\' Output.sensor_name '_' bias_string 'V_' flux_string 'mA_Pockels_output.mat'], 'Output');
    figure(4)
    savefig([Output.sensor_name '_' bias_string 'V_' flux_string 'mA_E-field_uncorrected'])
    saveas(gcf,[Output.sensor_name '_' bias_string 'V_' flux_string 'mA_E-field_uncorrected'], 'png')
    if strcmp(calib_status,'No')
        save([Dir '\' Output.sensor_name '_Pockels_calib_file.mat'], 'Calib');
    end
    calib_status='yes';
end

%% --judge integrated E-field-----------
if N>2
    bias_all=unique(parameter(:,1));
    Int_E_dark=Integral_Efield_all(parameter(:,2)==0);

    bias_500=bias_all(bias_all>=500);
    Int_E_dark_500=Int_E_dark(bias_all>=500);

    fitResults = polyfit(bias_500, Int_E_dark_500, 1);
    slope = fitResults(1)
    deviation = (slope-1);
    if abs(deviation)>0.3
        fprintf(2,'Warning:the integral of E-field poorly matches bias, a better calibration may be needed');
    end
else
    deviation=0;
end
save([Dir '\' Output.sensor_name '_integral_E_field.mat'], 'Integral_Efield_all');

%% --------------------------------------

figure(2)
savefig([Output.sensor_name '_all_E_field_profiles'])
saveas(gcf,[Output.sensor_name '_all_E_field_profiles'], 'png')
figure(3)
savefig([Output.sensor_name '_Integral_E_field_vs_bias'])
saveas(gcf,[Output.sensor_name '_Integral_E_field_vs_bias'], 'png')
figure(1)
hold off
figure(2)
hold off
figure(3)
hold off
if abs(deviation)>0.3
    waitfor(msgbox(['Process is complete. The integral of E-field poorly matches bias, a better calibration may be needed']));
else
waitfor(msgbox(['Process is complete']));
end
clearvars -except Output f1 f2 f3 f4 Calib slope

%% Revision history
%2023-08-23 Version 5.2 by Yuxin. Corrected the parameters in the Pockels equation and the calculation of r41, in calibration_ini.m 
%2023-01-27 Version 5.1 by Yuxin. the colrmap range in Fig.4 E-field map is
% made adjusted automatically based on the bias.
%2023-01-03 Version 5.0 by Yuxin. image distortion is corrected for
%calculated E-field. by function:Func_find_sensor_edges_and_distortion_correction.m inside Calculate_Efield.m
%(I_over_I0 is not stored in Output anymore.)
%2022-11-07 Version 4.0 by Yuxin. allows negative E-field in the calculation 
% to tolarent noise in low field region and thus removes the fake above zero 
% E-field in the zero field region. Modifications are in "Calculate_Efield.m"
%2022-09-26 Version 3.1 by Yuxin. Add lines to skip broken data.
%2022-09-08 Version 3.0 by Yuxin. Automatically load all measurements.
%Modifications impelmented in Read_Pockels_data_group_extraction.m and add
%the calibration voltage into Calib.mat file
%2022-06-15 Version 2.6 by Yuxin. In calibration, median of series of
%integers can result in decimal numbers, which is corrected.
%2022-04-04 Version 2.5 by Yuxin. (1)Save figures also in .png format
%      (2)modifiy the E-field profile figures to be: cathode at 0, vertical
%      axis in unit of (100V/mm)
%2022-03-11 Version 2.4 by Yuxin. For "calibration_ini.m", added a step to
% manually select a region in the "0V parallel" image to remove apparent regions
% with light through other than the sensor. thus, the chance to wrongly find
% sensor edge is reduced.
%2021-12-09 Version 2.3 by Yuxin. Change the sensor thickness as an input parameter
%2021-09-08 Version 2.2 by Yuxin. Improvement of judging where the signal is
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