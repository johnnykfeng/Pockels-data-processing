%-------------------------------------------------------------------------
%for grou extraction of all Pockels profile for the same sensor with
%different bias and current conditons
%2021-06-21 Created by Yuxin Song based on Niloofar's single profile code
%2021-06-29 Yuxin modified (1)fix bugs (2) include the 2D E-field in
%Output(3) automatically save key figures
%2021-07-05 Yuxin modifed (1)fix bugs (2) adjust the E-field image
%Ver:1.2
%
%Procedure:
%(1)Prepare a file called "P.txt". Each row is for one measurement
%condition. the first column is for bias and the second for X-ray tube
%current. the code will ask you if you have prepared this file and wait.
%Example"
% 700 0
% 700 25
% 500 0
% 500 25
% 600 0
% 600 25
% 800 0
% 800 25
% 900 0
% 900 25
% the first row will also be used for the common calibration for all
% measurements, so it is better to use 700V 0mA data
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

waitfor(msgbox(['Select the 0V, parallel data']));
[file_0V_parallel,path_0V_parallel] = uigetfile('*.txt');
Output.data_directory=path_0V_parallel;

buffer1 = fileread(char([path_0V_parallel  file_0V_parallel]));
data1 = textscan(buffer1,'%f %f %f', 'headerlines', 19) ;
I_parallel_nobias_matrix(:,1)=data1{1,1};
I_parallel_nobias_matrix(:,2)=data1{1,2};
I_parallel_nobias_matrix(:,3)=data1{1,3};
%[I_parallel_nobias_matrix(:,1) I_parallel_nobias_matrix(:,2) I_parallel_nobias_matrix(:,3)] = textread([path_0V_parallel  file_0V_parallel], '%f %f %f', 'headerlines', 19) ;

waitfor(msgbox(['Select the 0V, crossed data']));
[file_0V_crossed,path_0V_crossed] = uigetfile('*.txt');

buffer2 = fileread(char([path_0V_crossed  file_0V_crossed]));
data2 = textscan(buffer2,'%f %f %f', 'headerlines', 19) ;
I_crossed_nobias_matrix(:,1)=data2{1,1};
I_crossed_nobias_matrix(:,2)=data2{1,2};
I_crossed_nobias_matrix(:,3)=data2{1,3};
%[I_crossed_nobias_matrix(:,1) I_crossed_nobias_matrix(:,2) I_crossed_nobias_matrix(:,3)] = textread([path_0V_crossed  file_0V_crossed], '%f %f %f', 'headerlines', 19) ;

%%
calib_status = questdlg('Do you want to use existing calibration file?');
if strcmp(calib_status,'Yes')
    waitfor(msgbox(['Select your calibration data']));
    [Calib_file, Calib_path] = uigetfile('*.mat');
    load([Calib_path '\' Calib_file]);
    
else if strcmp(calib_status,'No')
        Calib.n_0=2.8;
        Calib.r_41=3.8*10^-12;
        Calib.d=11.175*10^-3;
        Calib.lambda_0=980*10^-9;
        Calib.alpha=(sqrt(3)*pi*Calib.n_0^3*Calib.r_41*Calib.d)/(2*Calib.lambda_0);
    end
end
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
for i=1:N
    bias_string=num2str(parameter(i,1));
    flux_string=num2str(parameter(i,2));
    %%%%%Gets the intensity profiles for 0v parallel, 0v crossed and biased crossed
    [Output]=Read_Pockels_data_group_extraction(sensor_name,parameter(i,1),parameter(i,2),I_parallel_nobias_matrix,I_crossed_nobias_matrix);
    %%%%%Calculates the E-field
    [Output]=Calculate_Efield(Output, Calib.alpha);
    figure(1)
    plot([1:520],Output.E_cross_section_average_corrected(:,:,1),'displayname',[bias_string 'V,' flux_string 'mA']);
    title(['Pockels average profile sensor-' Output.sensor_name '@' bias_string 'V,' flux_string 'mA'])
    xlabel('thickness (mm)')
    ylabel('E field (V/m)')
    
    %%%%%%%%%%%%%%%%%%%%%%
    %    Choosing sensor thickness
    %Each camera pixel = 9.3um
    if strcmp(calib_status,'No')
        waitfor(msgbox(['Select the cathode and anode ']));
        [thickness,y]=ginput(2);
        Calib.cathode=thickness(1);
        Calib.anode=thickness(2);
    end
    
    
    
    sin_max = questdlg('Does the graph hit the maximum of the sine wave?');
    switch sin_max
        case 'Yes'
            [x_max,y_max]=ginput(1);
            %             close(f1)
            %%%%%%Correcting for sin
            Nom3(:,:,1)=(Output.I_crossed_biased_field(:,:,1)-Output.I_crossed_nobias_field);
            A3=Nom3(:,:,1)<0;
            Nom3(A3)=0;
            I_over_I0=(Nom3(:,:,1)./Output.I_parallel_nobias_field);
            I_over_I0(I_over_I0>1)=1;
            E_biased_corrected_recons(:,round(x_max):520)=(1/Calib.alpha) *(max(mean(asin(sqrt(I_over_I0))))*2-asin(sqrt(I_over_I0(:,round(x_max):end,1))));
            
            Output.E_cross_section_average_corrected(:,round(x_max):end,1)=mean(E_biased_corrected_recons(:,round(x_max):end,1),'omitnan');
            
            figure
            plot([1:520]*9.3e-3,Output.E_cross_section_average_corrected(:,:,1));
            title(['CORRECTED-Pockels average profile sensor-' Output.sensor_name ', bias=' Output.bias, ', Xray tube current=' Output.flux])
            xlabel('thickness (mm)')
            ylabel('E field (V/m)')
            
            %Finds the area under the E-field between chosen cathode and anode
            y=[Output.E_cross_section_average_corrected(:,round(Calib.cathode),1) Output.E_cross_section_average_corrected(:,round(Calib.anode),1)];
            thickness=[Calib.cathode Calib.anode];
            Output.integral_Efield=trapz(thickness*9.3*10^-6,y);
            
            %%%%%%%%%%%%%%%%%%%%%%
        case 'No'
            %close(f1)
            %
            %Finds the area under the E-field between chosen cathode and anode
            y=[Output.E_cross_section_average_corrected(:,round(Calib.cathode),1) Output.E_cross_section_average_corrected(:,round(Calib.anode),1)];
            thickness=[Calib.cathode Calib.anode];
            Output.integral_Efield=trapz(thickness*9.3*10^-6,y);
            
    end
    if strcmp(calib_status,'No')
        %Finds the difference between bias applied and integral of Efield between chosen cathode and anode
        %         diff=str2num(Output.bias)-Output.integral_Efield;
        diff=Output.bias-Output.integral_Efield;
        
        %adjusts the r parameter and re-calculates E-field until the area
        %under E-field and applied bias are less than 30v different
        while abs(diff)>30
            if diff>0
                Calib.r_41=Calib.r_41-.1e-12;
                Calib.alpha=(sqrt(3)*pi*Calib.n_0^3*Calib.r_41*Calib.d)/(2*Calib.lambda_0);
                [Output]=Calculate_Efield(Output, Calib.alpha);
                y=[Output.E_cross_section_average_corrected(:,round(Calib.cathode),1) Output.E_cross_section_average_corrected(:,round(Calib.anode),1)];
                Output.integral_Efield=trapz(thickness*9.3*10^-6,y);
                Calib.r_41=Calib.r_41;
                %                 diff=str2num(Output.bias)-Output.integral_Efield;
                diff=Output.bias-Output.integral_Efield;
            else
                Calib.r_41=Calib.r_41+.1e-12;
                alpha=(sqrt(3)*pi*Calib.n_0^3*Calib.r_41*Calib.d)/(2*Calib.lambda_0);
                [Output]=Calculate_Efield(Output, Calib.alpha);
                y=[Output.E_cross_section_average_corrected(:,round(Calib.cathode),1) Output.E_cross_section_average_corrected(:,round(Calib.anode),1)];
                Output.integral_Efield=trapz(thickness*9.3*10^-6,y);
                Calib.r_41=Calib.r_41;
                %                 diff=str2num(Output.bias)-Output.integral_Efield;
                diff=Output.bias-Output.integral_Efield;
            end
            
        end
    end
    if Output.flux>0
        figure(2)
        plot([1:520],Output.E_cross_section_average_corrected(:,:,1),'displayname',[bias_string 'V_' flux_string 'mA']);
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
        plot([1:520],Output.E_cross_section_average_corrected(:,:,1),'--','displayname',[bias_string 'V_' flux_string 'mA']);
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
    Output.version='1.2';
    %--------------
    %     waitfor(msgbox('Please enter the directory where you would like to save the results'))
    %     Dir = uigetdir('X:\');
    Dir=path_biased_crossed;
    save([Dir '\' Output.sensor_name '_' bias_string 'V_' flux_string 'mA_Pockels_output.mat'], 'Output');
    figure(4)
    savefig([Output.sensor_name '_' bias_string 'V_' flux_string 'mA_E-field_uncorrected'])
    if strcmp(calib_status,'No')
        save([Dir '\' Output.sensor_name '_Pockels_calib_file.mat'], 'Calib');
    end
    calib_status='yes';
    %     newname=['E_field_' Output.sensor_name '_' bias_string 'V_' flux_string 'mA'];
    %     eval(sprintf('%s = %g',newname,Output.E_cross_section_average_corrected));
    %     newname=['integral_E_field_' Output.sensor_name '_' bias_string 'V_' flux_string 'mA'];
    %     eval(sprintf('%s = %g',newname,Output.integral_Efield));
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
clearvars -except Output*