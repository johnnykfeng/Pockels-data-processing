function [Output, Calib]=Pockels_Profile

%%%%%% Use existing calibration file or generate a calibration
%%%%%% file?
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
%%%%%Gets the intensity profiles for 0v parallel, 0v crossed and biased crossed
[Output]=Read_Pockels_data;

%%%%%Calculates the E-field
[Output]=Calculate_Efield(Output, Calib.alpha);

f1=figure ;
plot([1:520],Output.E_cross_section_average_corrected(:,:,1));
title(['Pockels average profile - sensor=' Output.sensor_name ', bias=' Output.bias, ', Xray tube current=' Output.flux])
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
        close(f1)
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
        title(['CORRECTED - Pockels average profile - sensor=' Output.sensor_name ', bias=' Output.bias, ', Xray tube current=' Output.flux])
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
    diff=str2num(Output.bias)-Output.integral_Efield;
    
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
            diff=str2num(Output.bias)-Output.integral_Efield;
        else
            Calib.r_41=Calib.r_41+.1e-12;
            alpha=(sqrt(3)*pi*Calib.n_0^3*Calib.r_41*Calib.d)/(2*Calib.lambda_0);
            [Output]=Calculate_Efield(Output, Calib.alpha);
            y=[Output.E_cross_section_average_corrected(:,round(Calib.cathode),1) Output.E_cross_section_average_corrected(:,round(Calib.anode),1)];
            Output.integral_Efield=trapz(thickness*9.3*10^-6,y);
            Calib.r_41=Calib.r_41;
            diff=str2num(Output.bias)-Output.integral_Efield;
        end
        
    end
end

%Saves the output files into a directory specified by user
waitfor(msgbox('Please enter the directory where you would like to save the results'))
Dir = uigetdir('X:\');
save([Dir '\' Output.sensor_name '_Pockels_output.mat'], 'Output');
if strcmp(calib_status,'No')
    save([Dir '\' Output.sensor_name '_Pockels_calib_file.mat'], 'Calib');
end
end