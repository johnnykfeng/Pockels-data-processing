function [Output, Calib, calib_status] = calibration_ini(parameter, sensor_name, sensor_thickness)
    %% DOCSTRING
    % Routine for calibration step of Pockels data analysis
    % 
    % Inputs:
    %   - parameter: 2 x N array of bias and x-ray currents used for measurements
    %   - sensor_name: user input sensor name
    %   - sensor_thickness: user input sensor thickness in mm
    % 
    % Outputs:
    %   - Output: struct containing the extracted data
    %   - Calib: struct containing variables for calibration
    %   - calib_status: boolean input for using pre-existing calibration
    % 
    calib_status = questdlg('Do you want to use existing calibration file?');

    if strcmp(calib_status, 'Yes')
        waitfor(msgbox(['Select your calibration data']));
        [Calib_file, Calib_path] = uigetfile('*.mat');
        load([Calib_path '\' Calib_file]);
        Output = nan;

    else

        if strcmp(calib_status, 'No')
            Calib.n_0 = 2.8;
            Calib.r_41 = 5.6 * 10^-12;
            Calib.d = 8.748 * 10^-3;
            Calib.lambda_0 = 980 * 10^-9; % wavelength
            Calib.alpha = (sqrt(3) * pi * Calib.n_0^3 * Calib.r_41 * Calib.d) / (2 * Calib.lambda_0);
            % Calib is an output of the function, but I don't think it's useful

            %% Input and plot 0V Parallel
            waitfor(msgbox(['Select the 0V, parallel data']));
            [file_0V_parallel, path_0V_parallel] = uigetfile('*.txt'); % the name of the file doesn't matter
            disp([file_0V_parallel, path_0V_parallel])
            Output.data_directory = path_0V_parallel;

            buffer1 = fileread(char([path_0V_parallel file_0V_parallel]));
            data1 = textscan(buffer1, '%f %f %f', 'headerlines', 19); % skip the first 19 header lines
            I_parallel_nobias_matrix(:, 1) = data1{1, 1}; % Rows
            I_parallel_nobias_matrix(:, 2) = data1{1, 2}; % Columns
            I_parallel_nobias_matrix(:, 3) = data1{1, 3}; % Pixel Values
            %[I_parallel_nobias_matrix(:,1) I_parallel_nobias_matrix(:,2) I_parallel_nobias_matrix(:,3)] = textread([path_0V_parallel  file_0V_parallel], '%f %f %f', 'headerlines', 19) ;

            I_parallel_nobias_field = reshape(I_parallel_nobias_matrix(:, 3), [696, 520]); % converts 1-D array to 2-D array of shape

            figure_parallel = figure('Name', '0V Parallel');
            movegui([1184 575]);
            imagesc(I_parallel_nobias_field)
            title(['Intensity @ parallel polarizers'])
            axes1 = gca;
            set(axes1, 'DataAspectRatio', [1 1 1], 'Layer', 'top');
            set(axes1, 'View', [-90 90]);
            box

            %% Input and plot 0V Crossed
            waitfor(msgbox(['Select the 0V, crossed data']));
            [file_0V_crossed, path_0V_crossed] = uigetfile('*.txt');

            buffer2 = fileread(char([path_0V_crossed file_0V_crossed]));
            data2 = textscan(buffer2, '%f %f %f', 'headerlines', 19);
            I_crossed_nobias_matrix(:, 1) = data2{1, 1};
            I_crossed_nobias_matrix(:, 2) = data2{1, 2};
            I_crossed_nobias_matrix(:, 3) = data2{1, 3};
            %[I_crossed_nobias_matrix(:,1) I_crossed_nobias_matrix(:,2) I_crossed_nobias_matrix(:,3)] = textread([path_0V_crossed  file_0V_crossed], '%f %f %f', 'headerlines', 19) ;

            I_crossed_nobias_field = reshape(I_crossed_nobias_matrix(:, 3), [696, 520]);

            figure_crossed = figure('Name', '0V Crossed');
            movegui([1184 66]);
            imagesc(I_crossed_nobias_field)
            title(['Intensity @ crossed polarizers'])
            axes1 = gca;
            set(axes1, 'DataAspectRatio', [1 1 1], 'Layer', 'top');
            set(axes1, 'View', [-90 90]);
            %         colorbar
            box

            Calib.I_parallel_nobias_field = I_parallel_nobias_field;
            Calib.I_crossed_nobias_field = I_crossed_nobias_field;

            %% select a region where the E-field is unform laterally and intensity is high for calibration
            waitfor(msgbox(['Please select one raw data for calibration purpose (Recommend 500~700V, no X-ray data)']))
            %-------------

            answer = ...
                inputdlg({'bias for calibration'}, ...
                'Calibration data selection', ...
                [1, 35], ...
                {'700'});

            bias_string = answer{1, 1};
            flux_string = '0'; %answer{2,1};
            calib_bias = str2double(bias_string);
            Calib.calibrationVoltage = calib_bias;
            calib_flux = 0; %str2num(flux_string);
            %-------------
            % load the raw data on the first row of P.txt for calibration.
            % usually bias between 500 to 700V is good to avoid reaching the
            % peak of sin() function
            [Output] = Read_Pockels_data_group_extraction(sensor_name, calib_bias, calib_flux);
            [Output] = Calculate_Efield(Output, Calib);

            % plot the E-field of the chosen calibration data
            figure_calib = figure('Name', 'Output.E_cross_section_average_corrected');
            plot(Output.E_cross_section_average_corrected, 'displayname', [bias_string 'V,' flux_string 'mA']);
            title(['Pockels average profile sensor-' Output.sensor_name '@' bias_string 'V,' flux_string 'mA'])
            xlabel('Cathode to anode distance (pixel)')
            ylabel('E field (V/m)')

            % plot the image of corrected E-field
            figure_calib_corrected = figure('Name', 'Output.E_field_biased_corrected');
            imagesc(Output.E_field_biased_corrected)
            % title(['E-field sensor-' Output.sensor_name '@' bias_string 'V,' flux_string 'mA'])
            axes1 = gca;
            set(axes1, 'DataAspectRatio', [1 1 1], 'Layer', 'top');
            set(axes1, 'View', [-90 90]);
            colorbar
            box

            E_field_biased_corrected = Output.E_field_biased_corrected;
            I0 = Calib.I_parallel_nobias_field;
            %% roughly select a region in 0V parallel image to remove obvious regions with light other than the sensor
            % added on 2022/03/11 by Yuxin
            % for version 2.4
            I0_background = min(min(I0));

            %% Select ROI for Calibration 1
            fig_select_ROI = figure('Name', 'Select ROI for Calibration');
            imagesc(I0)
            waitfor(msgbox(['Select the region of interest']));
            [points, y] = ginput(2);
            points = round(points);
            region_left = points(1);
            Calib.rough_region_left = region_left;
            region_right = points(2);
            Calib.rough_region_right = region_right;
            I0(:, region_right:end) = I0_background;
            imagesc(I0)
            I0(:, 1:region_left) = I0_background;
            imagesc(I0)

            %%
            dimension = size(E_field_biased_corrected);
            N0 = dimension(1);
            L = dimension(2);
            %-----find threshold for judging signal
            I0_center = I0(round(N0 / 2), :);
            I0_center_smooth = movmean(I0_center, 20);
            thre = 1/4 * max(I0_center_smooth);

            %--------------------------------------
            % find edge
            for i = 1:N0
                Iline = I0(i, :);
                %             I_signal_rough=Iline(Iline>200);
                I_signal_rough = Iline(Iline > thre);
                Signal_median = median(I_signal_rough);
                Th = Signal_median / 2; %threshold for judging a edge
                signal_range = find(Iline > Th);

                if ~isempty(signal_range)
                    Edge_left(i) = signal_range(1);
                    Edge_right(i) = signal_range(end);
                else
                    Edge_left(i) = 1;
                    Edge_right(i) = L;
                end

            end

            % find outliers
            Edge_left_median = round(median(Edge_left));
            Edge_right_median = round(median(Edge_right));
            T = Edge_right - Edge_left; % thickness of sensor in unite of pixel
            T_median = median(T);
            Calib.scale = sensor_thickness / T_median; %[mm/pixel] the corresponding thickness of one pixel
            outlier = find((abs(T - T_median) / T_median > 0.2) == 1);

            for i = 1:N0
                Eline = E_field_biased_corrected(i, Edge_left(i):Edge_right(i));
                x = (0:T(i)) / T(i) * sensor_thickness * 1E-3;
                E_int(i) = trapz(x, Eline);
            end

            %% Select ROI for Calibration 2
            % f6 = figure(6);
            fig_select_ROI_b = figure('Name', 'Select ROI for Calibration');
            plot(E_int)
            waitfor(msgbox(['Select the region of interest']));
            [points, y] = ginput(2);
            points = round(points);
            Calib.selected_region_edge1 = points(1);
            Calib.selected_region_edge2 = points(2);
            E_selected = Output.E_field_biased_corrected(points(1):points(2), :); %select a region the integation of E-field roughly equals to the bias applied
            E_profile_selected = mean(E_selected, 'omitnan'); %

            f7 = figure(7);
            movegui(f7, [1184 66]);
            % movegui(f3,[40 66]);%movegui(f3,[40 392]);
            imagesc(E_selected)
            title('selected E-field sensor-')
            axes1 = gca;
            set(axes1, 'DataAspectRatio', [1 1 1], 'Layer', 'top');
            set(axes1, 'View', [-90 90]);
            colorbar
            box

            Edge_left_sel = Edge_left(points(1):points(2));
            Edge_right_sel = Edge_right(points(1):points(2));
            Edge_left_sel_median = median(Edge_left_sel);
            Edge_right_sel_median = median(Edge_right_sel);
            T_sel = Edge_right_sel - Edge_left_sel; % thickness of sensor in unite of pixel
            T_sel_median = median(T_sel);

            xs = (0:(Edge_right_sel_median - Edge_left_sel_median)) / T_sel_median * sensor_thickness * 1E-3; %[m]
            E_integral_selected = trapz(xs, E_profile_selected(Edge_left_sel_median:Edge_right_sel_median)); %integral of E-field in the seleteced region
            Output.integral_Efield = E_integral_selected;

            x_all = [1:L] * (xs(2) - xs(1)) * 1e3; %[mm]
            Calib.x_sensor = (0:(Edge_right_median - Edge_left_median)) / T_sel_median * sensor_thickness; %[mm]
            Calib.x_all = x_all; %[mm]
            %         Calib.xs=xs;
            %% adjusts the alpha parameter to match the integral E-field with bias applied
            ratio = E_integral_selected / Output.bias;
            Calib.alpha = Calib.alpha * ratio;
            %adjusts the r parameter and re-calculates E-field until the area
            Calib.r_41 = Calib.alpha * (2 * Calib.lambda_0) / (sqrt(3) * pi * Calib.n_0^3 * Calib.d);
            [Output] = Calculate_Efield(Output, Calib);
            E_selected = Output.E_field_biased_corrected(points(1):points(2), :); %select a region the integation of E-field roughly equals to the bias applied
            E_profile_selected = mean(E_selected, 'omitnan'); %
            E_integral_selected = trapz(xs, E_profile_selected(Edge_left_sel_median:Edge_right_sel_median)); %integral of E-field in the seleteced region
            Output.integral_Efield = E_integral_selected;
            Calib.cathode = Edge_left_median;
            Calib.anode = Edge_right_median;

            f5 = figure('Name', 'Integrated E-field of selected calibration');
            movegui(f5, [1184 575]); %movegui(f2,[612 903]);
            %         E_selected=Output.E_field_biased_corrected(Calib.selected_region_edge1:Calib.selected_region_edge2,:);
            %         E_profile_selected=mean(E_selected,'omitnan');%
            plot(Calib.x_all, Output.E_cross_section_average_corrected, 'displayname', 'all');
            hold on
            plot(Calib.x_all, E_profile_selected, 'displayname', 'selected')
            legend
            y = Output.E_cross_section_average_corrected(round(Calib.cathode):round(Calib.anode));
            xs2 = (0:(Edge_right_median - Edge_left_median)) / T_median * sensor_thickness;
            Output.integral_Efield = trapz(xs2, y);
            int_E_all = num2str(Output.integral_Efield);
            int_E_selected = num2str(E_integral_selected);
            title(['Integration of E-field: all=' int_E_all 'V, Selected=' int_E_selected 'V'])
            xlabel('Distance (mm)')
            ylabel('E field (V/m)')
            box
            waitfor(msgbox(['Calibration done']));
            close 5 6 7 8 9 10 11
            % close figure_parallel figure_crossed fig_select_ROI fig_select_ROI_b
        end

    end

end
