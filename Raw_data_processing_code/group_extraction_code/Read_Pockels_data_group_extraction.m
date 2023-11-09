function [Output] = Read_Pockels_data_group_extraction(sensor_name, bias, X_ray_tube_current)
    %% DOCSTRING
    % Read_Pockels_data_group_extraction reads Pockels data from *.txt file based on filename
    % then extracts the pixel values as an array size [696, 520, 1]
    % 
    % Inputs:
    %   - sensor_name: name of the sensor
    %   - bias: bias voltage of measurement
    %   - X_ray_tube_current: X-ray tube current value
    % 
    % Outputs:
    %   - Output: struct containing the extracted data
    %   Output contain these fields: [sensor_name, bias, flux, I_crossed_biased_field]
    %   Output.I_crossed_biased_field is the main data array of interest
    % 
    Output.sensor_name = sensor_name;
    Output.bias = bias;
    Output.flux = X_ray_tube_current;
    bias_string = num2str(bias);
    flux_string = num2str(X_ray_tube_current);

    %% Matching and handling filenames
    %---------2022-09-08-new for Verson 3.0---------------------------------
    filename_pattern = ['*' bias_string '*' flux_string '*.txt']; 
    file_match = dir(fullfile(filename_pattern)); % list of files that match filename_pattern
    duplicate = size(file_match); % duplicate should be 1 most of the time

    if duplicate(1) > 1
    % loop through each duplicate filename and use the one with minimal length
        %     fprintf(2,'Warning:there are more than 1 data set with the same bias and current/');
        for i = 1:duplicate(1)
            names = file_match(i).name;
            namelength(i) = strlength(names);
        end
        filename = file_match(namelength == min(namelength)).name;

    else
        filename = file_match.name;
    end

    file_content = fileread(filename);
    sensor_name = Output.sensor_name;

    %% Parsing .txt file into matrix
    %---------2022-09-08-------------------------------------
    file_data = textscan(file_content, '%f %f %f', 'headerlines', 19);
    I_crossed_biased_matrix(:, 1) = file_data{1, 1}; %Row
    I_crossed_biased_matrix(:, 2) = file_data{1, 2}; %Col
    I_crossed_biased_matrix(:, 3) = file_data{1, 3}; %Value

    if isempty(I_crossed_biased_matrix)
    % if the file_data is unreadable, then return NaN matrix
        waitfor(msgbox(['The data for ' bias_string 'V' flux_string 'A is broken!!! press OK to skip']))
        Output.I_crossed_biased_field = NaN(696, 520);
    else
        % Reshapes the 3rd column into a 2D matrix
        I_crossed_biased_field(:, :, 1) = reshape(I_crossed_biased_matrix(:, 3, 1), [696, 520]);
        Output.I_crossed_biased_field(:, :, 1) = (single(I_crossed_biased_field(:, :, 1)));

    end

end
