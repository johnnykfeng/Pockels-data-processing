% Script to launch GUI for P_file_generator

% Combined prompts for all inputs
prompt = {'Enter start bias value:', ...
    'Enter end bias value:', ...
    'Enter flux values (comma-separated):',...
    'Enter filename'};
dlgtitle = 'Input Parameters';
dims = [1 50]; % Adjust dimensions as needed
definput = {'100', '900', '0, 5, 25', 'P.txt'}; % Default inputs for start_bias, end_bias, and flux_values
answer = inputdlg(prompt, dlgtitle, dims, definput);

% Check if the user provided the input or canceled the operation
if isempty(answer)
    disp('User cancelled the operation.');
    return;
end

% Extract and convert inputs
start_bias = str2double(answer{1});
end_bias = str2double(answer{2});
flux_values_str = strrep(answer{3}, ' ', ''); % Remove any spaces from flux_values
flux_values = str2num(flux_values_str); % Convert to numeric array
filename = answer{4};

% Validate inputs
if isnan(start_bias) || isnan(end_bias) || isempty(flux_values)
    disp('Invalid inputs. Please enter numeric values.');
    return;
end

% Call the P_file_generator function
P_file_generate(start_bias, end_bias, flux_values, filename);

% Display a pop-up success message
msg_str = sprintf('Success! The %s has been generated.', filename);
msgbox(msg_str, 'Operation Completed');


function parameters = P_file_generate( ...
    start_bias, end_bias, flux_values, filename)

    % Generate bias values
    bias_values = start_bias:100:end_bias;
    
    % Initialize the output matrix
    parameters = zeros(length(bias_values) * length(flux_values), 2);
    
    % fileID = fopen('P.txt','w');
    
    % Loop to fill the matrix
    index = 1;
    for i = 1:length(bias_values)
        for j = 1:length(flux_values)
            row = [bias_values(i), flux_values(j)];
            % fprintf(fileID, num2str(row)); % print row to txt file
            % fprintf('\n');
            parameters(index, :) = row; % store row in matrix
            index = index + 1;
        end
    end
    % fclose(fileID);
    writematrix(parameters, filename, 'Delimiter', 'tab')
end
