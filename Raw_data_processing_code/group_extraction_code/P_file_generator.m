function parameters = P_file_generator(start_bias, end_bias, flux_values)

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
    writematrix(parameters, 'P.txt', 'Delimiter', 'tab')
end




