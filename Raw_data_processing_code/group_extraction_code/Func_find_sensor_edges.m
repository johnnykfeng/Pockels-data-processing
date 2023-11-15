function [Edge_left, Edge_right] = Func_find_sensor_edges(I0)
    %FUNC_FIND_SENSOR_EDGES Detects the left and right edges of signals in an image.
    %
    % Inputs:
    %   I0: A 2D matrix representing the input image. The matrix elements are
    %       expected to be intensity values of the image.
    %
    % Outputs:
    %   Edge_left: 1D array indicating left edge positions of signals per row.
    %   Edge_right: 1D array indicating right edge positions of signals per row.
    %
    % Algorithm:
    %   1. Determine the size of the input image and extract the center row.
    %   2. Smooth the center row using a moving mean filter with a window size of 20.
    %   3. Set the signal detection threshold to one-fourth of the maximum value
    %      in the smoothed center row.
    %   4. For each row in the image:
    %       a. Find intensity values above the threshold.
    %       b. Determine the median of these values and set a new threshold to half of this median.
    %       c. Find the range in the row where intensity values exceed this new threshold.
    %       d. Record the start and end indices of this range as the left and right edges.
    %   5. If no signal is detected in a row, default the edges to the start and end of the row.
    %
    % Note:
    %   This function is designed for images where signals are represented by higher
    %   intensity values against a darker background. It may not work as intended for
    %   images with different characteristics.
    disp('Calling Func_find_sensor_edges...')

    dimension = size(I0);
    image_width = dimension(1); % 696
    image_length = dimension(2); % 520

    %----- threshold for judging signal
    I0_center = I0(round(image_width / 2), :);
    I0_center_smooth = movmean(I0_center, 20);
    threshold = 1/4 * max(I0_center_smooth);

    % initialize array
    Edge_left = zeros(1, image_width);
    Edge_right = zeros(1, image_width);

    for i = 1:image_width
        Iline = I0(i, :); % I0 per line
        I_signal_rough = Iline(Iline > threshold);
        Signal_median = median(I_signal_rough);
        Th = Signal_median / 2; % threshold for judging a edge
        signal_range = find(Iline > Th);

        if ~isempty(signal_range)
            Edge_left(i) = signal_range(1);
            Edge_right(i) = signal_range(end);
        else
            Edge_left(i) = 1;
            Edge_right(i) = image_length;
        end

    end

end