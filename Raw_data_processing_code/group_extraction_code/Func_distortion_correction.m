function [image_corrected, Edge_cathode_fit, Edge_anode_fit] = Func_distortion_correction(I0, region_left, region_right, E_biased_raw)
    % disp('calling Func_distortion_correction ...')

    I0(:, region_right:end) = NaN;
    I0(:, 1:region_left) = NaN;

    [Edge_cathode, Edge_anode] = Func_find_sensor_edges(I0);

    % smooth edges
    [M, N] = size(I0);
    x = 1:M;
    x_clean = x;
    [cleanedData, index] = cleanData(Edge_cathode);
    x_clean(index) = [];
    p = polyfit(x_clean, cleanedData, 3);
    % Edge_cathode = polyval(p, x_clean);
    Edge_cathode_fit = p(1) * x.^3 + p(2) * x.^2 + p(3) * x + p(4);
    Edge_cathode_fit(Edge_cathode_fit < 1) = 1;

    [cleanedData, index] = cleanData(Edge_anode);
    x = 1:M;
    x_clean = x;
    x_clean(index) = [];
    p = polyfit(x_clean, cleanedData, 3);
    % Edge_anode = polyval(p, x_clean);
    Edge_anode_fit = p(1) * x.^3 + p(2) * x.^2 + p(3) * x + p(4);
    Edge_anode_fit(Edge_anode_fit > N) = N;

    % distortion correction for each vertical line
    image = E_biased_raw;
    % figure (85)
    % imagesc(image)
    % axes1=gca;
    % set(axes1,'DataAspectRatio',[1 1 1],'Layer','top');
    % set(axes1,'View',[-90 90]);

    % e1 = mean(Edge_cathode_fit);
    % e2 = mean(Edge_anode_fit);
    % d=round(e2-e1);%average thickness
    sensor_thickness_fit = mean(Edge_anode_fit - Edge_cathode_fit); %average thickness
    % X_standard=round(e1):round(e2);

    % center=mean(e1+e2);
    image_corrected = nan(M, N); % create empty matrix of the same dimensions
    delta_x = 50; %expand the region to show in image

    if sensor_thickness_fit + 2 * delta_x >= N
        %     delta_x=floor((N-d)/2);
        delta_x = 0;
    elseif min(Edge_cathode_fit) < delta_x
        delta_x = 0;
    elseif max(Edge_anode_fit) > N - delta_x
        delta_x = 0;
    end

    X_standard = (round(mean(Edge_cathode_fit)) - delta_x):(round(mean(Edge_anode_fit)) + delta_x);

    % Creates the image_corrected line by line
    for i = 1:M
        di = Edge_anode_fit(i) - Edge_cathode_fit(i);
        x0 = 0:di + 2 * delta_x;
        bottom_image_edge = (Edge_cathode_fit(i) - delta_x); %#ok<NOPRT>
        top_image_edge = (Edge_anode_fit(i) + delta_x); %#ok<NOPRT>
        % y0=image(i, (Edge_cathode_fit(i)-delta_x):(Edge_anode_fit(i) + delta_x));
        y0 = image(i, bottom_image_edge:top_image_edge);

        x1 = (1:length(X_standard)) / length(X_standard) * (sensor_thickness_fit + 2 * delta_x);
        y1 = interp1(x0, y0, x1);
        image_corrected(i, X_standard(1):X_standard(end)) = y1;
    end

    % figure (86)
    % imagesc(image_corrected)
    % axes1=gca;
    % set(axes1,'DataAspectRatio',[1 1 1],'Layer','top');
    % set(axes1,'View',[-90 90]);
end

function [cleanedData, index] = cleanData(data)
    % flatten the curve
    x = 1:length(data);
    p = polyfit(x, data, 3);
    yFit = polyval(p, x);
    yFlat = data - yFit;

    % Compute the median absolute deviation (MAD) of the data
    medianData = median(yFlat);
    madData = median(abs(yFlat - medianData));

    % Compute the z-scores of the data
    zScores = (yFlat - medianData) / madData;

    % Find the data points that have z-scores greater than 3
    index = find(abs(zScores) > 5);

    % Remove the outlying data points
    cleanedData = data;
    cleanedData(index) = [];
end