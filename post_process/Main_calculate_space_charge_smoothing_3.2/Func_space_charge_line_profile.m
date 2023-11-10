% calculate space charge profile based on Pockels line profile
% 2021-10-20 Created by Yuxin
% "order": (1~5) is the order of polynominal fitting, for a linear fitting, order=1
% "cut": (default 0) data points to be removed adjacent to the two contacts for the fitting

function [rho, xx] = Func_space_charge_line_profile(Output, Calib)
%% DOCSTRING
% Func_space_charge_line_profile - Calculate and visualize space charge density profile
%   The generated figures display the E-field and space charge profiles for a specific
%   bias voltage and flux current, with legends and axis labels for reference.
%
%   Inputs:
%   - Output: Structure containing E-field data and bias information.
%   - Calib: Structure containing sensor calibration parameters.
%
%   Outputs:
%   - rho: Space charge density profile in units of [e/cm^3].
%   - xx: Anode-cathode distance in millimeters, representing the x-axis of the sensor region.
%
%   The function performs the following steps:
%   1. Imports relevant constants from the 'constants.m' file.
%   2. Calculates the space charge density profile (rho) using E-field data and sensor
%      calibration parameters.
%   3. Generates rolling mean of the E-field data and calculates the correlation coefficient.
%   4. Creates two plots to visualize the E-field and space charge density profiles.
%%
    constants  % import variables from constants.m
    eps_CZT = C.eps_0 * C.epsilon;
    bias = num2str(Output.bias);
    flux = num2str(Output.flux);

    %% moving average
    E = Output.E_cross_section_average_corrected(Calib.cathode:Calib.anode);
    xx = Calib.x_sensor; % [mm] x-axis of the sensor region selected

    E_fit = movmean(E, 10); % make rolling mean of E
    % E_div2=diff(E_fit2)./(10^-5);% E_div is in SI unit [V/m^2]
    E_div = diff(E_fit) ./ ((xx(2) - xx(1)) * 1e-3); % E_div is in SI unit [V/m^2]

    rho_raw=-E_div .* eps_CZT ./ 10^6 / C.q; %[e/cm^3] space charge density profile
    rho = movmean(rho_raw, 10);

    R = corrcoef(E_fit, E);
    coef = R(2); % Correlation coefficient as an indicator for the fitting quality
    f91 = figure(91);
    movegui(f91, [40 575]); %movegui(f1,[40 903]);
    f92 = figure(92);
    movegui(f92, [612 575]); %movegui(f2,[612 903]);

    figure(91)
    plot(xx, E, '--');
    hold on;
    plot(xx, E_fit);
    title(['E-field - ' bias 'V, ' flux 'mA']);
    xlabel("anode-cathode distance [mm]");
    legend('raw E-field', 'rolling mean');
    grid;

    figure(92)
    plot(xx(1:end - 1), rho);
    hold on;
    title(['Space charge - ' bias 'V, ' flux 'mA']);
    xlabel("anode-cathode distance [mm]");
    grid;

    pause(1)
    close 91 92

end
