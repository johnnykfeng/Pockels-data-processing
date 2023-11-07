% extract mu.tau-h based on net space charge profile
% based on method 2: calculate the mean net space charge under X-ray at 700V, fit
% to theoretical calculation of trapped hole defined by mutau
% 2021-11-09 Created by Yuxin
function [mu_tau_h]=Func_mu_tau_h_method2(rho_net_mean,current)
%%
load 'muTau_vs_space_charge.mat'
% 'muTau_vs_space_charge.mat' was obtained from code:Main_mutau_sweep.m
% y=current/25*y;% the trapped hole density needs to be adjusted based on the X-ray flux(tube current)
mu_tau_h=x(abs(rho_net_mean-y)==min(abs(rho_net_mean-y)));
if isempty(mu_tau_h)
    mu_tau_h=NaN;
else
    if mu_tau_h<=min(x) | mu_tau_h>=max(x)
        fprintf(2,'Warning:the mu.tau-h is beyond the lower limit, the value can be inaccurate');
        mu_tau_h=NaN;
    end
end
end