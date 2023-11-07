% to calculate space charge map based on E-field map
% 2021-11-26 created by Yuxin
function [rho]=Func_space_charge_map(Output, Calib)
order=5;%fitting order for polyfit
C.eps_0=8.8541878128E-12;% [F⋅m−1]permittivity of vacumm
C.q=1.602176634e-19; %[C]
C.epsilon=10.9;
eps_CZT=C.eps_0*C.epsilon;
[E_clean,Edge_cathode,Edge_anode]=Func_clean_E_field(Output,Calib);
E_smooth=movmedian(E_clean,20);
[W,H]=size(E_clean);

rho=zeros(W,H);
xx=Calib.x_all;% [mm] x-axis of the sensor region selected
for i=1:W
    E_line=(E_smooth(i,:));
%     E_line=movmean(E_line,10);
    idx = ~isnan(E_line);
    x_sensor=xx(idx);
    E_line_sensor=E_line(idx);
    fitResults = polyfit(x_sensor,E_line_sensor,order);% order-th order polynominal fit to the E-field profile
    AA=zeros(1,6);
    for j=1:order+1
        AA(end-j+1)=fitResults(end-j+1);
    end

    E_fit=AA(1)*x_sensor.^5+AA(2)*x_sensor.^4+AA(3)*x_sensor.^3+AA(4)*x_sensor.^2+AA(5)*x_sensor+AA(6);
    E_div=diff(E_fit)./(10^-5);% E_div is in SI unit [V/m^2]
    E_div2=[E_div 0];
    rho(i,idx)=-E_div2.*eps_CZT./10^6/C.q; %[e/cm^3] space charge density profile
end

bias_string=num2str(Output.bias);
flux_string=num2str(Output.flux);
figure
imagesc(rho)
title(['Space charge map ' Output.sensor_name '@' bias_string 'V,' flux_string 'mA'])
axes1=gca;
set(axes1,'DataAspectRatio',[1 1 1],'Layer','top');
set(axes1,'View',[-90 90]);
colorbar
box
end
