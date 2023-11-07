% calculate space charge profile based on Pockels line profile
% 2021-10-20 Created by Yuxin
% "order": (1~5) is the order of polynominal fitting, for a linear fitting, order=1
% "cut": (default 0) data points to be removed adjacent to the two contacts for the fitting

% function [rho,xx,coef]=Func_space_charge_line_profile(Output,Calib,order,cut)
constants
eps_CZT=C.eps_0*C.epsilon;
E=Output.E_cross_section_average_corrected(Calib.cathode:Calib.anode);
xx=Calib.x_sensor;% [mm] x-axis of the sensor region selected
f1=figure(31);
plot(E)
xlabel('a.u.')
ylabel('E=field [V/m]')

%%
waitfor(msgbox(['Select the region of interest']));
[points,y]=ginput(2);
if points(1)<1
    points(1)=1;
end
points=round(points);
E_selected=E(points(1):points(2));%select a region the integation of E-field roughly equals to the bias applied
xx2=xx(points(1):points(2));
order=1;
fitResults = polyfit(xx2,E_selected,order);% order-th order polynominal fit to the E-field profile
AA=zeros(1,6);
for i=1:order+1
    AA(end-i+1)=fitResults(end-i+1);
end

E_fit=AA(1)*xx2.^5+AA(2)*xx2.^4+AA(3)*xx2.^3+AA(4)*xx2.^2+AA(5)*xx2+AA(6);
f2=figure(32);
plot(xx2,E_selected)
hold on
plot(xx2,E_fit,'k--')

E_div=diff(E_fit)./(10^-5);% E_div is in SI unit [V/m^2]
rho=-E_div(1).*eps_CZT./10^6/C.q %[e/cm^3] space charge density profile


%% plot figure and check fitting qualtiy
% check=1;
% while check==1
%     
%     R = corrcoef(E_fit,E);
%     coef=R(2);% Correlation coefficient as an indicator for the fitting quality
%     f91=figure(91);
%     movegui(f91,[40 575]);%movegui(f1,[40 903]);
%     f92=figure(92);
%     movegui(f92,[612 575]);%movegui(f2,[612 903]);
%     figure(91)
%     plot(xx,E)
%     hold on
%     plot(xx,E_fit)
%     
%     figure(92)
%     plot(xx(1:end-1),rho)
%     
%     
%     if coef<0.95
%         % waitfor(msgbox(['check the fitting quality']))
%         coef_str=num2str(coef);
%         redo = questdlg(['The fitting coefficient ('  coef_str ') is low. Do you want to select a region and redo the fitting?']);
%         switch redo
%             case 'Yes'
%                 % select a region and redo the curve fitting
%                 figure(91)
%                 waitfor(msgbox(['Select the region of interest']));
%                 [m,n]=ginput(2);
%                 points=round((m-xx(1))/(xx(end)-xx(1))*length(xx));
%                 
%                 
%                 E=E(points(1):points(2));%
%                 xx=xx(points(1):points(2));
%                 fitResults = polyfit(xx,E,order);% order-th order polynominal fit to the E-field profile
%                 AA=zeros(1,6);
%                 for i=1:order+1
%                     AA(end-i+1)=fitResults(end-i+1);
%                 end
%                 
%                 E_fit=AA(1)*xx.^5+AA(2)*xx.^4+AA(3)*xx.^3+AA(4)*xx.^2+AA(5)*xx+AA(6);
%                 
%                 
%             case 'No'
%                 check=0;
%         end
%     else
%         check=0;
%     end
%     E_div=diff(E_fit)./(10^-5);% E_div is in SI unit [V/m^2]
%     rho=-E_div.*eps_CZT./10^6/C.q; %[e/cm^3] space charge density profile
% end
% 
% 
% pause(2)
% close 91 92

% end