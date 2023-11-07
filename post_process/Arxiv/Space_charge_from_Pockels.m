% Pocels average profile sensor 7534-88-10, bias=700V
% clear
load 'data.mat'
constants
eps_CZT=C.eps_0*C.epsilon;
% %---settings----
% dataspan1=128;
% dataspan2=317;
% smoothstep=10;
%---------------
x=1:183;
%dark
y_dark=3.6745.*x.^2+226.37*x+314031;% 2nd order polynomial assumption
% y_dark=898.8*x+293410;% linear assumption
rho_dark=-diff(y_dark)./(10^-5)*eps_CZT./10^6/C.q; %[e/cm^3];

% under X-ray
y_25mA=3.6397.*x.^2-939.94*x+381744;
rho_x_ray=-diff(y_25mA)./(10^-5)*eps_CZT./10^6/C.q; %[e/cm^3];
%
rho_net=rho_x_ray-rho_dark;%net charge change when X-ray is on


figure(2)
hold on
plot(rho_x_ray,'displayname','X-ray')
plot(rho_dark,'displayname','dark')
plot(rho_net,'displayname','net increase')
plot([0 200],[0 0],'--k','displayname','0-reference')
hold off
xlabel('position [10um]')
ylabel('Space charge [e/cm^3]')
legend
title('7628-68-19@700V (negative lag)')
grid
box
%%
figure(1)
hold on
plot(profile_25mA,'displayname','25mA')
plot(y_25mA,'k--','displayname','')
plot(profile_0mA,'displayname','dark')
plot(y_dark,'k--','displayname','')
hold off
xlabel('position [10um]')
ylabel('E-field [V/m]')
legend
grid
box
