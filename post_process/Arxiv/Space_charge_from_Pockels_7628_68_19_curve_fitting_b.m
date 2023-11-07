% Pocels average profile sensor 7628-68-19, bias=700V
clear
load 'pockels_7628_68_19_good_lag.mat'
constants
eps_CZT=C.eps_0*C.epsilon;
%---settings----
dataspan1=128;
dataspan2=317;
smoothstep=10;
%---------------
% under X-ray
y1=ydata(2,:);
y1=y1(dataspan1:dataspan2);
y1=smoothdata(y1,'movmedian',smoothstep);
rho1=-diff(y1)./(10^-5)*eps_CZT./10^6/C.q; %e/cm^3;
rho1_smooth=smoothdata(rho1,'movmedian',smoothstep);
% dark
y2=ydata(1,:);
y2=y2(dataspan1:dataspan2);
y2=smoothdata(y2,'movmedian',smoothstep);
rho2=-diff(y2)./(10^-5)*eps_CZT./10^6/C.q; %e/cm^3;;
rho2_smooth=smoothdata(rho2,'movmedian',smoothstep);

rho_neg=rho2_smooth;
rho_pos=rho1_smooth-rho_neg;

% E-field component caused by positive space charge
y3=y1-y2;


figure(2)
hold on
% plot(rho1,'displayname','dark')
% plot(rho2,'displayname','X-ray')
plot(rho1_smooth,'displayname','X-ray')
plot(rho2_smooth,'displayname','dark')
plot(rho_pos,'displayname','positive')
plot([0 200],[0 0],'--k','displayname','0-reference')
hold off
xlabel('position [10um]')
ylabel('Space charge [e/cm^3]')
legend
grid
box

% figure(1)
% hold on
% x12=[dataspan1:dataspan2];
% plot(xdata(1,:),ydata(1,:),'displayname','dark')
% plot(x12,y1,'displayname','X-ray')
% plot(xdata(2,:),ydata(2,:),'displayname','X-ray')
% plot(x12,y2,'displayname','dark')
% hold off
% grid
% box
%%
figure(1)
hold on
x12=[dataspan1:dataspan2];
% plot(xdata(1,:),ydata(1,:),'displayname','dark')
plot(y1,'displayname','X-ray')
% plot(xdata(2,:),ydata(2,:),'displayname','X-ray')
plot(y2,'displayname','dark')
hold off
xlabel('position [10um]')
ylabel('E-field [V/m]')
legend
grid
box
%%
% figure(1)
% hold on
% x12=[dataspan1:dataspan2];
% % plot(xdata(1,:),ydata(1,:),'displayname','dark')
% plot(ydata(1,:),'displayname','X-ray 1')
% % plot(xdata(2,:),ydata(2,:),'displayname','X-ray')
% plot(ydata(2,:),'displayname','dark 2')
% plot(ydata(3,:),'displayname','dark 3')
% plot(ydata(4,:),'displayname','dark 4')
% hold off
% xlabel('position [10um]')
% ylabel('E-field [V/m]')
% legend
% grid
% box