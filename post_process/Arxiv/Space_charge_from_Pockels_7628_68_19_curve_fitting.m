% Pocels average profile sensor 7628-68-19, bias=700V
clear
constants
eps_CZT=C.eps_0*C.epsilon;

x=0:200;
% under X-ray
y1=4.0289*x.^2-525.1*x+373295;
% dark
% y2=2.22967*x.^2+456.13*x+304744;
y2=900.45*x+290177;


% under X-ray
rho1=-diff(y1)./(10^-5)*eps_CZT./10^6/C.q; %e/cm^3;
% dark
rho2=-diff(y2)./(10^-5)*eps_CZT./10^6/C.q; %e/cm^3;;


rho_neg=rho2;
rho_pos=rho1-rho_neg;

figure(2)
hold on
% plot(rho1,'displayname','dark')
% plot(rho2,'displayname','X-ray')
plot(rho1,'displayname','X-ray')
plot(rho2,'displayname','dark')
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
% x12=[dataspan1:dataspan2];
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