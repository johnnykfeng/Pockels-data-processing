function [x,V,Energy]=Func_potential_profile_figureVersion(Calib)
% calculate the potential profile based on E-field profile
% open E-field profile figure first, then run this code
% 2021-08-13 created by Yuxin

%% read figure data
a=get(gca,'Children');
xdata_cell = get(a, 'XData');
ydata_cell = get(a, 'YData');
name_cell = get(a, 'DisplayName');
N=size(xdata_cell,1);% amount of curves in the figure
x=cell2mat(xdata_cell);% x-axis data
E=cell2mat(ydata_cell);% y-axis data

% clear xdata_cell ydata_cell
%%
E_clean=E;
E_clean(:,1:Calib.cathode)=0;
E_clean(:,Calib.anode:end)=0;
[W,H]=size(E);
V=zeros(W,H);
Energy=zeros(W,H);
figure(21)
hold on
title(['Potential profile [V]'])
xlabel('Position [mm]')
ylabel('Potential [V]')
figure(22)
hold on
title(['Band energy [eV]'])
xlabel('Position [mm]')
ylabel('Band energy [V]')
for i=1:W
%     V(i,1)=1;
    for j=2:Calib.anode
        V(i,j)=E_clean(i,j)+V(i,j-1);
    end
    V(i,:)=V(i,:)*Calib.scale*1e-3;
    Energy(i,:)=-V(i,:);
    name=cell2mat(name_cell(i));
    X=x(i,:)-(Calib.cathode-1)*Calib.scale;%*1e-3;
    figure(21)
    plot(X,V(i,:),'displayname', name)
    figure(22)
    plot(X,Energy(i,:),'displayname', name)
end
% xlim([0,])


box
end
