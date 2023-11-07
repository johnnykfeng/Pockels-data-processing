function [V]=Func_potential_map(Output, Calib)
% calculate the potential map based on E-field map
% 2021-08-06 created by Yuxin
[E_clean,Edge_cathode,Edge_anode]=Func_clean_E_field(Output,Calib);
[W,H]=size(E_clean);
E_clean2=E_clean;
E_clean2(isnan(E_clean2))=0;
V=zeros(W,H);
for i=1:W
    V(i,1)=1;
    for j=2:Edge_anode(i)
        V(i,j)=E_clean2(i,j)+V(i,j-1);
    end
end
V=V*Calib.scale*1e-3;
bias_string=num2str(Output.bias);
flux_string=num2str(Output.flux);
figure
imagesc(V)
title(['Potential map ' Output.sensor_name '@' bias_string 'V,' flux_string 'mA'])
axes1=gca;
set(axes1,'DataAspectRatio',[1 1 1],'Layer','top');
set(axes1,'View',[-90 90]);
colorbar
box
end
