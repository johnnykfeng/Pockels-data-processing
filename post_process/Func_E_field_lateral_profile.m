% calculate and plot the lateral E-field profile
function [Profile]=Func_E_field_lateral_profile(Output,edge)
E=Output.E_field_biased_corrected;
Profile=mean(E,2,'omitnan');
bias_string=num2str(Output.bias);
flux_string=num2str(Output.flux);
figure(1)
plot(Profile,'displayname',[Output.sensor_name '@' bias_string 'V,' flux_string 'mA']);
title(['lateral average profile sensor-' Output.sensor_name '@' bias_string 'V,' flux_string 'mA'])
xlabel('thickness')
ylabel('E field (V/m)')
%%
Prof_smoo=movmean(Profile,10);
x=1:length(Prof_smoo);
x=x-edge;
figure(2)
plot(x,Prof_smoo/max(Prof_smoo),'displayname',[Output.sensor_name '@' bias_string 'V,' flux_string 'mA']);
plot([0 0],[0 1.1],'k--','displayname','edge')
title(['lateral average profile sensor-' Output.sensor_name '@' bias_string 'V,' flux_string 'mA'])
xlabel('thickness')
ylabel('E field (V/m)')
box
legend
end
