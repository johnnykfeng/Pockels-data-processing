% select a region in the E-field map that the integration of the E-field is
% roughly equal to the bias applied. then, recalculate the E-field profile
% and integration

% close all
f1=figure(1);
movegui(f1,[40 575]);%movegui(f1,[40 903]);
imagesc(Output.E_field_biased_corrected)
% title(['E-field sensor-' Output.sensor_name '@' bias_string 'V,' flux_string 'mA'])
axes1=gca;
set(axes1,'DataAspectRatio',[1 1 1],'Layer','top');
set(axes1,'View',[-90 90]);
colorbar
box

E=Output.E_field_biased_corrected;

version=str2double(Output.version);
if version<2
    I0=Output.I_parallel_nobias_field;
else
    I0=Calib.I_parallel_nobias_field;
end

dimention=size(E);
N=dimention(1);
L=dimention(2);
% find edge
% for i=1:N
%     Iline=I0(i,:);
%     I_signal_rough=Iline(Iline>200);
%     Signal_median=median(I_signal_rough);
%     Th=Signal_median/2;%threshold for judging a edge
%     signal_range=find(Iline>Th);
%     if ~isempty(signal_range)
%         Edge_left(i)=signal_range(1);
%         Edge_right(i)=signal_range(end);
%     else
%         Edge_left(i)=1;
%         Edge_right(i)=L;
%     end
% end
%-----find threshold for judging signal
I0_center=I0(round(N/2),:);
I0_center_smooth=movmean(I0_center,20);
thre=1/4*max(I0_center_smooth);
%--------------------------------------
% find edge
for i=1:N
    Iline=I0(i,:);
    %             I_signal_rough=Iline(Iline>200);
    I_signal_rough=Iline(Iline>thre);
    Signal_median=median(I_signal_rough);
    Th=Signal_median/2;%threshold for judging a edge
    signal_range=find(Iline>Th);
    if ~isempty(signal_range)
        Edge_left(i)=signal_range(1);
        Edge_right(i)=signal_range(end);
    else
        Edge_left(i)=1;
        Edge_right(i)=L;
    end
end
% find outliers
Edge_left_median=median(Edge_left);
Edge_right_median=median(Edge_right);
T=Edge_right-Edge_left;% thickness of sensor in unite of pixel
T_median=median(T);
outlier=find((abs(T-T_median)/T_median>0.2)==1);
d=(Calib.anode-Calib.cathode)*Calib.scale;%thickness of sensor
for i=1:N
    Eline=E(i,Edge_left(i):Edge_right(i));
    x=(0:T(i))/T(i)*d*1e-3;
    E_int(i)=trapz(x,Eline);
end
f2=figure(2);
movegui(f2,[612 575]);%movegui(f2,[612 903]);
plot(E_int)
xlabel('lateal pixel')
xlabel('Vertical integration of E-field [V]')
%%
waitfor(msgbox(['Select the region of interest']));
[points,y]=ginput(2);
if points(1)<1
    points(1)=1;
end
E_selected=Output.E_field_biased_corrected(points(1):points(2),:);%select a region the integation of E-field roughly equals to the bias applied
E_profile_selected=mean(E_selected,'omitnan');%
f3=figure(3);
movegui(f3,[40 66]);%movegui(f3,[40 392]);
imagesc(E_selected)
title('selected E-field sensor-')
axes1=gca;
set(axes1,'DataAspectRatio',[1 1 1],'Layer','top');
set(axes1,'View',[-90 90]);
colorbar
box


xs=(0:(Edge_right_median-Edge_left_median))/T_median*d*1e-3;
E_integral_selected=trapz(xs,E_profile_selected(Edge_left_median:Edge_right_median))%integral of E-field in the seleteced region

f4=figure(4);
movegui(f4,[612  66]);%movegui(f4,[612  392]);
hold on
% plot(Output.E_cross_section_average_corrected,'displayname','all')
plot(Calib.x_all-Calib.x_all(Calib.cathode),Output.E_cross_section_average_corrected,'displayname','all');
% plot(E_profile_selected,'displayname','selected')
plot(Calib.x_all-Calib.x_all(Calib.cathode),E_profile_selected,'displayname','selected');

legend
int_E_all=num2str(Output.integral_Efield);
int_E_selected=num2str(E_integral_selected);
title(['integration of E-field: all=' int_E_all 'V, selected=' int_E_selected 'V'])
% E_max=(1/Calib.alpha) * (asin(1));% max possible E-field calculated from light intensity
% E_max_vec=ones(1,L)*E_max;
% plot(E_max_vec,'k--','displayname','E-max')
xlabel('thickness (mm)')
ylabel('E field (V/m)')
box