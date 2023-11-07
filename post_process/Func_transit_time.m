function [t_trans]=Func_transit_time(Output,Calib)
%  -----------------------------------------------------------------------
% | Calculate the transit time of electrons from cathode to anode based on|
% | E-field profile                                                       |
% | 2021-07-29 Yuxin created code                                         |
%  -----------------------------------------------------------------------
%-------inputs------
mue=940;%[cm2/(Vs)] mobility of electron
smooth_step=10;%smooth step size for smoothing the E-field-profile
% thickness=2;%[mm] thickness of sensor
thickness=Output.thickness;%[mm] thickness of sensor
%-------------------

[E_clean,Edge_cathode,Edge_anode]=Func_clean_E_field(Output,Calib);
%------------------------------------------
% version=str2double(Output.version);
% if version<2
    scale=thickness/median(Edge_anode-Edge_cathode);%[mm/pixel] real thickness of one pixel in E-field map
% else
%     scale=Calib.scale;
% end
%------------------------------------------
[W,~]=size(E_clean);
for i=1:W
    E=E_clean(i,(Edge_cathode(i):Edge_anode(i)));
    E_smooth=movmean(E,smooth_step);
    x=[1:length(E)]*scale;%[mm]
    v=mue*1e-4*E_smooth;%[m/s]
    t=scale*1e-3./v;%[s]
    t_trans(i)=sum(t,'omitnan');%[s]transit time along the E-field profile
%     %------------------
%     figure(11)
%     hold on
%     plot(x,E)
%     plot(x,E_smooth)
%     %------------------
end

E_ideal=Output.bias/(thickness*1e-3); % ideal uniform E-field based on bias and thickness
t_transist_ideal=thickness*1e-3/(mue*1e-4*E_ideal); % ideal transit time
t_transist_ideal2=ones(1,W)*t_transist_ideal;
y=[1:W]*scale;%[mm]
figure, hold on
plot(y,t_trans)
plot(y,t_transist_ideal2,'k--')
xlabel('Width [mm]')
ylabel('Transit time [s]')
box
grid










end