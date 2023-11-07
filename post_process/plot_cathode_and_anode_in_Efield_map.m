% plot cathode and anode positions in E-field map
% Calib.mat is loaded and E-field map figure is open
hold on
plot([Calib.cathode,Calib.cathode],[0,696],'--w')
plot([Calib.anode,Calib.anode],[0,696],'--w')