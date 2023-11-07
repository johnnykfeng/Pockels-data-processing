% calculate space charge density based on two points from E-field profile
% coordinate format [x(mm) y(100V/mm)]
% (100V/mm)=(1E5 V/m)
% point1 --> point2, left--> right
%---------input-----------------
%point=[  x     E-field]
point1=[1.05797 2.13702];
point2=[1.23188 1.92574];
%-------------------------------
E_div=(point2(2)-point1(2))*1e5./((point2(1)-point1(1))*1e-3);

C.eps_0=8.8541878128E-12;% [F⋅m−1]permittivity of vacumm
C.q=1.602176634e-19; %[C]
C.epsilon=10.9;
eps_CZT=C.eps_0*C.epsilon;


% rho=-E_div./(10^-5)*eps_CZT./10^6/C.q; %[e/cm^3];
rho=-E_div.*eps_CZT./10^6/C.q %[e/cm^3];
