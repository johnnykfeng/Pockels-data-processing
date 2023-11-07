function [rho]=Func_rho_from_divE(E_div)
C.eps_0=8.8541878128E-12;% [F⋅m−1]permittivity of vacumm
C.q=1.602176634e-19; %[C]
C.epsilon=10.9;
eps_CZT=C.eps_0*C.epsilon;


rho=-E_div./(10^-5)*eps_CZT./10^6/C.q; %[e/cm^3];
end