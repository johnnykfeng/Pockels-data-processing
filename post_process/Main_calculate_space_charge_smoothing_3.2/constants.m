%% ---------contants----
T = 300; %[K]
C.eps_0 = 8.8541878128E-12; % [F⋅m−1]permittivity of vacumm
C.q = 1.602176634e-19; %[C]
C.k = 1.380649e-23; % [J/K]
C.h = 6.62607015E-34; %[J.s]
C.c = 299792458; %[m/s]
C.m0 = 9.1093837015e-31; %[kg] mass of electron
C.pi = 3.141592653589793;
C.h_bar = 1.054571817e-34; %[J.s]
C.Vt = C.k * T / C.q; %0.0259;%[V] thermal voltage
C.NA = 6.02214076E23; %[mol^−1] Avogadro constant
%% material parameters (HF CZT)
C.epsilon = 10.9;
C.Eg = 1.57; %[eV] bandgap
C.me = 0.11; %effective mass of electron
C.mh = 0.73; %effective mass of hole
C.mue = 940; %mobility of electron
C.muh = 114; %mobility of hole
C.NC = 2 * (2 * pi * C.me * C.m0 * C.k .* T / (C.h^2)).^1.5; %[m-3]density of states in conduction band
C.NV = 2 * (2 * pi * C.mh * C.m0 * C.k .* T / (C.h^2)).^1.5; %[m-3]density of states in valence band
C.rho = 5.8E3; %[kg/m^3] density
C.ma = 117.6539; %[u] atomic mass of CZT
C.mw = 0.23531; %[kg/mol] atomic mass of CZT
C.V_th_n = sqrt(C.k * T / (C.me * C.m0)); %[m/s] Thermal velocity of electron
C.V_th_p = sqrt(C.k * T / (C.mh * C.m0)); %[m/s] Thermal velocity of hole
