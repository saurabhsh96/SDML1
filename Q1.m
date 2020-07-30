%Spectral Greens function for stratified media
clear;
%% Defining Inputs

%Material
er = 10;
h = 2e-3;

%EM
freq = 10e9;
freq1 = 20e9;

c = 3e8;
lam = c/freq;
lam1 = c/freq1;

k0 = 2*pi/lam;
k01 = 2*pi/lam1;

kx = 0:0.5:5*k0;
kx1 = 0:0.5:5*k01;
ky = 0;

kRho = sqrt(kx.^2 + ky.^2);
kRho1 = sqrt(kx1.^2 + ky.^2);

%Wave impedances
eps_0 = 8.854187817e-12;
mu_0 = 1.2566370614e-6;
zeta0 = (sqrt(mu_0/(eps_0*1)));
zetaS = (sqrt(mu_0/(eps_0*er)));

%Observation point
z = h + 0.0001;

%% Calling TRXline function and SGF
% Solution of equivalent transmission line
[vTM, vTE, iTM, iTE] = trxline_GroundSlab(k0, er, h, zeta0, zetaS, kRho, z);

[vTM1, vTE1, iTM1, iTE1] = trxline_GroundSlab(k01, er, h, zeta0, zetaS, kRho1, z);

% Dyadic Green's function - x component 
[Gxx, Gyx, Gzx] = SGFej(k0, er, kx, ky, vTM, vTE, iTM, iTE, zeta0, zetaS);
Gej = (abs(Gxx));
%Gej = abs(sqrt( Gxx.^2 + Gyx.^2 + Gzx.^2));

[Gxx1, Gyx1, Gzx1] = SGFej(k01, er, kx1, ky, vTM1, vTE1, iTM1, iTE1, zeta0, zetaS);
Gej1 = (abs(Gxx1));
%Gej1 = abs(sqrt( Gxx1.^2 + Gyx1.^2 + Gzx1.^2));

figure(1);
titl = 'X-component of spectral field at z = h^{+} w.r.t. k_x';
plot(kx./k0, Gej, 'LineWidth', 1.5, 'DisplayName', '10 GHz');
hold on;
plot(kx1./k01, Gej1, 'LineWidth', 1.5, 'DisplayName', '20 GHz');
legend show;
grid on;
xlabel('k_x/k_0');
ylabel('|G_{xx}^{ej}| (in V)');
title(titl);

%lim = [0 1300];
%name = 'Gxx';
%plotReq(kx./k0, Gej, 'k_x/k_0', '|Gxx|', titl, lim, name);