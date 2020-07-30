%Q3 Infinite medium
clear;

%% Inputs

%Material
er = [2.5 4 12];
h = 5e-3;

%EM
freq = 30e9;
c = 3e8;
lam = c/freq;
k0 = 2*pi/lam;
kx = linspace(0, 2*k0, 1500);
ky = 0;
kRho = sqrt(kx.^2 + ky.^2);

%Wave impedances
eps_0 = 8.854187817e-12;
mu_0 = 1.2566370614e-6;
zeta0 = round(sqrt(mu_0/(eps_0*1)));
zetaS = round(sqrt(mu_0./(eps_0*er)));

%Position
z = h;

%% Calling req functions

vTM = zeros(size(er, 2), size(kx, 2));
vTE = zeros(size(er, 2), size(kx, 2));
iTM = zeros(size(er, 2), size(kx, 2));
iTE = zeros(size(er, 2), size(kx, 2));

Gxx = zeros(size(er, 2), size(kx, 2));
Gyx = zeros(size(er, 2), size(kx, 2));
Gzx = zeros(size(er, 2), size(kx, 2));
Gej = zeros(size(er, 2), size(kx, 2));

for ind=1:size(er, 2)
    [vTM(ind, :), vTE(ind, :), iTM(ind, :), iTE(ind, :)] = ...
        trxline_SuperStrate3(k0, er(ind), h, zeta0, zetaS(ind), kRho, z);
    [Gxx(ind, :), Gyx(ind, :), Gzx(ind, :)] = ...
        SpectralGFem(k0, er(ind), kx, ky, vTM(ind, :), ...
        vTE(ind, :), iTM(ind, :), iTE(ind, :), zeta0, zetaS(ind));
    Gej(ind, :) = abs(Gyx(ind,:));
end

for ind=1:size(er,2)
    plot(kx./k0, Gej(ind, :), 'LineWidth', 1.5, ...
        'DisplayName', ['\epsilon_r = ', num2str(er(ind))]);
    hold on;
end
title('Y-component of spectral field at z=h^+ w.r.t k_x');
xlabel('k_x/k_0');
ylabel('|G^{em}_{yx}| (in V)');
legend show;