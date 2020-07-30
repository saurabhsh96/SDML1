%Q2 Assign 1
clear;
%% Defining inputs

%Material
er = 12;
hs = 2.1e-3;
h = 15e-3;
loc = 0;

%EM
freq = [9e9 9.5e9 10e9 10.5e9 11e9];
c = 3e8;
lam = c./freq;
k0 = 2*pi./lam;
kx = zeros(size(freq));
res = 1000; %resolution
ky = zeros(size(freq, 1), res);
kRho = zeros(size(freq, 1), res);
for ind=1:size(freq, 2)
    ky(ind,:) = linspace(0, 1*k0(ind), res);
    kRho(ind,:) = sqrt(ky(ind,:).^2 + kx(ind).^2);
end

%Wave impedances
eps_0 = 8.854187817e-12;
mu_0 = 1.2566370614e-6;
zeta0 = (sqrt(mu_0/(eps_0*1)));
zetaS = (sqrt(mu_0/(eps_0*er)));

%Position
z = h+hs+.1;
%% Calling requ functions

vTM = zeros(size(freq, 2), size(ky, 2));
vTE = zeros(size(freq, 2), size(ky, 2));
iTM = zeros(size(freq, 2), size(ky, 2));
iTE = zeros(size(freq, 2), size(ky, 2));

Gxx = zeros(size(freq, 2), size(ky, 2));
Gyx = zeros(size(freq, 2), size(ky, 2));
Gzx = zeros(size(freq, 2), size(ky, 2));
Gej = zeros(size(freq, 2), size(ky, 2));

for ind=1:size(freq, 2)
    %[vTM(ind, :), vTE(ind, :), iTM(ind, :), iTE(ind, :)] = ...
    %    trxline_SuperStrate(k0(ind), er, h, hs, zeta0, zetaS, kRho(ind,:), z);
    [vTM(ind, :), vTE(ind, :), iTM(ind, :), iTE(ind, :)] = ...
        trxline_SuperstrateA(k0(ind),er,h,hs,kRho(ind,:),z);
    [Gxx(ind, :), Gyx(ind, :), Gzx(ind, :)] = ...
        SpectralGFem(k0(ind), er, kx(ind), ky(ind, :), vTM(ind, :), ...
        vTE(ind, :), iTM(ind, :), iTE(ind, :), zeta0, zetaS);
    Gej(ind, :) = abs(Gyx(ind,:));
end

for ind=1:size(freq,2)
    plot(ky(ind,:)./k0(ind), Gej(ind, :), 'LineWidth', 1.5, ...
        'DisplayName', ['freq = ', num2str(freq(ind)./10^9), ' GHz']);
    hold on;
end
title('Y-component of spectral field at z=h+h_s^+ w.r.t k_y');
xlabel('k_y/k_0');
ylabel('|G^{em}_{yx}| (in V)');
%ylim([0 4]);
legend show;