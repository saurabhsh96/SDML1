%GEM
function [Gxx, Gyx, Gzx] = SpectralGFem3(k0, er, kx, ky, vTM, vTE, iTM, iTE, zeta0, zetaS)
    kRho = sqrt(kx.^2 + ky.^2);
    Gxx = (vTM - vTE).*kx.*ky./(kRho.^2); 
    Gyx = (vTE.*kx.^2 + vTM.*ky.^2)./(kRho.^2);
    Gzx = -zeta0.*ky.*iTM./k0;
end 