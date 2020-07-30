%TRX Superstrate
function [vTM, vTE, iTM, iTE] = trxline_SuperStrate3(k0, er, h, zeta0, zetaS, kRho, z)
    %In Slab
    ks = sqrt(er).*k0;
    
    %kZ
    kz0 = -1j.*sqrt(-((k0.^2)-(kRho.^2)));
    kzs = -1j.*sqrt(-((ks.^2)-(kRho.^2)));

    %Air TE TM impedance
    Z0TE = (zeta0.*k0)./kz0;
    Z0TM = (zeta0.*kz0)./k0;

    %Slab TE TM impedance
    ZsTE = (zetaS.*ks)./kzs;
    ZsTM = (zetaS.*kzs)./ks;
    
    %Gamma
    gammaTE = (ZsTE - Z0TE)./(ZsTE + Z0TE).*exp(-2.*1j.*kz0.*h);
    gammaTM = (ZsTM - Z0TM)./(ZsTM + Z0TM).*exp(-2.*1j.*kz0.*h);

    %Constants
    V0TE = 1./(1+gammaTE);
    V0TM = 1./(1+gammaTM);
    
    VsTE = V0TE.*(exp(-1j.*kz0.*h) + gammaTE.*exp(1j.*kz0.*h)).*exp(1j.*kzs.*h);
    VsTM = V0TM.*(exp(-1j.*kz0.*h) + gammaTM.*exp(1j.*kz0.*h)).*exp(1j.*kzs.*h); 
    
    kRhoConst = 1;
    if(z<h)
        vTM = V0TM.*(exp(-1j.*kz0.*z) + gammaTM.*exp(1j.*kz0.*z)).*kRhoConst;
        iTM = (V0TM./Z0TM).*(exp(-1j.*kz0.*z) - gammaTM.*exp(1j.*kz0.*z)).*kRhoConst;
        vTE = V0TE.*(exp(-1j.*kz0.*z) + gammaTE.*exp(1j.*kz0.*z)).*kRhoConst;
        iTE = (V0TE./Z0TE).*(exp(-1j.*kz0.*z) - gammaTE.*exp(1j.*kz0.*z)).*kRhoConst;
    else
        vTM = VsTM.*(exp(-1j.*kzs.*z)).*kRhoConst;
        iTM = (VsTM./ZsTM).*(exp(-1j.*kzs.*z)).*kRhoConst;
        vTE = VsTE.*(exp(-1j.*kzs.*z)).*kRhoConst;
        iTE = (VsTE./ZsTE).*(exp(-1j.*kzs.*z)).*kRhoConst;
    end
end