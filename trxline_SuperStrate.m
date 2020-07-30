%TRX Superstrate
function [vTM, vTE, iTM, iTE] = trxline_SuperStrate(k0, er, h, hs, zeta0, zetaS, kRho, z)
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
    
    %ZL from voltage
    ZlTM = ZsTM.*(Z0TM + 1j.*ZsTM.*tan(kzs.*hs))./...
        (ZsTM + 1j.*Z0TM.*tan(kzs.*hs));
    ZlTE = ZsTE.*(Z0TE + 1j.*ZsTE.*tan(kzs.*hs))./...
        (ZsTE + 1j.*Z0TE.*tan(kzs.*hs));
    
    %Gamma
    gamma1TE = (ZlTE - Z0TE)./(ZlTE + Z0TE).*exp(-2.*1j.*kz0.*h);
    gamma1TM = (ZlTM - Z0TM)./(ZlTM + Z0TM).*exp(-2.*1j.*kz0.*h);
    
    gamma2TE = (Z0TE - ZsTE)./(Z0TE + ZsTE).*exp(-2.*1j.*kzs.*(h+hs));
    gamma2TM = (Z0TM - ZsTM)./(Z0TM + ZsTM).*exp(-2.*1j.*kzs.*(h+hs));
    
    %Constants
    %Substrate and Air
    ASTM = 1./(1+gamma1TM);
    ASTE = 1./(1+gamma1TE);
    
    slabTM = ASTM.*(exp(-1j.*kz0.*h) + gamma1TM.*exp(1j.*kz0.*h))./...
        (exp(-1j.*kzs.*h) + gamma2TM.*exp(1j.*kzs.*h));
    slabTE = ASTE.*(exp(-1j.*kz0.*h) + gamma1TE.*exp(1j.*kz0.*h))./...
        (exp(-1j.*kzs.*h) + gamma2TE.*exp(1j.*kzs.*h));
    
    airTM = slabTM.*exp(1j.*kz0.*(h+hs)).*(exp(-1j.*kzs.*(h+hs))...
        + gamma2TM.*exp(1j.*kzs.*(h+hs)));
    airTE = slabTE.*exp(1j.*kz0.*(h+hs)).*(exp(-1j.*kzs.*(h+hs))...
        + gamma2TE.*exp(1j.*kzs.*(h+hs)));
    
    %kRhoConst = (1j.*kRho)./(1j.*kRho);
    kRhoConst = 1;
    %Why no negative sign?
    if (z <= h) %between ground and slab
        vTM = ASTM.*(exp(-1j.*kz0.*z) + gamma1TM.*exp(1j.*kz0.*z)).*kRhoConst;
        iTM = (ASTM./Z0TM).*(exp(-1j.*kz0.*z) - gamma1TM.*exp(1j.*kz0.*z)).*kRhoConst;
        vTE = ASTE.*(exp(-1j.*kz0.*z) + gamma1TE.*exp(1j.*kz0.*z)).*kRhoConst;
        iTE = (ASTE./Z0TE).*(exp(-1j.*kz0.*z) - gamma1TE.*exp(1j.*kz0.*z)).*kRhoConst;    
    
    elseif (z > h & z<=h+hs) %Air
        vTM = slabTM.*(exp(-1j.*kzs.*z) + gamma2TM.*exp(1j.*kzs.*z)).*kRhoConst;
        iTM = (slabTM./ZsTM).*(exp(-1j.*kzs.*z) - gamma2TM.*exp(1j.*kzs.*z)).*kRhoConst;
        vTE = slabTE.*(exp(-1j.*kzs.*z) + gamma2TE.*exp(1j.*kzs.*z)).*kRhoConst;
        iTE = (slabTE./ZsTE).*(exp(-1j.*kzs.*z) - gamma2TE.*exp(1j.*kzs.*z)).*kRhoConst;
    
    else 
        vTM = airTM.*(exp(-1j.*kz0.*z)).*kRhoConst;
        iTM = (airTM./Z0TM).*(exp(-1j.*kz0.*z)).*kRhoConst;
        vTE = airTE.*(exp(-1j.*kz0.*z)).*kRhoConst;
        iTE = (airTE./Z0TE).*(exp(-1j.*kz0.*z)).*kRhoConst;
    
    end
end