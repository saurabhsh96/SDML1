%TRX Line for ground slab
function [vTM, vTE, iTM, iTE] = trxline_GroundSlab(k0, er, h, zeta0, zetaS, kRho, z)
    
    %In Slab
    ks = sqrt(er)*k0;
    
    %kZ
    kz0 = -1j*sqrt(-((k0^2)-(kRho.^2)));
    kzs = -1j*sqrt(-((ks^2)-(kRho.^2)));

    %Air TE TM impedance
    Z0TE = (zeta0.*k0)./kz0;
    Z0TM = (zeta0.*kz0)./k0;

    %Slab TE TM impedance
    ZsTE = (zetaS.*ks)./kzs;
    ZsTM = (zetaS.*kzs)./ks;

    %Zup and Zdown
    ZupTE = Z0TE;
    ZupTM = Z0TM;

    ZdownTE = 1j*ZsTE.*tan(kzs.*h);
    ZdownTM = 1j*ZsTM.*tan(kzs.*h);

    %Term in the Slab Voltage and current
    %Voltage
    TE_Constant = (ZupTE.*ZdownTE)./(ZupTE + ZdownTE);
    TM_Constant = (ZupTM.*ZdownTM)./(ZupTM + ZdownTM);

    %Current In
    TE_inSlab = TE_Constant./ZsTE;
    TM_inSlab = TM_Constant./ZsTM;

    %Current Out
    TE_Air = TE_Constant./Z0TE;
    TM_Air = TM_Constant./Z0TM;

    const = 1;
    if (z < h) %Slab
        vTM = (TM_Constant.*sin(kzs.*z))./(sin(kzs.*h)).*const;
        iTM = (TM_inSlab.*1j.*cos(kzs.*z))./(sin(kzs.*h)).*const;
        vTE = (TE_Constant.*sin(kzs.*z))./(sin(kzs.*h)).*const;
        iTE = (TE_inSlab.*1i.*cos(kzs.*z))./(sin(kzs.*h)).*const;    
    else  %Air
        vTM = TM_Constant.*exp(1j*kz0.*h).*exp(-1j.*kz0.*z).*const;
        iTM = TM_Air.*exp(1j*kz0.*h).*exp(-1j.*kz0.*z).*const;
        vTE = TE_Constant.*exp(1j*kz0.*h).*exp(-1j.*kz0.*z).*const;
        iTE = TE_Air.*exp(1j*kz0.*h).*exp(-1j.*kz0.*z).*const;
    end
end