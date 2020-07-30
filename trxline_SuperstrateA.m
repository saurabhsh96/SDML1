function [vTM, vTE, iTM, iTE] = trxline_SuperstrateA(k0,er,h,hs,kro,z)

%Characteristic impedance
zeta = 120*pi;
zetas = zeta/sqrt(er);

%kros = sqrt(er).*kro;

% Propagation along z
ks = sqrt(er)*k0;
kz = -1j*sqrt(-((k0^2)-(kro.^2)));
kzs = -1j*sqrt(-((ks^2)-(kro.^2)));

%Impedances:
Z0TE = (zeta.*k0)./kz;
Z0TM = (zeta.*kz)./k0;

ZsTE = (zetas.*ks)./kzs;
ZsTM = (zetas.*kzs)./ks;

ZLTE = ZsTE.*((Z0TE + (1j*ZsTE.*tan(kzs*hs)))./(ZsTE + (1j*Z0TE.*tan(kzs*hs))));
ZLTM = ZsTM.*((Z0TM + (1j*ZsTM.*tan(kzs*hs)))./(ZsTM + (1j*Z0TM.*tan(kzs*hs))));

Tau1TE = ((ZLTE - Z0TE).*exp(-2j.*kz.*h))./(ZLTE + Z0TE);
Tau1TM = ((ZLTM - Z0TM).*exp(-2j.*kz.*h))./(ZLTM + Z0TM);

V01TE = 1./(1+Tau1TE);
V01TM = 1./(1+Tau1TM);

Tau2TE = ((Z0TE - ZsTE).*exp(-2i.*kzs.*(h+hs)))./(Z0TE + ZsTE);
Tau2TM = ((Z0TM - ZsTM).*exp(-2i.*kzs.*(h+hs)))./(Z0TM + ZsTM);

VsTE = (V01TE.*(exp(-1i*kz.*h) + (Tau1TE.*exp(1i*kz.*h))))./(exp(-1i*kzs.*h) + (Tau2TE.*exp(1i*kzs.*h)));
VsTM = (V01TM.*(exp(-1i*kz.*h) + (Tau1TM.*exp(1i*kz.*h))))./(exp(-1i*kzs.*h) + (Tau2TM.*exp(1i*kzs.*h)));

V02TE = VsTE.*exp(1i.*kz.*(h+hs)).*(exp(-1i.*kzs.*(h+hs)) + (Tau2TE.*exp(1i.*kzs.*(h+hs))));
V02TM = VsTM.*exp(1i.*kz.*(h+hs)).*(exp(-1i.*kzs.*(h+hs)) + (Tau2TM.*exp(1i.*kzs.*(h+hs))));

if (z<=h)
    vTE = (V01TE.*(exp(-1i*kz.*z) + (Tau1TE.*exp(1i*kz.*z))));
    iTE = vTE./Z0TE;
    
    vTM = (V01TM.*(exp(-1i*kz.*z) + (Tau1TM.*exp(1i*kz.*z))));
    iTM = vTM./Z0TM;
    
elseif (z>h && z<h+hs)
        vTE = (VsTE.*(exp(-1i*kzs.*z) + (Tau2TE.*exp(1i*kzs.*z))));
        iTE = vTE./ZsTE;
        
        vTM = (VsTM.*(exp(-1i*kzs.*z) + (Tau2TM.*exp(1i*kzs.*z))));
        iTM = vTM./ZsTM;
        
else
    vTE = (V02TE.*exp(-1i.*kz.*z));
    iTE = (vTE./Z0TE);

    vTM = (V02TM.*exp(-1i.*kz.*z));
    iTM = (vTM./Z0TM);

end

end