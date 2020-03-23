function [vte, vtm, ite, itm] = txline_Slot(z0_TE, z0_TM, er, k0, krho, hd, mode)
zeta = 120 * pi;
ksub = k0 .* sqrt(er);
kzs = -1j.*sqrt(-(ksub.^2 - krho.^2));
zeta_s = zeta ./ sqrt(er);
zs_TE = zeta_s .* ksub./kzs;
zs_TM = zeta_s .* kzs./ksub;

if mode == "PPW"
    
    zup_TE = 0;
    zup_TM = 0;
    zdn_TE = 1j .* zs_TE .* tan(kzs .* hd);
    zdn_TM = 1j .* zs_TM .* tan(kzs .* hd);
    
    vte = 1;
    vtm = 1;
    ite = 1./zdn_TE;
    itm = 1./zdn_TM;
    
elseif mode == "PEC"
    
    vte = 1;
    vtm = 1;
    ite = 1./z0_TE;
    itm = 1./z0_TM;
    
end

end