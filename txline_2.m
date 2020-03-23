function [vte, vtm, ite, itm] = txline_2(Gamma1_TE, Gamma1_TM, Gamma2_TE, Gamma2_TM, h, hs, kzs, kz0, zin_TE, zin_TM, zL_TE, zL_TM, z02_TE, z02_TM, z);

v01p_TE = 1 ./ (1 + Gamma1_TE);
v01p_TM = 1 ./ (1 + Gamma1_TM);

vsp_TE = v01p_TE .* (exp(-1j .* kz0 .* h) + Gamma1_TE .* exp(1j .* kz0 .* h)) ./ (exp(-1j .* kzs .* h) + Gamma2_TE .* exp(1j .* kzs .* h));
vsp_TM = v01p_TM .* (exp(-1j .* kz0 .* h) + Gamma1_TM .* exp(1j .* kz0 .* h)) ./ (exp(-1j .* kzs .* h) + Gamma2_TM .* exp(1j .* kzs .* h));

v02p_TE = vsp_TE .* (exp(-1j .* kzs .* (h + hs)) + Gamma2_TE .* (exp(1j .* kzs .* (h + hs)))) .* exp(1j .* kz0 .* (h + hs));
v02p_TM = vsp_TM .* (exp(-1j .* kzs .* (h + hs)) + Gamma2_TM .* (exp(1j .* kzs .* (h + hs)))) .* exp(1j .* kz0 .* (h + hs));


if z > 0 && z < h
    %For TE
    
    vte = v01p_TE .* (exp(-1j .* kz0 .* z) + Gamma1_TE .* exp(1j .* kz0 .* z));
    ite = v01p_TE ./ zin_TE .* (exp(-1j .* kz0 .* z) - Gamma1_TE .* exp(1j .* kz0 .* z));
    
    %For TM
    
    vtm = v01p_TM .* (exp(-1j .* kz0 .* z) + Gamma1_TM .* exp(1j .* kz0 .* z));
    itm = v01p_TM ./ zin_TM .* (exp(-1j .* kz0 .* z) - Gamma1_TM .* exp(1j .* kz0 .* z));
    
end

if z > h && z < (h + hs)
    %For TE
    vte = vsp_TE .* (exp(-1j .* kzs .* z) + Gamma2_TE .* exp(1j .* kzs .* z));
    ite = vsp_TE ./ zL_TE .* (exp(-1j .* kzs .* z) - Gamma2_TE .* exp(1j .* kzs .* z));
    
    %FOr TM
    vtm = vsp_TM .* (exp(-1j .* kzs .* z) + Gamma2_TM .* exp(1j .* kzs .* z));
    itm = vsp_TM ./ zL_TM .* (exp(-1j .* kzs .* z) - Gamma2_TM .* exp(1j .* kzs .* z));
end

if z > (h + hs)
    %For TE
    vte = v02p_TE .* exp(-1j .* kz0 .* z);
    ite = v02p_TE ./ z02_TE .* exp(-1j .* kz0 .* z);
    
    %For TM
    vtm = v02p_TM .* exp(-1j .* kz0 .* z);
    itm = v02p_TM ./ z02_TM .* exp(-1j .* kz0 .* z);
end

end