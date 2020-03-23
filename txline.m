function [vtm, vte, itm, ite] =  txline(zup_TE, zdn_TE, zup_TM, zdn_TM, h, z, kz0, kzs, z0_TE, zs_TE, z0_TM, zs_TM)

    if z > h
       % for TE
       v0_TE_p = exp(1j .* kz0 .* h) .* (zup_TE .* zdn_TE)./(zup_TE + zdn_TE);
       vte = v0_TE_p .* exp(-1j .* kz0 .* z);
       ite = vte ./ z0_TE;
       
       % for TM
       v0_TM_p = exp(1j .* kz0 .* h) .* (zup_TM .* zdn_TM)./(zup_TM + zdn_TM);
       vtm = v0_TM_p .* exp(-1j .* kz0 .* z);
       itm = vtm ./ z0_TM;
       
       
    end
    
    if z < h && z >= 0
        
        % for TE 
        Vs_TE_p = (zup_TE * zdn_TE) / ((zup_TE + zdn_TE) * (exp(1j * kzs * h) - exp(-1j * kzs * h)));
        vte = Vs_TE_p * (exp(1j .* kzs) + exp(-1j * kzs));
        ite = Vs_TE_p / zs_TE * (exp(1j * kzs) - exp(-1j * kzs));
        % for TM
        Vs_TM_p = (zup_TM * zdn_TM) / ((zup_TM + zdn_TM) * (exp(1j * kzs * h) - exp(-1j * kzs * h)));
        vtm = Vs_TM_p * (exp(1j * kzs) + exp(-1j * kzs));
        itm = Vs_TM_p / zs_TM * (exp(1j * kzs) - exp(-1j * kzs));
        
    end
end