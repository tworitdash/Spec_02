function [zup_TE, zdn_TE, zup_TM, zdn_TM] = zupdn(z0_TE, zs_TE, z0_TM, zs_TM, h, kzs) 
    zup_TE = z0_TE;
    zdn_TE = 1j .* zs_TE .* tan(kzs * h);
    
    zup_TM = z0_TM;
    zdn_TM = 1j .* zs_TM .* tan(kzs * h);
    
    
end