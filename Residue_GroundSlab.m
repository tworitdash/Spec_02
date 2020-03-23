function [VtmR, ItmR] = Residue_GroundSlab(k0_, epsilon_r, h, k_sw, z)

    
    zeta = 120 .* pi;
    zeta_s = zeta ./ sqrt(epsilon_r);
    
    k_sub_ = k0_ .* sqrt(epsilon_r);
            
    kz0_ = -1j * sqrt(-(k0_.^2 - k_sw.^2));
    kzs_ = -1j * sqrt(-(k_sub_.^2 - k_sw.^2));
     
    z0_TM = zeta * kz0_ ./ k0_;
    zs_TM = zeta_s * kzs_ ./ k_sub_;
    
     
    zup = z0_TM;
    zdn = 1j .* zs_TM .* tan(kzs_ * h);
    
    del_k = k0_ ./ 500;
    
    mode = "TM";
    
    D_dkro_i = (D_kro(k0_, mode, epsilon_r, h, zeta, zeta_s, k_sw + del_k/2) - D_kro(k0_, mode, epsilon_r, h, zeta, zeta_s, k_sw - del_k/2)) ./ del_k;
     
    
    
    if z > 0 && z <= h
       
        
        VtmR = (zup .* zdn .* sin(kzs_ .* z)) ./ (D_dkro_i .* sin(kzs_ .* h));
        ItmR = (zup .* zdn .* 1j .* cos(kzs_ .* z)) ./ (D_dkro_i .* sin(kzs_ .* h) .* zs_TM);
        
    end
    
    if z > h
        
        VtmR = (zup .* zdn)./(D_dkro_i) .* exp(1j .* kz0_ .* h) .* exp(-1j .* kz0_ .* z);
        ItmR = VtmR ./ z0_TM;
        
    end

end