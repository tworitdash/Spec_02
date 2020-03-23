function [Psw_TM] = PswTMelem(k0_, epsilon_r, h, k_sw, ItmR, z)

    dz = z(2) - z(1);
    zeta = 120 * pi;
    zeta_s = zeta ./ sqrt(epsilon_r);
    k_sub_ = k0_ .* sqrt(epsilon_r);
    
    Iz_int_s = zeta_s ./ k_sub_ .* abs(ItmR(z<=h)).^2 .* dz;
    
    Iz_s = nansum(Iz_int_s);
    
    Iz_int_0 = zeta ./ k0_ .* abs(ItmR(z>h)).^2 .* dz;
    
    Iz_0 = nansum(Iz_int_0);
    
    Iz = Iz_s + Iz_0;
    
    I_phi_del = pi;
    Psw_TM = 1/2 .* k_sw.^2 .* 1/(2 * pi) .* Iz .* I_phi_del;
    
end