function [Erho, Ez, Hphi] = TMSwFields(k0_, k_sw, epsilon_r, VtmR, ItmR, rho, phi, z, h)
    k_eq = k0_ .* sqrt((1 + epsilon_r)./2);
    
    l = pi/k0_;
    
    w = l ./ 10;
    
    zeta = 120 .* pi;
    zeta_s = zeta ./ sqrt(epsilon_r);
    k_sub_ = k0_ .* sqrt(epsilon_r);
    
    C = 1j .* sqrt(k_sw ./ (2 .* pi)) .* exp(1j .* pi ./ 4);
    kx = k_sw .* cos(phi);
    ky = k_sw .* sin(phi);
    
    
    
    [Jx, Jy, Jz] = JFT_freq(l, w, k_eq, kx, ky);
    
    Erho = VtmR .* Jx .* C .* cos(phi) .* exp(-1j .* k_sw .* rho) ./ sqrt(rho);
    if z < h && z > 0
        Ez = -(zeta_s .* k_sw ./k_sub_) .* ItmR .* Jx .* C .* cos(phi) .* exp(-1j .* k_sw .* rho) ./ sqrt(rho);
    else
        Ez = -(zeta .* k_sw ./k0_) .* ItmR .* Jx .* C .* cos(phi) .* exp(-1j .* k_sw .* rho) ./ sqrt(rho);
    end
    
    Hphi = ItmR .* Jx .* C .* cos(phi) .* exp(-1j .* k_sw .* rho) ./ sqrt(rho);

end