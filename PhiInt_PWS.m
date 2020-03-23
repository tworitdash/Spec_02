function Iphi_pws = PhiInt_PWS(k0_, k_sw, l, w, epsilon_r)

    phi = eps:pi/180:2*pi;
    dphi = pi/180;
    
    kx = k_sw .* cos(phi);
    ky = k_sw .* sin(phi);
    k_eq = k0_ .* sqrt((1 + epsilon_r)./2);

    [Jx, Jy, Jz] = JFT_freq(l, w, k_eq, kx, ky);
    
    Iphi_pws_int = Jx.^2 .* (cos(phi)).^2 .* dphi;
    Iphi_pws = nansum(Iphi_pws_int);
    
end