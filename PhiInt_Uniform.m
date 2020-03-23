function Iphi_uni = PhiInt_Uniform(k0_, k_sw, l, w)

    phi = eps:pi/180:2*pi;
    dphi = pi/180;
    
    kx = k_sw .* cos(phi);
    ky = k_sw .* sin(phi);
   

    [Jx, Jy, Jz] = JFT_uni(l, w, kx, ky);
    
    Iphi_uni_int = Jx.^2 .* (cos(phi)).^2 .* dphi;
    Iphi_uni = nansum(Iphi_uni_int);
    
end