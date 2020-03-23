function [Jx, Jy, Jz] = JFT_freq(l, w, k0, kx, ky, kz)
    
    T = sinc(ky.*w/2/pi);
    L = 2 .* k0 .* (cos(kx.*l/2) - cos(k0.*l/2)) ./ ((k0.^2 - kx.^2) .* sin(k0.*l/2));

    Jx = L.*T;      Jy = 0; Jz = 0;

end