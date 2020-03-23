function [Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz] = Dyad(k0, kx, ky, kz)
    
    Dxx = k0.^2 - kx.^2;  Dxy = - kx .* ky;  Dxz = - kx .* kz;
    Dyx = - ky .* kx;  Dyy = k0.^2 - ky.^2;  Dyz = - ky .* kz;
    Dzx = - kz .* kx;  Dzy = - kz .* ky;  Dzz = k0.^2 - kz.^2;
    
   

end