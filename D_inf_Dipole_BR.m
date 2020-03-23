function D_inf = D_inf_Dipole_BR(kx, ky0, k0, dy, wd, my_lim, th, hd)
D_inf = zeros(size(th,1), size(th,2));
zeta = 120*pi;
for my = -my_lim:my_lim
    kym = ky0 - 2*pi*my/dy;
    kz = -1j.*sqrt(-(k0.^2-kx.^2-kym.^2));
    [SGFxx, ~, ~, ~, ~, ~, ~, ~, ~] = SGF(zeta, k0, kx, kym, kz);
    D_inf = D_inf + SGFxx.*(1-exp(-1j*kz*2*hd)).*besselj(0,kym*wd/2);
end
D_inf = D_inf/dy;
end

