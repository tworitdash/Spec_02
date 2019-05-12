function [G_xx, G_yx, G_zx] = Green_em(vtm, vte, itm, ite, kx, ky, k_rho_2, zeta, k0)

    G_xx = (vtm - vte).*kx.*ky./k_rho_2;
    G_yx = (vte .* kx.^2 + vtm .* ky.^2) ./ k_rho_2;
    G_zx = - zeta .* ky/k0 .* itm;

end