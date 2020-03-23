function [G_xx, G_yx, G_zx] = Green_hm(vtm, vte, itm, ite, kx, ky, krho, zeta, k0)

    G_xx = -(ite .* kx.^2 + itm .* ky.^2) ./ krho.^2;
    G_yx = (itm - ite).*kx.*ky ./ krho.^2;
    G_zx = kx./(zeta .* k0) .* vte;

end