function [zin_TE, zin_TM] = Zin(z0_TE, zs_TE, z0_TM, zs_TM, h, kz0)
    % For TE
   zin_TE = z0_TE .* (zs_TE + 1j .* z0_TE .* tan(kz0 .* h)) ./ (z0_TE + 1j .* zs_TE .* tan(kz0 .* h));
   zin_TM = z0_TM .* (zs_TM + 1j .* z0_TM .* tan(kz0 .* h)) ./ (z0_TM + 1j .* zs_TM .* tan(kz0 .* h));

end