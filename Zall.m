function [zin_TE, zin_TM, zL_TE, zL_TM, z02_TE, z02_TM] = Zall(z0_TE, zs_TE, z0_TM, zs_TM, h, hs, kzs, kz0)

% Impedance at z = h + hs
z02_TE = z0_TE;
z02_TM = z0_TM;

% Impedance at z = h
zL_TE = zs_TE .* (z02_TE + 1j .* zs_TE .* tan(kzs .* hs)) ./ (zs_TE + 1j .* z02_TE .* tan(kzs .* hs));
zL_TM = zs_TM .* (z02_TM + 1j .* zs_TM .* tan(kzs .* hs)) ./ (zs_TM + 1j .* z02_TM .* tan(kzs .* hs));

%Impedance at z = 0
z01_TE = z0_TE;
z01_TM = z0_TM;

zin_TE = z01_TE .* (zL_TE + 1j .* z01_TE .* tan(kz0 .* h)) ./ (z01_TE + 1j .* zL_TE .* tan(kz0 .* h)) ;
zin_TM = z01_TM .* (zL_TE + 1j .* z01_TM .* tan(kz0 .* h)) ./ (z01_TM + 1j .* zL_TM .* tan(kz0 .* h));

end