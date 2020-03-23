function [Gamma1_TE, Gamma1_TM, Gamma2_TE, Gamma2_TM] = Gamma(zs_TE, zs_TM, zL_TE, zL_TM, z0_TE, z0_TM, h, hs, kz0, kzs)

%Gamma1
Gamma1_TE = (zL_TE - z0_TE) ./ (zL_TE + z0_TE) .* exp(-2j .* kz0 .* h);
Gamma1_TM = (zL_TM - z0_TM) ./ (zL_TM + z0_TM) .* exp(-2j .* kz0 .* h);

%Gamma2
Gamma2_TE = (z0_TE - zs_TE) ./ (z0_TE + zs_TE) .* exp(-2j .* kzs .* (h + hs));
Gamma2_TM = (z0_TM - zs_TM) ./ (z0_TM + zs_TM) .* exp(-2j .* kzs .* (h + hs));

end