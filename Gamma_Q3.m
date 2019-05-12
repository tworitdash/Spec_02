function [Gamma_TE, Gamma_TM] = Gamma_Q3(zs_TE, zs_TM, z0_TE, z0_TM, h, kz0)

    % For TE
    Gamma_TE = (zs_TE - z0_TE)./(z0_TE + zs_TE) .* exp(-2j .* kz0 .* h);
    Gamma_TM = (zs_TM - z0_TM)./(z0_TM + zs_TM) .* exp(-2j .* kz0 .* h);
    

end