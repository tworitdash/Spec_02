function [vte, vtm, ite, itm] = txline_3(Gamma_TE, Gamma_TM, h, kzs, kz0, z0_TE, z0_TM, zs_TE, zs_TM, z)
    % V_0^+ for TE
    v0p_TE = 1./(1 + Gamma_TE);
    % V_0^+ for TM
    v0p_TM = 1./(1 + Gamma_TM);
    
    % V_s^+ for TE
    vsp_TE = v0p_TE .* (exp(-1j .* kz0 .* h) + Gamma_TE .* exp(1j .* kz0 .* h)) .* exp(1j .* kzs .* h);
    % V_s^+ for TM
    vsp_TM = v0p_TM .* (exp(-1j .* kz0 .* h) + Gamma_TM .* exp(1j .* kz0 .* h)) .* exp(1j .* kzs .* h);

    if z > 0 && z < h
        %For TE
        vte = v0p_TE .* (exp(-1j .* kz0 .* z) + Gamma_TE .* exp(1j .* kz0 .* z));
        ite = v0p_TE ./ z0_TE .* (exp(-1j .* kz0 .* z) + Gamma_TE .* exp(1j .* kz0 .* z));
        % For TM
        vtm = v0p_TM .* (exp(-1j .* kz0 .* z) + Gamma_TM .* exp(1j .* kz0 .* z));
        itm = v0p_TM ./ z0_TM .* (exp(-1j .* kz0 .* z) + Gamma_TM .* exp(1j .* kz0 .* z));
        
    end
    
    if z > h
        % For TE
        vte = vsp_TE .* exp(-1j .* kzs .* z);
        ite = vsp_TE ./ zs_TE .* exp(-1j .* kzs .* z);
        
        % For TM
        vtm = vsp_TM .* exp(-1j .* kzs .* z);
        itm = vsp_TM ./ zs_TM .* exp(-1j .* kzs .* z);
    end

end