function [Dir, P_rad] = Directivity(r_obs, th, ph, E_th, E_ph, zeta, dth, dph, k0)

    V_far_th = E_th * r_obs./(exp(-1j * k0 * r_obs));
    V_far_ph = E_ph * r_obs./(exp(-1j * k0 * r_obs));
    
    C_rad = 1/(2 * zeta);
    V_far_abs = sqrt(abs(V_far_th).^2 + abs(V_far_ph).^2);
    U = C_rad * (V_far_abs).^2;

    P_rad_int = C_rad * (V_far_abs).^2 .* sin(th) * dth * dph;

    P_rad = sum(sum(P_rad_int.'));
    
    Dir = 4 * pi * U ./ P_rad;
end
