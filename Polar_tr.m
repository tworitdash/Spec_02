function [E_r, E_th, E_ph] = Polar_tr(E_abs, th, ph)

    E_r = E_abs .* sin(th) .* sin(ph);
    E_th = E_abs .* cos(th) .* sin(ph);
    E_ph = E_abs .* cos(ph);


end