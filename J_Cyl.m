function [J_r, J_theta, J_phi] = J_Cyl(D, k0, th, ph)

    J_r = pi * (D/2)^2 * besselj(1, (k0 * D/2 .* sin(th))) ./ (k0 * (D/2) .* sin(th)) .* sin(th) .* sin(ph);
    J_theta  = pi * (D/2)^2 * besselj(1, (k0 * D/2 .* sin(th))) ./ (k0 * (D/2) .* sin(th)) .* cos(th) .* sin(ph);
    J_phi = pi * (D/2)^2 * besselj(1, (k0 * D/2 .* sin(th))) ./ (k0 * (D/2) .* sin(th)) .* cos(ph);
end