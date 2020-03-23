function [Jx, Jy, Jz] = J_Circular(d, k0, th)

    Jx = 0;
    Jy = pi * (d/2)^2 * besselj(1, (k0 * d/2 .* sin(th))) ./ (k0 * (d/2) .* sin(th));
    %Jy = pi * (d/2)^2 * besselj(1, (k0 * d/2 .* sin(th))) ./ (k0 * (d/2) .* sin(th));
    Jz = 0;
end