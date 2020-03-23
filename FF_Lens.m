function [Eth, Ephi] = FF_Lens(c2, order, theta, phi)

    Eth = c2 .* cos(theta).^order .* cos(phi);
    Ephi = -c2 .* cos(theta).^order .* sin(phi);

end