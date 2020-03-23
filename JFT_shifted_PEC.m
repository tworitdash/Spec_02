function [Jx, Jy, Jz] = JFT_shifted_PEC(l, w, h, k0, kx, ky, kz)
    
    T = sinc(ky*w/2/pi);
    L = 2 .* k0 .* (cos(kx.*l/2)) ./ ((k0.^2 - kx.^2)) .* (2 * 1j * sin(kz .* h)); %only valid for lambda/2 dipoles

    Jx = L.*T;      Jy = 0; Jz = 0;

end