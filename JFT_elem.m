function [Jx_elem, Jy_elem, Jz_elem] = JFT_elem(l, w, kx, ky)
    
    Jx_elem = l .* sinc(kx.*l/2/pi) .* sinc(ky.*w/2/pi);
    Jy_elem = 0;
    Jz_elem = 0;

end