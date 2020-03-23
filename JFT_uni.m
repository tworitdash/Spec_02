function [Jx_uni, Jy_uni, Jz_uni] = JFT_uni(l, w, kx, ky)
    
    Jx_uni = l .* sinc(kx.*l/2/pi) .* sinc(ky.*w/2/pi);
    Jy_uni = 0;
    Jz_uni = 0;

end