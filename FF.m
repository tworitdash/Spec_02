function [Ex, Ey, Ez] = FF(c2, SGFxx, SGFxy, SGFxz, SGFyx, SGFyy, SGFyz, SGFzx, SGFzy, SGFzz, Jx, Jy, Jz)
    
    Ex = c2 .* (SGFxx.*Jx + SGFxy.*Jy + SGFxz.*Jz);
    Ey = c2 .* (SGFyx.*Jx + SGFyy.*Jy + SGFyz.*Jz);
    Ez = c2 .* (SGFzx.*Jx + SGFzy.*Jy + SGFzz.*Jz);

end