function [SGFxx, SGFxy, SGFxz, SGFyx, SGFyy, SGFyz, SGFzx, SGFzy, SGFzz] = SGF(Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz, c)
    
    
    
    SGFxx = c.* Dxx; SGFxy = c.* Dxy; SGFxz = c.* Dxz;
    SGFyx = c.* Dyx; SGFyy = c.* Dyy; SGFyz = c.* Dyz;
    SGFzx = c.* Dzx; SGFzy = c.* Dzy; SGFzz = c.* Dzz;

end