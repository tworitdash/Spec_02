function [D_kro] = D_kro_leak(k0_fi, mode, epsilon_r, h, hs, zeta, zeta_s, k_guess)

     
     k_sub_fi = k0_fi .* sqrt(epsilon_r);
            
     kz0_fi = 1j * sqrt(-(k0_fi.^2 - k_guess.^2));
     kzs_fi = 1j * sqrt(-(k_sub_fi.^2 - k_guess.^2));
     
     if mode == "TE"
         z0 = zeta * k0_fi ./ kz0_fi;
         zs = zeta_s * k_sub_fi ./ kzs_fi;
         zL = zs .* (z0 + 1j .* zs .* tan(kzs_fi .* hs)) ./ (zs + 1j .* z0 .* tan(kzs_fi .* hs));
         D_kro = zL + 1j .* z0 .* tan(kz0_fi .* h);
     else
          z0 = zeta * kz0_fi ./ k0_fi;
          zs = zeta_s * kzs_fi ./ k_sub_fi;
          zL = zs .* (z0 + 1j .* zs .* tan(kzs_fi .* hs)) ./ (zs + 1j .* z0 .* tan(kzs_fi .* hs));
          D_kro = zL + 1j .* z0 .* tan(kz0_fi .* h);
     end
     
    
   
     
     
end