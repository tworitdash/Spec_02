function [krho, f_axis] = finddrop(k0, epsilon_r, krho_g, h, zeta, zeta_s, mode)

    
       k_guess_norm = krho_g ./ k0; 
       
       c = 3*10^8;
       
       f0 = c * k0 /(2 * pi);
       
       f_axis = f0:-0.1e9:1e9;
       
       krho = zeros(size(f_axis, 1), size(f_axis, 2));
       
        for i = 1:size(f_axis, 2)
            
            k0_fi = (2 .* pi .* f_axis(i))./c;
            
            k_guess = k_guess_norm .* k0_fi;
           
            [D_kro_i] = D_kro(k0_fi, mode, epsilon_r, h, zeta, zeta_s, k_guess);
            
            
            
            del_k = k0_fi / 500;
            
            D_dkro = (D_kro(k0_fi, mode, epsilon_r, h, zeta, zeta_s, k_guess + del_k/2) - D_kro(k0_fi, mode, epsilon_r, h, zeta, zeta_s, k_guess - del_k/2)) ./ del_k;
            
            corr = D_kro_i./D_dkro;
            
            krho(i) = (k_guess - corr)/k0_fi;
            
            k_guess_norm = krho(i);
            
            if krho(i) < 1 || ~isreal(krho)
                krho(i) = NaN;
            end
            
        end
        
        
      
end