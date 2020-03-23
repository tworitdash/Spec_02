function [D_inf] = D_inf_func_er(kx, th, ph, w, k0, dy, hd, er)

            %kx0 = k0 .* sin(th) .* cos(ph);
            zeta = 120 .* pi;
            ky0 = k0 .* sin(th) .* sin(ph);
            
            [kx0, my] = meshgrid(kx, -30:1:30);
            
            kym = ky0 - (2*pi*my)/dy;
            
            krho = sqrt(kx0.^2 + kym.^2);
            kz0 = (-1j).*sqrt(-(k0.^2 - krho.^2));
            
            z = hd + eps;
            
            zeta_s = zeta./sqrt(er); 
            
            
            ksub = k0 .* sqrt(er);
            
            kzs = (-1j)*sqrt(-(ksub.^2 - krho.^2));
            
            % For TM
            z0_TM = zeta * kz0 ./ k0;
            zs_TM = zeta_s * kzs ./ ksub;

            % For TE

            z0_TE = zeta * k0 ./ kz0;
            zs_TE = zeta_s * ksub ./ kzs;

            [zup_TE, zdn_TE, zup_TM, zdn_TM] = zupdn(z0_TE, zs_TE, z0_TM, zs_TM, hd, kzs);


            [vtm, vte, itm, ite] = txline(zup_TE, zdn_TE, zup_TM, zdn_TM, hd, z, kz0, kzs, z0_TE, zs_TE, z0_TM, zs_TM);


            [SGFxx, ~, ~] = Green(vtm, vte, itm, ite, kx0, kym, krho.^2, zeta, k0);


            
            %kzm = (-1j)*sqrt(-(k0.^2 - kx0.^2 - kym.^2));

       
            const = 1/dy;
            
            D_int = const .* SGFxx .* besselj(0, (kym .* w)/2); %.* (1 - exp(-1j .* kz0 .* 2 .* hd));
            
            D_inf = sum(D_int, 1);
            
end