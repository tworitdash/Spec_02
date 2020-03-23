function [D_inf] = D_inf_func(kx, th, ph, w, k0, dy, hd)

            %kx0 = k0 .* sin(th) .* cos(ph);
            zeta = 120 .* pi;
            ky0 = k0 .* sin(th) .* sin(ph);
            
            [kx0, my] = meshgrid(kx, -20:1:20);
            
            kym = ky0 - (2*pi*my)/dy;
            
            kzm = (-1j)*sqrt(-(k0.^2 - kx0.^2 - kym.^2));
    
            c = (-zeta./(2 * k0 * kzm)); %constant term in the equations

        %% Calculation of the Dyad
            [Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz] = Dyad(k0, kx0, kym, kzm);

        %% Calculation of Spectral Green's function

            [SGFxx, ~, ~, ~, ~, ~, ~, ~, ~] = SGF(Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz, c);

       
            const = 1/dy;
            
            D_int = const .* SGFxx .* besselj(0, (kym .* w)/2) .* (1 - exp(-1j .* kzm .* 2 .* hd));
            
            D_inf = sum(D_int, 1);
            
end