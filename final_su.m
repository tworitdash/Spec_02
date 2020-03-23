%% Antenna dimesnions

% for dipole
clear all;
close all;

c0 = 3e8;
dth = pi/180;
dph = pi/180;
er = 1;

Nx = 7;
Dummy = 2;

%f = linspace(2e9, 15e9, 1000);
%k = (2*pi*f)./c0;
f0 = 10e9;

lambda0 = c0./f0;
lambda_ = lambda0 ./ sqrt(er);



w = 0.2 .* lambda0;
deld = 0.25 .* lambda0;
dx =  0.5 * lambda0;
dy =  0.5 * lambda0;
hd = 0.25 .* lambda_;
    
    
[mx, my] = meshgrid(-30:1:30, -30:1:30);
%[th, ph] = meshgrid(eps:dth:pi/2, eps:dph:2*pi);
    
    th0 = pi/6;
    ph0 = 1e-7;
    
    k0 = 2*pi/lambda0;
    
    kx = k0.* sin(th0) * cos(ph0);
    ky = k0.* sin(th0) * sin(ph0);
    
    
    %kz = sqrt(k0.^2 - krho.^2);
    
    kxm = k0 - (2.*pi*mx)/dx;
    kym = k0 - (2.*pi*my)/dy;
    
    krho = sqrt(kxm.^2 + kym.^2);
    
    kz0 = (-1j).*sqrt(-(k0.^2 - krho.^2));
            
            z = hd+eps;
            zeta = 120*pi;
            
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


            [SGFxx, ~, ~] = Green(vtm, vte, itm, ite, kxm, kym, krho.^2, zeta, k0);


            
            %kzm = (-1j)*sqrt(-(k0.^2 - kx0.^2 - kym.^2));

       
            const = 1/dy;
            
            D_int = const .* SGFxx .* besselj(0, (kym .* w)/2); %.* (1 - exp(-1j .* kz0 .* 2 .* hd));
            
            D_inf = sum(D_int, 1);
    
    
    %D_inf = D_inf_func_er(kxm(1, :), th0, ph0, w, k0, dy, hd, er); 
            
    
 figure(2);
 plot(kxm(1, :)/k0, abs(1./D_inf), 'LineWidth', 2);
 grid on;
% 
% xlabel('Number of active elements', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Angle of scan blindness (\theta_{blind}) on \phi = 0 plane', 'FontSize', 12, 'FontWeight', 'bold');
% title('Angle of scan blindness vs Number of elements', 'FontSize', 12, 'FontWeight', 'bold');



