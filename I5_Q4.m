%% Antenna dimesnions

% for dipole
clear all;
close all;

c0 = 3e8;
dth = pi/180;
dph = pi/180;
er = 2;

f = 2e9:0.05e9:15e9;

lambda0 =c0./10^10;

w = 0.2 .* lambda0;
deld = 0.25 .* lambda0;
dx = 0.5 .* lambda0;
dy = 0.5 .* lambda0;
h = 0.31 .* lambda0;



th = pi/4;
ph = eps; % For E plane; change it to 90 for H plane

for j = 1:size(th, 2)

y_in_f = zeros(size(f, 1), size(f, 2));
z_in_f = zeros(size(f, 1), size(f, 2));
yin_higher = zeros(size(f, 1), size(f, 2));
y00 = zeros(size(f, 1), size(f, 2));
zin = zeros(size(f, 1), size(f, 2));
yin = zeros(size(f, 1), size(f, 2));
z00 = zeros(size(f, 1), size(f, 2));
Gamma = zeros(size(f, 1), size(f, 2));

%% wave properties 

    %f = 10e9;


   
%% Spectral Green's function
    

for p=1:size(f, 2)
    
    [mx, my] = meshgrid(-20:1:20, -30:1:30);
    %D_inf = zeros(size(mx));
    
 
     
        lambda = 3e8./f(p);
    
        k0 = 2 * pi ./ lambda;
    
        zeta = 120*pi;
   

  

            kx0 = k0 .* sin(th(j)) .* cos(ph);
            ky0 = k0 .* sin(th(j)) .* sin(ph);
            %kz0 = k0 .* cos(th(j));

            kxm = kx0 - (2*pi*mx)/dx;
            kym = ky0 - (2*pi*my)/dy;
            %kzm = (-1j)*sqrt(-(k0.^2 - kxm.^2 - kym.^2));
            
            krho2 = sqrt(kxm.^2 + kym.^2);
            k_rho_2 =  krho2.^2;
            
            
            
            kz0 = (-1j)*sqrt(-(k0.^2 - krho2.^2));
    
            z = h+eps;
            
            zeta_s = zeta./sqrt(er); 
            
            
            ksub = k0 .* sqrt(er); %why the hell man? change it to er :P
            
            kzs = (-1j)*sqrt(-(ksub.^2 - krho2.^2));

% For TM
z0_TM = zeta * kz0 ./ k0;
zs_TM = zeta_s * kzs ./ ksub;

% For TE

z0_TE = zeta * k0 ./ kz0;
zs_TE = zeta_s * ksub ./ kzs;

[zup_TE, zdn_TE, zup_TM, zdn_TM] = zupdn(z0_TE, zs_TE, z0_TM, zs_TM, h, kzs);


[vtm, vte, itm, ite] = txline(zup_TE, zdn_TE, zup_TM, zdn_TM, h, z, kz0, kzs, z0_TE, zs_TE, z0_TM, zs_TM);


[Gxx_1, Gyx_1, Gzx_1] = Green(vtm, vte, itm, ite, kxm, kym, k_rho_2, zeta, k0);


            const = 1/dy;
            
            D_int = const .* (Gxx_1) .* besselj(0, (kym .* w)/2);
            
            D_inf = sum(D_int, 1);
            
    
            
            const2 = (-1/dx);
            
            zin_int = const2 .* sinc((kxm(1,:) .* deld)/2/pi).^2./(D_inf);
            
            yin(p) = sum(zin_int);
            zin(p) = 1./yin(p);
            
            z0 = 400;
            
            Gamma(p) = (zin(p) - z0) / (zin(p) + z0);
            
            
            
end

% figure(1);
% 
% plot(f*10^(-9), real(zin), 'LineWidth', 3);
% grid on;
% hold on;
% plot(f*10^(-9), imag(zin), 'LineWidth', 3);
% hold on;

figure(1);

    plot(f*10^(-9), 20*log10(abs(Gamma)), 'LineWidth', 3);
    grid on;
    hold on;



end

% xlabel('Frequency(GHz)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Re Z_{in} , Im Z_{in} (\Omega) ', 'FontSize', 12, 'FontWeight', 'bold');
% title('Active input impedance', 'FontSize', 12, 'FontWeight', 'bold');

xlabel('Frequency(GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Active \Gamma (dB)', 'FontSize', 12, 'FontWeight', 'bold');
title(['Active reflection coefficient at z0 = ', num2str(z0), '\Omega'], 'FontSize', 12, 'FontWeight', 'bold');


%legend({'Real Z_{in}(0, 0)','Imag Z_{in}(0, 0)', 'Real Z_{in}(45, 0) E plane','Imag Z_{in}(45, 0) E plane', 'Real Z_{in}(45, 90) H plane','Imag Z_{in}(45, 90) H Plane'},'Location','northeast', 'FontSize', 12, 'FontWeight', 'bold');
%ylim([-2000 3000]);
%print('A5_Q1_dipole_br', '-dpng');
%xlim([-180 180]);

legend({'\Gamma (0, 0)', '\Gamma (45, 0) E Plane', '\Gamma (45, 90) H Plane'},'Location','northeast', 'FontSize', 12, 'FontWeight', 'bold');


