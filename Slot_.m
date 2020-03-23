%% Antenna dimesnions
clear all;
close all;

c0 = 3e8;
dth = pi/180;
dph = pi/180;
er = 1;
f = 2e9:0.05e9:10e9;

lambda0 = c0./10^10;

w = 0.2 .* lambda0;
deld = 0.25 .* lambda0;
dx = 0.5 .* lambda0;
dy = 0.5 .* lambda0;
hd = 0.31 .* lambda0;

% w = 0.1 .* lambda0;
% deld = 0.1 .* lambda0;
% dx = 0.5 .* lambda0;
% dy = 0.5 .* lambda0;
% hd = 0.25 .* lambda0;


th = pi/4;
ph = pi/2; % For E plane; change it to 90 for H plane

for j = 1:size(th, 2)

zin = zeros(size(f, 1), size(f, 2));
Gamma = zeros(size(f, 1), size(f, 2));
    
for p=1:size(f, 2)
    
    [mx, my] = meshgrid(-20:1:20, -20:1:20);
    
        lambda = c0./f(p);
    
        k0 = 2 * pi ./ lambda;
    
        zeta = 120*pi;
   
        
  

            kx0 = k0 .* sin(th(j)) .* cos(ph);
            ky0 = k0 .* sin(th(j)) .* sin(ph);
            
            %kz0 = k0 .* cos(th(j));

            kxm = kx0 - (2*pi*mx)/dx;
            kym = ky0 - (2*pi*my)/dy;
            
            %krho = sqrt(kx0.^2 + kym.^2);
            
            krho2 = sqrt(kxm.^2 + kym.^2);
            
            kz0 = (-1j)*sqrt(-(k0.^2 - krho2.^2));
            
            
    
            z0_TM = zeta .* kz0./k0;
            z0_TE = zeta .* k0./kz0;
           
          
           [vte_PPW, vtm_PPW, ite_PPW, itm_PPW] = txline_Slot(z0_TE, z0_TM, er, k0, krho2, hd, 'PPW');

           [vte_PEC, vtm_PEC, ite_PEC, itm_PEC] = txline_Slot(z0_TE, z0_TM, er, k0, krho2, hd, 'PEC');
           
           [Gxx_PPW, Gyx_PPW, Gzx_PPW] = Green_hm(vtm_PPW, vte_PPW, itm_PPW, ite_PPW, kxm, kym, krho2, zeta, k0);
           
           [Gxx_PEC, Gyx_PEC, Gzx_PEC] = Green_hm(vtm_PEC, vte_PEC, itm_PEC, ite_PEC, kxm, kym, krho2, zeta, k0);
           
           const = 1/dy;
            
            D_int = const .* (Gxx_PPW + Gxx_PEC) .* besselj(0, (kym .* w)/2);
            
            D_inf = sum(D_int, 1);
            const2 = (-1/dx);
            zin_int = const2 .* sinc((kxm(1,:) .* deld)/2/pi).^2 ./(D_inf);
            
            zin(p) = sum(zin_int);
            
            z0 = 400;
            
            Gamma(p) = (zin(p) - z0) / (zin(p) + z0);
            
            
            
end
%     figure(1);
% 
%     plot(f/f(end), real(zin), 'LineWidth', 3);
%     grid on;
%     hold on;
%     plot(f/f(end), imag(zin), 'LineWidth', 3);
%     hold on;   
% 
    figure(1);

    plot(f*10^(-9), 20*log10(abs(Gamma)), 'LineWidth', 3);
    grid on;
    hold on;




end


% xlabel('Frequency(GHz)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Re Z_{in} , Im Z_{in} (\Omega) ', 'FontSize', 12, 'FontWeight', 'bold');
% title('Active input impedance', 'FontSize', 12, 'FontWeight', 'bold');
% legend({'Real Z_{in}(0, 0)','Imag Z_{in}(0, 0)', 'Real Z_{in}(45, 0)','Imag Z_{in}(45, 0)'},'Location','northeast', 'FontSize', 12, 'FontWeight', 'bold');

% xlabel('f/f_0', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Re Z_{in} , Im Z_{in} (\Omega) ', 'FontSize', 12, 'FontWeight', 'bold');
% title('Active input impedance at \theta_0 = 0', 'FontSize', 12, 'FontWeight', 'bold');
% legend({'Real Z_{in}','Imag Z_{in}'},'Location','northeast', 'FontSize', 12, 'FontWeight', 'bold');
% print('A5_Q1_slot_br', '-depsc');

xlabel('Frequency(GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Active \Gamma (dB)', 'FontSize', 12, 'FontWeight', 'bold');
title(['Active reflection coefficient at z0 = ', num2str(z0), '\Omega'], 'FontSize', 12, 'FontWeight', 'bold');

% legend({'\Gamma (0, 0)', '\Gamma (45, 0) H Plane', '\Gamma (45, 90) E Plane'},'Location','northeast', 'FontSize', 12, 'FontWeight', 'bold');

%ylim([-2000 3000]);
%print('Admittance_in', '-dpng');
%xlim([-180 180]);



