%%
close all;
clear all;

epsilon_r = linspace(4, 25, 25);
%epsilon_r = [1.5];
%epsilon_r = 25;
f = linspace(3e9, 15e9, 1000);
Dir = zeros(size(f));
BW = zeros(size(epsilon_r));
LegendInfo = zeros(size(epsilon_r));

for i = 1:size(epsilon_r, 2)
    
for k = 1:size(f, 2)


c = 3e8;
lambda = c./f(k);
zeta = 120*pi;

h = 15e-3;


l = 15e-3;
w = 1.5e-3;

r_obs = 1000 .* lambda;

k0 = 2 * pi / lambda;

hs = lambda./(4 .* sqrt(epsilon_r(i)));

zeta_s = zeta/sqrt(epsilon_r(i));

ksub = k0 * sqrt(epsilon_r(i));


%% Directivity
dth = pi/180;
dph = dth;

[th, ph] = meshgrid(eps:dth:pi/2, eps:dph:2*pi);

kx = k0 .* sin(th) .* cos(ph);
ky = k0 .* sin(th) .* sin(ph);
kz = k0 .* cos(th);

krho = sqrt(kx.^2 + ky.^2);

kz0 = -1j * sqrt(-(k0.^2 - krho.^2));

kzs = -1j * sqrt(-(ksub.^2 - krho.^2));



z = h + hs + eps;

% For TM
z0_TM = zeta * kz0 ./ k0;
zs_TM = zeta_s * kzs ./ ksub;

% For TE

z0_TE = zeta * k0 ./ kz0;
zs_TE = zeta_s * ksub ./ kzs;

[zin_TE, zin_TM, zL_TE, zL_TM, z02_TE, z02_TM] = Zall(z0_TE, zs_TE, z0_TM, zs_TM, h, hs, kzs, kz0);

[Gamma1_TE, Gamma1_TM, Gamma2_TE, Gamma2_TM] = Gamma(zs_TE, zs_TM, zL_TE, zL_TM, z0_TE, ...
    z0_TM, h, hs, kz0, kzs);

[vte, vtm, ite, itm] = txline_2(Gamma1_TE, Gamma1_TM, Gamma2_TE, Gamma2_TM, h, hs,...
    kzs, kz0, zin_TE, zin_TM, zL_TE, zL_TM, z02_TE, z02_TM, z);


%Spectral Green's function

[G_xx, G_yx, G_zx] = Green_em(vtm, vte, itm, ite, kx, ky, krho.^2, zeta, k0);

[Mx, My, Mz] = JFT_freq(l, w, k0, kx, ky);


E_far_x = 1j .* kz .* Mx .* G_xx .* exp(-1j .* kz .* abs(z)) .* exp(-1j .* k0 .* r_obs) ./ (2 * pi * r_obs);
E_far_y = 1j .* kz .* Mx .* G_yx .* exp(-1j .* kz .* abs(z)) .* exp(-1j .* k0 .* r_obs) ./ (2 * pi * r_obs);
E_far_z = 1j .* kz .* Mx .* G_zx .* exp(-1j .* kz .* abs(z)) .* exp(-1j .* k0 .* r_obs) ./ (2 * pi * r_obs);


E_abs = sqrt(abs(E_far_x).^2 + abs(E_far_y).^2 + abs(E_far_z).^2);

E_th = E_far_x .* cos(th) .* cos(ph) + E_far_y .* cos(th) .* sin(ph) - E_far_z .* sin(th);
E_ph = -sin(ph) .* E_far_x + cos(ph) .* E_far_y; 



V_far_th = E_th * r_obs./(exp(-1j * k0 * r_obs));
V_far_ph = E_ph * r_obs./(exp(-1j * k0 * r_obs));
    
C_rad = 1/(2 * zeta);
V_far_abs = sqrt(abs(V_far_th).^2 + abs(V_far_ph).^2);
U = C_rad * (V_far_abs).^2;

P_rad_int = C_rad * (V_far_abs).^2 .* sin(th) * dth * dph;

P_rad = nansum(nansum(P_rad_int));
    
Dir_i = 4 * pi * U ./ P_rad;

Dir(k) = Dir_i(1, 1);

legendInfo{i} = ['\epsilon_r = ',  num2str(epsilon_r(i))];


end

figure(3);
hold on;

plot(f*10^(-9), 10*log10(abs(Dir)), 'LineWidth', 2);

grid on;
xlabel('frequency(GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Directivity(dB Scale)', 'FontSize', 12, 'FontWeight', 'bold');
title(['D(\theta = 0, \phi = 0)'], 'FontSize', 12, 'FontWeight', 'bold');
legend(legendInfo, 'FontSize', 12, 'FontWeight', 'bold');




%% Bandwidth 

ind = interp1(Dir,1:length(Dir),max(Dir)./2,'nearest');
fl = f(ind);

D_temp = Dir;

D_temp(ind:ind+1) = eps.*[1, 2];

ind2 = interp1(D_temp,1:length(D_temp),max(D_temp)./2,'nearest');
fh = f(ind2);

BW(i) = 200 .* abs(fh - fl)./(fh + fl);



end

%print(['A4_Q1_4'], '-depsc');

% figure(1);
% 
% plot(epsilon_r, BW, 'LineWidth', 2, 'Color', [0.6350 0.0780 0.1840]);
% 
% grid on;
% xlabel('\epsilon_r', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('BW(%)', 'FontSize', 12, 'FontWeight', 'bold');
% title(['relative BW vs permittivity \epsilon_r'], 'FontSize', 12, 'FontWeight', 'bold');

