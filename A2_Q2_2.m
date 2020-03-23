%%
clear all;
close all;

f = 10e9;
lambda = 3e8./f;
zeta = 120*pi;
epsilon_r = linspace(1, 25, 1000);
h = 15e-3;

Dir_i = zeros(size(epsilon_r, 1), size(epsilon_r, 2));





for i = 1:size(epsilon_r, 2)
k0 = 2 * pi / lambda;
l = lambda./2;
w = lambda./20;
r_obs = 1000 .* lambda;
dth = pi/180;
dph = pi/180;

hs = lambda./(4 .* sqrt(epsilon_r(i)));

[th, ph] = meshgrid(eps:dth:pi/2, eps:dph:2*pi);

kx = k0 .* sin(th) .* cos(ph);
ky = k0 .* sin(th) .* sin(ph);
kz = k0 .* cos(th);

k_rho_2 = kx.^2 + ky.^2;

kz0 = -1j * sqrt(-(k0.^2 - kx.^2 - ky.^2));

zeta_s = zeta/sqrt(epsilon_r(i));



ksub = k0 * sqrt(epsilon_r(i));
keq = (k0 + ksub)./2;

kzs = -1j * sqrt(-(ksub.^2 - kx.^2 - ky.^2));



z = h + hs + eps;

% For TM
z0_TM = zeta * kz0 ./ k0;
zs_TM = zeta_s * kzs ./ ksub;

% For TE

z0_TE = zeta * k0 ./ kz0;
zs_TE = zeta_s * ksub ./ kzs;

[zin_TE, zin_TM, zL_TE, zL_TM, z02_TE, z02_TM] = Zall(z0_TE, zs_TE, z0_TM, zs_TM, h, hs, kzs, kz);

[Gamma1_TE, Gamma1_TM, Gamma2_TE, Gamma2_TM] = Gamma(zs_TE, zs_TM, zL_TE, zL_TM, z0_TE, z0_TM, h, hs, kz, kzs);

[vte, vtm, ite, itm] = txline_2(Gamma1_TE, Gamma1_TM, Gamma2_TE, Gamma2_TM, h, hs, kzs, kz, zin_TE, zin_TM, zL_TE, zL_TM, z02_TE, z02_TM, z);


%Spectral Green's function

[G_xx, G_yx, G_zx] = Green_em(vtm, vte, itm, ite, kx, ky, k_rho_2, zeta, k0);

[Mx, My, Mz] = JFT_freq(l, w, k0, kx, ky);


E_far_x = 1j .* kz .* Mx .* G_xx .* exp(-1j .* kz .* abs(z - h)) .* exp(-1j .* k0 .* r_obs) ./ (2 * pi * r_obs);
E_far_y = 1j .* kz .* Mx .* G_yx .* exp(-1j .* kz .* abs(z - h)) .* exp(-1j .* k0 .* r_obs) ./ (2 * pi * r_obs);
E_far_z = 1j .* kz .* Mx .* G_zx .* exp(-1j .* kz .* abs(z - h)) .* exp(-1j .* k0 .* r_obs) ./ (2 * pi * r_obs);


E_abs = sqrt(abs(E_far_x).^2 + abs(E_far_y).^2 + abs(E_far_z).^2);

E_th = E_far_x .* cos(th) .* cos(ph) + E_far_y .* cos(th) .* sin(ph) - E_far_z .* sin(th);
E_ph = -sin(ph) .* E_far_x + cos(ph) .* E_far_y; 

% figure(i);
% 
% plot(th(1, :)*180/pi, db(abs(E_th(1, :))/max(abs(E_th(1, :)))), 'LineWidth', 2);
% grid on;
% hold on;
% 
% plot(th(91, :)*180/pi, db(abs(E_ph(91, :))/max(abs(E_ph(91, :)))), 'LineWidth', 2);
% 
% ylim([-50 0]);
% 
% xlabel('\theta(deg)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('E(\theta)_{\phi} = 0), E(\phi)_{\phi = 0})', 'FontSize', 12, 'FontWeight', 'bold');
% title('Normalized Far Electric Field at different cuts', 'FontSize', 12, 'FontWeight', 'bold');
% legend({['E_{\theta}(\phi = 0)', num2str(f(i)*10^(-9)), 'GHz'], ['E_{\phi}(\phi = 90)', num2str(f(i)*10^(-9)), 'GHz']},'Location','northeast', 'FontSize', 12, 'FontWeight', 'bold');


V_far_th = E_th * r_obs./(exp(-1j * k0 * r_obs));
V_far_ph = E_ph * r_obs./(exp(-1j * k0 * r_obs));
    
C_rad = 1/(2 * zeta);
V_far_abs = sqrt(abs(V_far_th).^2 + abs(V_far_ph).^2);
U = C_rad * (V_far_abs).^2;

P_rad_int = C_rad * (V_far_abs).^2 .* sin(th) * dth * dph;

P_rad = nansum(nansum(P_rad_int));
    
Dir_i(i) = 4 * pi * U(1, 1) ./ P_rad;


end
figure(1);

plot(epsilon_r, abs(Dir_i), 'LineWidth', 2, 'Color', [0.6350, 0.0780, 0.1840]);

grid on;
xlabel('\epsilon_r', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Directivity(Linear Scale) at \theta = 0', 'FontSize', 12, 'FontWeight', 'bold');
title('Directivity', 'FontSize', 12, 'FontWeight', 'bold');
xlim([1 25]);



%legend({['D(\phi = 0)', num2str(f(i)*10^(-9)), 'GHz'], ['D(\phi = 90)', num2str(f(i)*10^(-9)), 'GHz']},'Location','northeast', 'FontSize', 12, 'FontWeight', 'bold');

print(['D_Q2_2'], '-depsc');

figure(2);
plot(epsilon_r, 10*log10(abs(Dir_i)), 'LineWidth', 2, 'Color', [0.6350, 0.0780, 0.1840]);

grid on;
xlabel('\epsilon_r', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Directivity(dB) at \theta = 0', 'FontSize', 12, 'FontWeight', 'bold');
title('Directivity', 'FontSize', 12, 'FontWeight', 'bold');
xlim([1 25]);
print(['D_Q2_2_dB'], '-depsc');