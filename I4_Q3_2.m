%%
close all;
clear all;

f = 30e9;
c = 3e8;
lambda = c./f;
zeta = 120*pi;
epsilon_r = 12;
h = 5e-3;


% krho_TE_e = zeros(size(f));
% krho_TM_e = zeros(size(f));
% 
% krho_TE_a = zeros(size(f));
% krho_TM_a = zeros(size(f));
%  

l = lambda./2;
w = lambda./20;

k0 = 2 * pi / lambda;

% hs = lambda1./(4 .* sqrt(epsilon_r(m))); No hs for resonant lens feed antenna


%lambda = 3e8./f(k);
r_obs = 1000 .* lambda;
dth = pi/180;
dph = pi/180;


[th, ph] = meshgrid(-pi/2-eps:dth:pi/2-eps, eps:dph:2*pi);
    


kx = k0 .* sin(th) .* cos(ph);
ky = k0 .* sin(th) .* sin(ph);
kz = k0 .* cos(th);
k_rho_2 = (kx.^2 + ky.^2);

zeta_s = zeta/sqrt(epsilon_r);

ksub = k0 * sqrt(epsilon_r);
%kzs = -1j * sqrt(-(ksub.^2 - kx.^2 - ky.^2));

z = h+eps;
d = lambda./2;


kxsub = k0 .* sqrt(epsilon_r) .* sin(th) .* cos(ph);
kysub = k0 .* sqrt(epsilon_r) .* sin(th) .* sin(ph);
kzsub = k0 .* sqrt(epsilon_r) .* cos(th);

kzs = -1j * sqrt(-(ksub.^2 - kxsub.^2 - kysub.^2));

kz0 = -1j * sqrt(-(k0.^2 - kxsub.^2 - kysub.^2));

k_rho_2_s = (kxsub.^2 + kysub.^2);

% For TM
z0_TM = zeta * kz0 ./ k0;
zs_TM = zeta_s * kzs ./ ksub;

% For TE

z0_TE = zeta * k0 ./ kz0;
zs_TE = zeta_s * ksub ./ kzs;

%[zin_TE, zin_TM] = Zin(z0_TE, zs_TE, z0_TM, zs_TM, h, kz0);

[Gamma_TE, Gamma_TM] = Gamma_Q3(zs_TE, zs_TM, z0_TE, z0_TM, h, kz0);

[vte, vtm, ite, itm] = txline_3(Gamma_TE, Gamma_TM, h, kzs, kz0, z0_TE, z0_TM, zs_TE, zs_TM, z);


%Spectral Green's function

[G_xx, G_yx, G_zx] = Green_em(vtm, vte, itm, ite, kxsub, kysub, k_rho_2_s, zeta_s, ksub);

[Mx, My, Mz] = JFT_shift(l, w, k0, kxsub, kysub, d);


E_far_x = 1j .* kzsub .* Mx .* G_xx .* exp(-1j .* kzsub .* abs(z)) .* exp(-1j .* ksub .* r_obs) ./ (2 * pi * r_obs);
E_far_y = 1j .* kzsub .* Mx .* G_yx .* exp(-1j .* kzsub .* abs(z)) .* exp(-1j .* ksub .* r_obs) ./ (2 * pi * r_obs);
E_far_z = 1j .* kzsub .* Mx .* G_zx .* exp(-1j .* kzsub .* abs(z)) .* exp(-1j .* ksub .* r_obs) ./ (2 * pi * r_obs);


E_abs = sqrt(abs(E_far_x).^2 + abs(E_far_y).^2 + abs(E_far_z).^2);

E_th = E_far_x .* cos(th) .* cos(ph) + E_far_y .* cos(th) .* sin(ph) - E_far_z .* sin(th);
E_ph = -sin(ph) .* E_far_x + cos(ph) .* E_far_y; 


figure(1);

plot(th(1, :)*180/pi, db(abs(E_th(91, :))/max(max(abs(E_ph(1, :))), max(E_th(91, :)))), 'LineWidth', 2);
grid on;
hold on;

plot(th(91, :)*180/pi, db(abs(E_ph(1, :))/max(max(abs(E_ph(1, :))), max(E_th(91, :)))), 'LineWidth', 2);

ylim([-50 0]);

xlabel('\theta(deg)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('E(\theta)_{\phi} = 90), E(\phi)_{\phi = 0})', 'FontSize', 12, 'FontWeight', 'bold');
title('Normalized Far Electric Field at different cuts', 'FontSize', 12, 'FontWeight', 'bold');
legend({['E_{\theta}(\phi = 90)', num2str(f*10^(-9)), 'GHz'], ['E_{\phi}(\phi = 0)', num2str(f*10^(-9)), 'GHz']},'Location','northeast', 'FontSize', 12, 'FontWeight', 'bold');

print('A4_Q3_FF', '-depsc');
% % Directivity
% 
% V_far_th = E_th * r_obs./(exp(-1j * k0 * r_obs));
% V_far_ph = E_ph * r_obs./(exp(-1j * k0 * r_obs));
%     
% C_rad = 1/(2 * zeta_s);
% V_far_abs = sqrt(abs(V_far_th).^2 + abs(V_far_ph).^2);
% U = C_rad * (V_far_abs).^2;
% 
% P_rad_int = C_rad * (V_far_abs).^2 .* sin(th) * dth * dph;
% 
% P_rad = nansum(nansum(P_rad_int));
%     
% Dir(k) = 4 * pi * U(1, 1) ./ P_rad;
% 
% legendInfo{i} = ['\epsilon_r = ',  num2str(epsilon_r(i))];
