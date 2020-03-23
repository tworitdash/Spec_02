%%
close all;
clear all;

f = 30e9;
lambda = 3e8./f;
zeta = 120*pi;
epsilon_r = 12;
h = 5e-3;



for i = 1:size(epsilon_r, 2)
k0 = 2 * pi / lambda;

ky = k0;
kx = linspace(eps, 2.*k0, 1000);
krho = sqrt(kx.^2 + ky.^2);

kz0 = -1j * sqrt(-(k0.^2 - kx.^2 - ky.^2));

zeta_s = zeta/sqrt(epsilon_r(i));



ksub = k0 * sqrt(epsilon_r(i));
kzs = -1j * sqrt(-(ksub.^2 - kx.^2 - ky.^2));



z = h+eps;

% For TM
z0_TM = zeta * kz0 ./ k0;
zs_TM = zeta_s * kzs ./ ksub;

% For TE

z0_TE = zeta * k0 ./ kz0;
zs_TE = zeta_s * ksub ./ kzs;

[zin_TE, zin_TM] = Zin(z0_TE, zs_TE, z0_TM, zs_TM, h, kz0);

[Gamma_TE, Gamma_TM] = Gamma_Q3(zs_TE, zs_TM, z0_TE, z0_TM, h, kz0);

[vte, vtm, ite, itm] = txline_3(Gamma_TE, Gamma_TM, h, kzs, kz0, zin_TE, zin_TM, zs_TE, zs_TM, z);


%Spectral Green's function

Gem_yx = (vte .* kx.^2 + vtm .* ky.^2) ./ (krho.^2);

plot(kx/k0, abs(Gem_yx), 'LineWidth', 2);
hold on;
%ylim([0 300]);
end

xlabel('k_x/k_0', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('|G^{em}_{yx}(V)|', 'FontSize', 12, 'FontWeight', 'bold');
title('Y component of Spectral Greens Function for a magnetic source at z = 0', 'FontSize', 12, 'FontWeight', 'bold');
legend({'\epsilon_r = 2.5', '\epsilon_r = 4', '\epsilon_r = 12'},'Location','north', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

print('Q3_Gem', '-depsc')