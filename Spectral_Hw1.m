%%
f = [10e9, 20e9];
lambda = 3e8./f;
zeta = 120*pi;
epsilon_r = 10;
h = 2e-3;

for i = 1:size(f, 2)

k0 = 2 * pi / lambda(i);

kx = linspace(0, 5*k0, 1000);
ky = 0;
kz0 = -1j * sqrt(-(k0.^2 - kx.^2 - ky.^2));
k_rho_2 = kx.^2 + ky.^2

zeta_s = zeta/sqrt(epsilon_r);



ksub = k0 * sqrt(epsilon_r);
kzs = -1j * sqrt(-(ksub.^2 - kx.^2 - ky.^2));



z = 2.2e-3;

% For TM
z0_TM = zeta * kz0 ./ k0;
zs_TM = zeta_s * kzs ./ ksub;

% For TE

z0_TE = zeta * k0 ./ kz0;
zs_TE = zeta_s * ksub ./ kzs;

[zup_TE, zdn_TE, zup_TM, zdn_TM] = zupdn(z0_TE, zs_TE, z0_TM, zs_TM, h, kzs);


[vtm, vte, itm, ite] = txline(zup_TE, zdn_TE, zup_TM, zdn_TM, h, z, kz0, kzs, z0_TE, zs_TE, z0_TM, zs_TM);


[G_xx, G_yx, G_zx] = Green(vtm, vte, itm, ite, kx, ky, k_rho_2, zeta, k0);





plot(kx/k0, abs(G_xx), 'LineWidth', 2);
ylim([0 10^5]);
hold on;

end

xlabel('k_x/k_0', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('|G^{ej}_{xx}(V)|', 'FontSize', 12, 'FontWeight', 'bold');
title('Spectral Green Function in x Direction', 'FontSize', 12, 'FontWeight', 'bold');
legend({'10GHz', '20GHz'},'Location','north', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

print('Q1_Spec_Gej', '-depsc');



%%

%%








%[SGFxx, SGFxy, SGFxz, SGFyx, SGFyy, SGFyz, SGFzx, SGFzy, SGFzz] = SGF(Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz, c);

%% Fourier Transform of Current distribution

%depends on the current distribution:


