%%
f = 9e9:0.5e9:11e9;
lambda = 3e8./f;
zeta = 120*pi;
epsilon_r = 12;
h = 15e-3;
hs = 2.1e-3;


for i = 1:size(f, 2)
k0 = 2 * pi / lambda(i);

kx = 0;
ky = linspace(0, k0, 1000);
krho = sqrt(kx.^2 + ky.^2);

kz0 = -1j * sqrt(-(k0.^2 - kx.^2 - ky.^2));

zeta_s = zeta/sqrt(epsilon_r);



ksub = k0 * sqrt(epsilon_r);
kzs = -1j * sqrt(-(ksub.^2 - kx.^2 - ky.^2));



z = 19.1e-3;

% For TM
z0_TM = zeta * kz0 ./ k0;
zs_TM = zeta_s * kzs ./ ksub;

% For TE

z0_TE = zeta * k0 ./ kz0;
zs_TE = zeta_s * ksub ./ kzs;

[zin_TE, zin_TM, zL_TE, zL_TM, z02_TE, z02_TM] = Zall(z0_TE, zs_TE, z0_TM, zs_TM, h, hs, kzs, kz0);

[Gamma1_TE, Gamma1_TM, Gamma2_TE, Gamma2_TM] = Gamma(zs_TE, zs_TM, zL_TE, zL_TM, z0_TE, z0_TM, h, hs, kz0, kzs);

[vte, vtm, ite, itm] = txline_2(Gamma1_TE, Gamma1_TM, Gamma2_TE, Gamma2_TM, h, hs, kzs, kz0, zin_TE, zin_TM, zL_TE, zL_TM, z02_TE, z02_TM, z);


%Spectral Green's function

Gem_yx = (vte .* kx.^2 + vtm .* ky.^2) ./ (krho.^2);

plot(ky/k0, abs(Gem_yx), 'LineWidth', 2);
hold on;
%ylim([0 20]);
end

xlabel('k_y/k_0', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('|G^{em}_{yx}(V)|', 'FontSize', 12, 'FontWeight', 'bold');
title('Y component of Spectral Greens Function for a magnetic source at z = 0', 'FontSize', 12, 'FontWeight', 'bold');
legend({'9GHz', '9.5GHz', '10GHz', '10.5GHz', '11GHz'},'Location','northeast', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

print('Q2_Gem', '-depsc')