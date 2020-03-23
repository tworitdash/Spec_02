%%
close all;
clear all;



f = 11e9;
c = 3e8;
lambda = c./f;
zeta = 120*pi;
epsilon_r = 12;
h = 15e-3;
hs = 2.1e-3;

l = lambda./2;
w = lambda./20;

k0 = 2 * pi / lambda;


zeta_s = zeta/sqrt(epsilon_r);

ksub = k0 * sqrt(epsilon_r);
%krho = linspace(k0, ksub, 1000);

% krho_real = linspace(eps, k0, 1000);
% krho_imag = linspace(-k0, eps, 1000);

[krho_real, krho_imag] = meshgrid(linspace(eps, k0, 1000), linspace(k0, eps, 1000));

krho = krho_real - 1j .* krho_imag;

% krho = zeros(size(krho_real, 2), size(krho_imag, 2));
% 
% for i = 1:size(krho_real, 2)
%     for k = 1:size(krho_imag, 2)
%         krho(i, k) = krho_real(i) + 1j .* krho_imag(k);
%     end
% end
%krho_i = combvec(krho_real, 1j.*krho_imag);

[D_krho_TE] = D_kro_leak(k0, "TE", epsilon_r, h, hs, zeta, zeta_s, krho);
[D_krho_TM] = D_kro_leak(k0, "TM", epsilon_r, h, hs, zeta, zeta_s, krho);

D_TE_max = max(abs(1./D_krho_TE(:)));
D_TM_max = max(abs(1./D_krho_TM(:)));

[Ind_x_TE, Ind_y_TE] = find(abs(1./D_krho_TE) == D_TE_max);
[Ind_x_TM, Ind_y_TM] = find(abs(1./D_krho_TM) == D_TM_max);

krho_g_TE = krho(Ind_x_TE, Ind_y_TE);

krho_g_TM = krho(Ind_x_TM, Ind_y_TM);

figure(1);

contourf(krho_real/k0, -krho_imag/k0, abs(1./D_krho_TE));
%imagesc(krho_real(1, :)/k0, -krho_imag(1, :)/k0, abs(1./D_krho_TE));

figure(2);

contourf(krho_real/k0, -krho_imag/k0, abs(1./D_krho_TM));


%% finddrop

[krho_TE, f_axis_1] = finddrop_I4(k0, epsilon_r, krho_g_TE, h, hs, zeta, zeta_s, "TE");
[krho_TM, f_axis_2] = finddrop_I4(k0, epsilon_r, krho_g_TM, h, hs, zeta, zeta_s, "TM");


%imagesc(krho_real(1, :)/k0, -krho_imag(:, 1)/k0, abs(1./D_krho_TM));

% figure(4);
% 
% surf(krho_real(1, :)/k0, krho_imag(:, 1)/k0, abs(1./D_krho_TE));shading flat;
% 
% figure(5);
% 
% surf(krho_real(1, :)/k0, krho_imag(:, 1)/k0, abs(1./D_krho_TM));shading flat;



%% Analytic solutions of Dispertion equations:




f = flip(f_axis_1);
lambda = c./f;
hb = h./lambda;
k0 = (2 .* pi)./lambda;

kz0_TE = (k0./(pi.*epsilon_r.*(2.*hb).^2)) .* ((2.*pi.*hb.*epsilon_r + 1j)./(1 + (1./(2.*pi.*hb.*epsilon_r).^2)));
kz0_TM = (k0./(4 .* hb)).*(1 + sqrt(1 + 8j .* (hb./(pi.*epsilon_r))));

krho_TE_a = sqrt(k0.^2 - kz0_TE.^2);

krho_TM_a = sqrt(k0.^2 - kz0_TM.^2);


figure(3);

plot(f_axis_1 * 10^(-9), real(krho_TE), 'LineWidth', 2);
hold on;
plot(f*10^(-9), real(krho_TE_a./k0), '--', 'LineWidth', 2);
hold on;
plot(f_axis_2 * 10^(-9), real(krho_TM), 'LineWidth', 2);
hold on;
plot(f*10^(-9), real(krho_TM_a./k0), '--', 'LineWidth', 2);

grid on;

xlabel('frequency(GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('k_{\rho g}(f_i)/k_0', 'FontSize', 12, 'FontWeight', 'bold');
title(['Guess pole point vs Frequency  at \epsilon_r = ', num2str(epsilon_r)], 'FontSize', 12, 'FontWeight', 'bold');
legend({'Re(k_{\rho g}(f_i)(TE))', 'Re(k_{\rho g}(f_i)(TE)) Analytical', 'Re(k_{\rho g}(f_i)(TM))', 'Re(k_{\rho g}(f_i)(TM)) Analytical'}, 'Location', 'southeast', 'FontSize', 12, 'FontWeight', 'bold');
%print('A4_Q1_real', '-depsc');

figure(4);

plot(f_axis_1 * 10^(-9), imag(krho_TE), 'LineWidth', 2);
hold on;
plot(f*10^(-9), imag(krho_TE_a./k0), '--', 'LineWidth', 2);
hold on;
plot(f_axis_2 * 10^(-9), imag(krho_TM), 'LineWidth', 2);
hold on;
plot(f*10^(-9), imag(krho_TM_a./k0), '--', 'LineWidth', 2);
grid on;

xlabel('frequency(GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('k_{\rho g}(f_i)/k_0', 'FontSize', 12, 'FontWeight', 'bold');
title(['Guess pole point vs Frequency  at \epsilon_r = ', num2str(epsilon_r)], 'FontSize', 12, 'FontWeight', 'bold');
legend({'Im(k_{\rho g}(f_i)(TE))', 'Im(k_{\rho g}(f_i)(TE)) Analytical', 'Im(k_{\rho g}(f_i)(TM))', 'Im(k_{\rho g}(f_i)(TM)) Analytical'}, 'Location', 'southeast', 'FontSize', 12, 'FontWeight', 'bold');
%print('A4_Q1_imag', '-depsc');





%z = h + hs;

% [zin_TE, zin_TM, zL_TE, zL_TM, z02_TE, z02_TM] = Zall(z0_TE, zs_TE, z0_TM, zs_TM, h, hs, kzs, kz0);
% 
% [Gamma1_TE, Gamma1_TM, Gamma2_TE, Gamma2_TM] = Gamma(zs_TE, zs_TM, zL_TE, zL_TM, z0_TE, z0_TM, h, hs, kz0, kzs);
% 
% [vte, vtm, ite, itm] = txline_2(Gamma1_TE, Gamma1_TM, Gamma2_TE, Gamma2_TM, h, hs, kzs, kz0, zin_TE, zin_TM, zL_TE, zL_TM, z02_TE, z02_TM, z);