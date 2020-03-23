%% 


f = 20e9; % max frequency to do the guess point analysis

lambda = 3e8./f;
zeta = 120*pi;
epsilon_r = 10;
h = 2e-3;


r_obs = 1000 .* lambda;

[th, ph] = meshgrid(-pi/2-eps:pi/180:pi/2-eps, eps:pi/180:2*pi);

k0 = 2 * pi ./ lambda;
ksub = k0 * sqrt(epsilon_r);

% kx = linspace(k0, ksub, 100);
% ky = 0;


%k_rho_2 = kx.^2 + ky.^2;

krho_i = linspace(k0, ksub, 100);

kz0 = -1j * sqrt(-(k0.^2 - krho_i.^2));

zeta_s = zeta/sqrt(epsilon_r);

keq = k0 .* sqrt((1 + epsilon_r)./2);

kzs = -1j * sqrt(-(ksub.^2 - krho_i.^2));


z = h/2;

% For TM
z0_TM = zeta * kz0 ./ k0;
zs_TM = zeta_s * kzs ./ ksub;

% For TE

z0_TE = zeta * k0 ./ kz0;
zs_TE = zeta_s * ksub ./ kzs;

[zup_TE, zdn_TE, zup_TM, zdn_TM] = zupdn(z0_TE, zs_TE, z0_TM, zs_TM, h, kzs);

[D_krho_TE] = D_kro(k0, "TE", epsilon_r, h, zeta, zeta_s, krho_i);
[D_krho_TM] = D_kro(k0, "TM", epsilon_r, h, zeta, zeta_s, krho_i);

% D_krho_TE = zup_TE + zdn_TE;
% D_krho_TM = zup_TM + zdn_TM;

% figure(1);
% 
% plot(krho_i/k0, abs(D_krho_TE), 'LineWidth', 2);
% hold on;
% plot(krho_i/k0, abs(D_krho_TM), 'LineWidth', 2);
% grid on;
% 
% xlabel('k_{\rho}/k_0', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Denominator D(k_{\rho})', 'FontSize', 12, 'FontWeight', 'bold');
% title(['Denominator function for finding poles at \epsilon_r = ', num2str(epsilon_r)], 'FontSize', 12, 'FontWeight', 'bold');
% legend({'D(TE)', 'D(TM)'}, 'Location', 'north', 'FontSize', 12, 'FontWeight', 'bold');

%print(['A3_Q1_K_guess_' , num2str(epsilon_r)], '-depsc');

krho_g_TE = 1.874 .* k0;
krho_g_TM = 2.616 .* k0;

[krho_TE, f_axis_1] = finddrop(k0, epsilon_r, krho_g_TE, h, zeta, zeta_s, "TE");
[krho_TM, f_axis_2] = finddrop(k0, epsilon_r, krho_g_TM, h, zeta, zeta_s, "TM");



% figure(2);
% 
% plot(f_axis_1 * 10^(-9), abs(krho_TE), 'LineWidth', 2);
% hold on;
% plot(f_axis_2 * 10^(-9), abs(krho_TM), 'LineWidth', 2);
% grid on;
% 
% ylim([1 3]);

% xlabel('frequency(GHz)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('k_{\rho g}(f_i)', 'FontSize', 12, 'FontWeight', 'bold');
% title(['Guess pole point vs Frequency  at \epsilon_r = ', num2str(epsilon_r)], 'FontSize', 12, 'FontWeight', 'bold');
% legend({'k_{\rho g}(f_i)(TE)', 'k_{\rho g}(f_i)(TM)'}, 'Location', 'north', 'FontSize', 12, 'FontWeight', 'bold');

%print(['A3_Q1_K_guess_2_', num2str(epsilon_r)], '-depsc');

f0 = flip(f_axis_1);

Psw = zeros(size(f0, 1), size(f0, 2));
Psw_PWS = zeros(size(f0, 1), size(f0, 2));
Psw_Uni = zeros(size(f0, 1), size(f0, 2));

for i = 1:size(f0, 2)
    
lambda0 = 3e8./f0(i);

k0_ = (2*pi*f0(i))./(3e8);

k_sw = krho_TM(f_axis_1 == f0(i)) .* k0_;

z = linspace(eps, 1, 1000);
phi = eps:pi/180:2*pi;

VtmR = zeros(size(z, 1), size(z, 2));
ItmR = zeros(size(z, 1), size(z, 2));

for j = 1:size(z, 2)

    [VtmR(j), ItmR(j)] = Residue_GroundSlab(k0_, epsilon_r, h, k_sw, z(j));

end

rho = 1;


[Psw(i)] = PswTMelem(k0_, epsilon_r, h, k_sw, ItmR, z);

l = 5.3e-3;
w = 0.5e-3;

I_phi_del = pi;
   

Iphi_pws = PhiInt_PWS(k0_, k_sw, l, w, epsilon_r);

Iphi_uni = PhiInt_Uniform(k0_, k_sw, l, w);

Psw_PWS(i) = Psw(i) .* Iphi_pws ./ I_phi_del;

Psw_Uni(i) = Psw(i) .* Iphi_uni ./ I_phi_del;



end
%% Q1.2 


f = f0;
P_rad_i = zeros(size(f0, 1), size(f0, 2));


P_rad_i_d = zeros(size(f0, 1), size(f0, 2));
eta_sw_elem = zeros(size(f0, 1), size(f0, 2));
eta_sw_PWS = zeros(size(f0, 1), size(f0, 2));
eta_sw_Uni = zeros(size(f0, 1), size(f0, 2));

for i = 1:size(f0, 2)
    
lambda = 3e8./f0(i);

zeta = 120*pi;
epsilon_r = 10;
h = 2e-3;
l = 5.3e-3;
w = 0.5e-3;

r_obs = 1000 .* lambda;

dth = pi/180;
dph = pi/180;

[th, ph] = meshgrid(eps:dth:pi/2-eps, eps:dph:2*pi);

k0 = 2 * pi / lambda;

kx = k0 .* sin(th) .* cos(ph);
ky = k0 .* sin(th) .* sin(ph);
kz = k0 .* cos(th);

kz0 = -1j * sqrt(-(k0.^2 - kx.^2 - ky.^2));
k_rho_2 = kx.^2 + ky.^2;

zeta_s = zeta/sqrt(epsilon_r);

ksub = k0 * sqrt(epsilon_r);

keq = (k0 + ksub)./2;

kzs = -1j * sqrt(-(ksub.^2 - kx.^2 - ky.^2));
%kzs = kz * sqrt(epsilon_r);


z = h+eps;

% For TM
z0_TM = zeta * kz ./ k0;
zs_TM = zeta_s * kzs ./ ksub;

% For TE

z0_TE = zeta * k0 ./ kz;
zs_TE = zeta_s * ksub ./ kzs;

[zup_TE, zdn_TE, zup_TM, zdn_TM] = zupdn(z0_TE, zs_TE, z0_TM, zs_TM, h, kzs);


[vtm, vte, itm, ite] = txline(zup_TE, zdn_TE, zup_TM, zdn_TM, h, z, kz, kzs, z0_TE, zs_TE, z0_TM, zs_TM);


[G_xx, G_yx, G_zx] = Green(vtm, vte, itm, ite, kx, ky, k_rho_2, zeta, k0);

%[Jx, Jy, Jz] = JFT_freq(l, w, keq, kx, ky); %pws

Jx = 1; % elementary
%[Jx, Jy, Jz] = JFT_uni(l, w, kx, ky); % uniform



E_far_x = 1j .* kz .* Jx .* G_xx .* exp(-1j .* kz .* abs(z - h)) .* exp(-1j .* k0 .* r_obs) ./ (2 * pi * r_obs);
E_far_y = 1j .* kz .* Jx .* G_yx .* exp(-1j .* kz .* abs(z - h)) .* exp(-1j .* k0 .* r_obs) ./ (2 * pi * r_obs);
E_far_z = 1j .* kz .* Jx .* G_zx .* exp(-1j .* kz .* abs(z - h)) .* exp(-1j .* k0 .* r_obs) ./ (2 * pi * r_obs);


E_abs = sqrt(abs(E_far_x).^2 + abs(E_far_y).^2 + abs(E_far_z).^2);

E_th = E_far_x .* cos(th) .* cos(ph) + E_far_y .* cos(th) .* sin(ph) - E_far_z .* sin(th);
E_ph = -sin(ph) .* E_far_x + cos(ph) .* E_far_y; 


% Power radiated:

%[Dir, P_rad] = Directivity(r_obs, th, ph, E_th, E_ph, zeta, dth, dph, k0);

V_far_th = E_th * r_obs./(exp(-1j * k0 * r_obs));
V_far_ph = E_ph * r_obs./(exp(-1j * k0 * r_obs));
    
C_rad = 1/(2 * zeta);
V_far_abs = sqrt(abs(V_far_th).^2 + abs(V_far_ph).^2);
U = C_rad * (V_far_abs).^2;

P_rad_int = C_rad * (V_far_abs).^2 .* sin(th) * dth * dph;

P_rad = sum(sum(P_rad_int));
    
Dir = 4 * pi * U ./ P_rad;

P_rad_i(i) = P_rad;

eta_sw_elem(i) = P_rad_i(i)./(P_rad_i(i) + Psw(i));

%eta_sw_PWS(i) = P_rad_i(i)./(P_rad_i(i) + Psw_PWS(i));
%eta_sw_Uni(i) = P_rad_i(i)./(P_rad_i(i) + Psw_Uni(i));

% For only one dipole

[th, ph] = meshgrid(eps:dth:pi-eps, eps:dph:2*pi);

kx = k0 .* sin(th) .* cos(ph);
ky = k0 .* sin(th) .* sin(ph);
kz = k0 .* cos(th);

c = (-zeta./(2 * k0 .* kz)); %constant term in the equations

%% Calculation of the Dyad
[Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz] = Dyad(k0, kx, ky, kz);

%% Calculation of Spectral Green's function

[SGFxx, SGFxy, SGFxz, SGFyx, SGFyy, SGFyz, SGFzx, SGFzy, SGFzz] = SGF(Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz, c);

[Jx_d, Jy_d, Jz_d] = JFT_freq(l, w, k0, kx, ky);
% Jx_d = 1;
% Jy_d = 0;
% Jz_d = 0;

c2 = 1j * kz * exp(-1j * k0 * r_obs) / (2 * pi * r_obs);

[Ex, Ey, Ez] = FF(c2, SGFxx, SGFxy, SGFxz, SGFyx, SGFyy, SGFyz, SGFzx, SGFzy, SGFzz, Jx_d, Jy_d, Jz_d);

E_abs_d = sqrt(abs(Ex.^2 + Ey.^2 + Ez.^2));

E_th_d = Ex .* cos(th) .* cos(ph) + Ey .* cos(th) .* sin(ph) - Ez .* sin(th);
E_ph_d = -sin(ph) .* Ex + cos(ph) .* Ey; 


% Power radiated:

%[Dir, P_rad] = Directivity(r_obs, th, ph, E_th, E_ph, zeta, dth, dph, k0);

V_th = E_th_d * r_obs./(exp(-1j * k0 * r_obs));
V_ph = E_ph_d * r_obs./(exp(-1j * k0 * r_obs));
    
C_rad = 1/(2 * zeta);
V_abs = sqrt(abs(V_th).^2 + abs(V_ph).^2);
U_d = C_rad * (V_abs).^2;

P_rad_int_d = C_rad * (V_abs).^2 .* sin(th) * dth * dph;

P_rad_d = nansum(nansum(P_rad_int_d));
    

P_rad_i_d(i) = P_rad_d;



end



% figure(4);
% plot(f * 10^(-9), 0.5 .* db(P_rad_i./P_rad_i_d), 'LineWidth', 2);
% hold on;
% plot(f * 10^(-9), 0.5 .* db(Psw./P_rad_i_d), 'LineWidth', 2);
% 
% xlabel('frequency(GHz)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('P_{rad} Normalized', 'FontSize', 12, 'FontWeight', 'bold');
% title('Normalized P_{rad} SW and Far field with free space dipole', 'FontSize', 12, 'FontWeight', 'bold');
% grid on;
% 
% print('A3_Q3_Prad_1', '-depsc');

%figure(5);
hold on;
plot(f * 10^(-9), eta_sw_elem, 'LineWidth', 2);
hold on;
%plot(f * 10^(-9), eta_sw_PWS, 'LineWidth', 2);
%plot(f * 10^(-9), eta_sw_Uni, 'LineWidth', 2);

xlabel('frequency(GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('\eta_{sw}', 'FontSize', 12, 'FontWeight', 'bold');
title('Surface wave efficiency', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

%print('A3_Q3_Eff', '-depsc');
%%



