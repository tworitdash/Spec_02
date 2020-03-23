%%
f = 100e9;
lambda = 3e8/f;
k0 = 2 * pi / lambda;
eps_r = 1;
zeta = 120*pi;
r_obs = 10000 * lambda;
%% Circular dimensions

D = 6 * lambda;  %diameter of feed
f_number = 3;

d = 50 * lambda; %diameter of reflector
focal_length = f_number * d;

%% Incident electric field
E_0_pw = 1;

dth = pi/2000; dph = dth;
% For entire region of theta and phi
theta_0 = 2 * atan(1/(4 * f_number));

[th, ph] = meshgrid(eps:dth:theta_0-eps*dth, eps:dph:2*pi); 

%[th, ph] = meshgrid(eps:dth:pi-dth*eps, 0:dph:2*pi); 



e_go_th = - (2 ./ (1 + cos(th))) .* E_0_pw .* sin(ph);
e_go_ph = - (2 ./ (1 + cos(th))) .* E_0_pw .* cos(ph);

e_go_abs = sqrt(abs(e_go_th).^2 + abs(e_go_ph).^2);

%% Far field of feed 


[theta, phi] = meshgrid(eps:dth:pi/2-eps*dth, eps:dph:2*pi); 

%[theta, phi] = meshgrid(eps:dth:theta_0-eps*dth, eps:dph:2*pi); 


%Airy = pi * (D/2)^2 * besselj(1, (k0 * D/2 .* sin(theta))) ./ (k0 * (D/2) .* sin(theta));

C_feed = focal_length ./ exp(-1j * k0 * focal_length);

kx = k0 .* sin(theta) .* cos(phi);
ky = k0 .* sin(theta) .* sin(phi);
kz = k0 .* cos(theta);

c = (-zeta./(2 * k0 .* kz)); %constant term in the equations

%% Calculation of the Dyad
[Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz] = Dyad(k0, kx, ky, kz);

%% Calculation of Spectral Green's function

[SGFxx, SGFxy, SGFxz, SGFyx, SGFyy, SGFyz, SGFzx, SGFzy, SGFzz] = SGF(Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz, c);

[Jx, Jy, Jz] = J_Circular(D, k0 , theta);

%% Far field

c2 = 4j * exp(-1j * k0 * focal_length) / (2 * pi * focal_length).*kz;

[Ex, Ey, Ez] = FF(c2, SGFxx, SGFxy, SGFxz, SGFyx, SGFyy, SGFyz, SGFzx, SGFzy, SGFzz, Jx, Jy, Jz);

%E_abs = sqrt((abs(Ex).^2 + abs(Ey).^2 + abs(Ez).^2));

%% Current Distribution on the aperture:

c3 = -2j .* k0 .* zeta .* (exp(-1j * k0 * focal_length))./(2 .* pi * focal_length);
E_r = 0;
Eth = c3 .* Jy .* cos(theta) .* sin(phi);
Eph = c3 .* Jy .* cos(phi);

%Eth = Airy .* cos(theta) .* sin(phi)./C_feed;
%Eph = Airy .* cos(phi)./C_feed;

E_abs = sqrt(abs(Eth).^2 + abs(Eph).^2);

figure(1);
plot(th(91, :)*(180/pi), db(e_go_abs(91, :)./max(e_go_abs(91, :))), 'LineWidth', 2);
hold on;
plot(theta(91, :)*(180/pi), db(E_abs(91, :)./max(E_abs(91, :))), 'LineWidth', 2);

xlabel('\theta(Deg)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Norm GO field, Norm |E^{ff}_{feed}| [dB]', 'FontSize', 12, 'FontWeight', 'bold');
title('Absolute GO field incident and E^{ff} of feed at \phi = 90 ', 'FontSize', 12, 'FontWeight', 'bold');
legend({'Norm GO field','Norm |E^{ff}_{feed}|'},'Location','northeast', 'FontSize', 12, 'FontWeight', 'bold');

%% efficiency calculation 
V_far_th = Eth * focal_length./(exp(-1j * k0 * focal_length));
V_far_ph = Eph * focal_length./(exp(-1j * k0 * focal_length));

%V_far_th = Airy .* cos(theta) .* sin(phi);

%V_far_ph =  Airy .* cos(phi);

V_far_abs = sqrt(abs(V_far_th).^2 + abs(V_far_ph).^2);

C_rad = 1/(2 * zeta);

P_rad_int = C_rad * (V_far_abs).^2 .* sin(theta) * dth * dph;

P_rad = sum(sum(P_rad_int));

C_vg = (2/zeta);
C_go = exp(1j * k0 * focal_length)./(focal_length);

V_go_th = e_go_th ./ C_go;
V_go_ph = e_go_ph ./ C_go;

V_g_int_th = C_vg .* V_far_th(:,theta(1,:)<=theta_0) .* V_go_th .* sin(th);
V_g_int_ph = C_vg .* V_far_ph(:,theta(1,:)<=theta_0) .* V_go_ph .* sin(th);

V_g_int = V_g_int_th + V_g_int_ph;
V_g = sum(sum(V_g_int)) * dph * dth;


P_rx = abs(V_g).^2 ./ (16 .* P_rad);

A_r = pi * (d^2/4);

P_inc = E_0_pw.^2./(2 * zeta) * A_r;

eta_ap = (P_rx)./ P_inc;


