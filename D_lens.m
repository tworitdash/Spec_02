%%
f = 100e9;
lambda = 3e8/f;
k0 = 2 * pi / lambda;
eps_r = 1;
zeta = 120*pi;
a =  4 * lambda;
%% Ellipse dimensions
epsilon_r = 11.9; 
e = 1/sqrt(epsilon_r);
c = a * e;
b = sqrt(a^2 - c^2);
D = 2 * b;
k_d = sqrt(epsilon_r) * k0;
%%

r_obs = 10000 * lambda;
%l = lambda/2;
%w = lambda/20;

%%

dth = pi/180; dph = dth * 72;
% For entire region of theta and phi
%[theta, phi] = meshgrid(1e-7:dth:pi-dth*1e-7, 0:dph:2*pi); 


% For polar plot at specific phi [edit the phi here based on the question]
%[th, ph] = meshgrid(1e-7:dth:pi-dth*1e-7, pi/2); 


%% Far field

c2 = exp(-1j * k_d * r_obs) / (r_obs);
order = 2;


%% change to rho and phi
drho = D/100; dphi = pi/ 180;

[rho, phi] = meshgrid(eps:drho:(D/2), eps:dphi:2*pi);

z = a * sqrt(1 - (rho/b).^2) + c;
theta = atan(rho./z);

[Eth, Ephi] = FF_Lens(c2, order, theta, phi);
E_abs = sqrt(abs(Eth).^2 + abs(Ephi).^2);


%% Directivity Calculation

zeta_d = zeta/sqrt(epsilon_r);

U_tot = E_abs.^2 ./ (2 * zeta_d);

%U_3d = U_tot(1, 1, :);
%U = squeeze(U_3d(1, :));

%P_rad = zeros(length(freq), 1);
%dtheta = (theta(1, end) - theta(1, 1)) ./ size(theta, 2);


%P_rad = Prad_Assign(U_tot, dtheta, dphi, theta, phi);

%D = (4 * pi) * (U./P_rad');

%% current distribution:

%transmission co-efficients
theta_i = acos((1 - e.*cos(theta))./(sqrt(1 + e^2 - 2*e.*cos(theta))));
theta_t = asin(sqrt(epsilon_r) .* sin(theta_i));

T_parallel = (2 * zeta .* cos(theta_i))./(zeta .* cos(theta_t) + zeta_d .* cos(theta_i));
T_perpendicular = (2 * zeta .* cos(theta_i))./(zeta_d .* cos(theta_t) + zeta .* cos(theta_i));


% Current distribution

r_ellipse = (a.*(1 - e^2))./(1 - e .* cos(theta));
S_theta = sqrt((cos(theta_t)./cos(theta_i)) .* ((e .* cos(theta)) - 1)./(e - cos(theta)));

cj = -(2/zeta) .* (S_theta./r_ellipse);

c3 = (r_obs) ./ exp(-1j * k_d * r_obs) ;

Eth_without_phase = Eth.*c3;
Eph_without_phase = Ephi.*c3;

%E_abs_without_phase = sqrt(abs(Eth_without_phase).^2 + abs(Eph_without_phase).^2);

J_x = cj .*(T_parallel .* (Eth_without_phase) .* cos(phi) - T_perpendicular .* (Eph_without_phase) .* sin(phi));
J_y = cj .*(T_parallel .* (Eth_without_phase) .* sin(phi) + T_perpendicular .* (Eph_without_phase) .* cos(phi));
J_z = cj .* 0;
J_abs = sqrt(abs(J_x).^2 + abs(J_y).^2);

%% Fourier Transform of the current distribution:

[theta_obs, phi_obs] = meshgrid(eps:dth:pi/2, eps:dph:2*pi);

[J_ft_x, J_ft_y, J_ft_z] = J_Aperture(J_x, J_y, J_z, drho, dphi, k0, rho, phi, theta_obs, phi_obs);

kx = k0 .* sin(theta_obs) .* cos(phi_obs);
ky = k0 .* sin(theta_obs) .* sin(phi_obs);
kz = k0 .* cos(theta_obs);


Constant_Term = (-zeta./(2 * k0 .* kz)); %constant term in the equations

[Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz] = Dyad(k0, kx, ky, kz);
[SGFxx, SGFxy, SGFxz, SGFyx, SGFyy, SGFyz, SGFzx, SGFzy, SGFzz] = SGF(Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz, Constant_Term);

%% Far field

%J_ax = cos(phi_obs) .* J_ft_rho - sin(phi_obs) .* J_ft_phi;
%J_ay = sin(phi_obs) .* J_ft_rho + cos(phi_obs) .* J_ft_phi;
%J_az = J_ft_z;

c_e = 1j * kz * exp(-1j * k0 * r_obs) / (2 * pi * r_obs);

[Ex, Ey, Ez] = FF(c_e, SGFxx, SGFxy, SGFxz, SGFyx, SGFyy, SGFyz, SGFzx, SGFzy, SGFzz, J_ft_x, J_ft_y, J_ft_z);

E_abs_ap = sqrt((abs(Ex).^2 + abs(Ey).^2 + abs(Ez).^2));

%% Directivity of Lens
Clens = (r_obs) ./ exp(-1j * k0 * r_obs);

E_lens_x_without_phase = Clens .* Ex;
E_lens_y_without_phase = Clens .* Ey;

E_lens_abs_without_phase = sqrt(abs(E_lens_x_without_phase).^2 + abs(E_lens_y_without_phase).^2);

U_lens = (1/(2*zeta)) * (E_lens_abs_without_phase).^2;

P_rad_i = U_lens .* sin(theta_obs);
    
P_rad = sum(sum(P_rad_i)) * dth * dph;


D_lens_ = 4*pi * U_lens ./ P_rad;

% P_rad of the feed

c_feed = exp(-1j * k_d * r_obs) / (r_obs);
order = 2;

[Eth_feed, Ephi_feed] = FF_Lens(c_feed, order, theta_obs, phi_obs);
E_abs_feed = sqrt(abs(Eth_feed).^2 + abs(Ephi_feed).^2);

E_abs_ = E_abs_feed / c_feed;
E_abs_without_phase = abs(E_abs_);

U_feed = (1/(2 *zeta_d)) .* (E_abs_without_phase).^2;

P_rad_feed_i = U_feed .* sin(theta_obs);
P_rad_feed = sum(sum(P_rad_feed_i)) .* dth .* dph;

eta_ = P_rad/P_rad_feed;

G_lens = D_lens_ .* eta_;


%% Plots
figure(1);

plot(theta_obs(1, :)*180/pi, 10*log10(D_lens_(1, :)));

hold on;

plot(theta_obs(1, :)*180/pi, 10*log10(G_lens(1, :)));




