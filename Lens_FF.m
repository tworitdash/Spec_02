%%
f = 100e9;
lambda = 3e8/f;
k0 = 2 * pi / lambda;
eps_r = 1;
zeta = 120*pi;
a = 4 * lambda;
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

dth = pi/180; dph = dth;
% For entire region of theta and phi
[th, ph] = meshgrid(1e-7:dth:pi-dth*1e-7, 0:dph:2*pi); 

% For polar plot at specific phi [edit the phi here based on the question]
%[th, ph] = meshgrid(1e-7:dth:pi-dth*1e-7, pi/2); 


%% Far field

c2 = exp(-1j * k_d * r_obs) / (r_obs);
order = 2;


%% change to rho and phi
drho = D/1000; dphi = pi/ 180;

[rho, phi] = meshgrid(eps:drho:(D/2), eps:dphi:2*pi);

z = a * sqrt(1 - (rho/b).^2) + c;
theta = atan(rho./z);

[Eth, Ephi] = FF_Lens(c2, order, theta, phi);
Eth_without_phase = Eth./c2;
Ephi_without_phase = Ephi./c2;

E_abs = sqrt(abs(Eth).^2 + abs(Ephi).^2);


%% Directivity Calculation

zeta_d = zeta/sqrt(epsilon_r);

U_tot = E_abs.^2 ./ (2 * zeta_d);

%U_3d = U_tot(1, 1, :);
%U = squeeze(U_3d(1, :));

%P_rad = zeros(length(freq), 1);
dtheta = (theta(1, end) - theta(1, 1)) ./ size(theta, 2);


P_rad = Prad_Assign(U_tot, dtheta, dphi, theta, phi);

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

J_x = cj .*(T_parallel .* (Eth/c2) .* cos(phi) - T_perpendicular .* (Ephi/c2) .* sin(phi));
J_y = cj .*(T_parallel .* (Eth/c2) .* sin(phi) + T_perpendicular .* (Ephi/c2) .* cos(phi));

J_abs = sqrt(abs(J_x).^2 + abs(J_y).^2);

%% Plots

%plot(th, D);
%polarplot(th, E_abs/max(E_abs))
surface(rho.*cos(phi), rho.*sin(phi), 20*log10(abs(J_x)) , 'linestyle' , 'none' );
%surf(kx/k0, ky/k0, E_abs);
E_plot = Eth(1, :);
E_max =  max(Eth(1, :));
figure;
plot(theta(1, :)*(180/pi), db(E_plot./E_max));

%% Far field of the aperture:



