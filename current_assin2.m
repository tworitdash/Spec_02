%%
f = 60e9;
lambda = 3e8/f;
k0 = 2 * pi / lambda;
eps_r = 1;
zeta = 120*pi;
%% Circular dimensions

%d = 3 * lambda;

D = 3 * lambda;
focal_length = 1;
d = focal_length/0.5;
a = D/2;

%%

r_obs = 10000 * lambda;
%l = lambda/2;
%w = lambda/20;

%%

%dth = pi/180; dph = dth;
% For entire region of theta and phi
%[th, ph] = meshgrid(1e-7:dth:pi-dth*1e-7, 0:dph:2*pi); 
%[th, ph] = meshgrid(-pi/2:dth:pi/2, 0:dph:2*pi); 

drho = d/1000; dphi = pi/ 180;

[rho, ph] = meshgrid((-d/2):drho:(d/2), eps:dphi:2*pi);

th = 2 * atan(rho/(2 * focal_length));

% For polar plot at specific phi [edit the phi here based on the question]
%[th, ph] = meshgrid(1e-7:dth:pi-dth*1e-7, pi/2); 

%wavenumbers in all directions (x, y and z)

kx = k0 .* sin(th) .* cos(ph);
ky = k0 .* sin(th) .* sin(ph);
kz = k0 .* cos(th);


c = (-zeta./(2 * k0 .* kz)); %constant term in the equations

%% Calculation of the Dyad
[Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz] = Dyad(k0, kx, ky, kz);

%% Calculation of Spectral Green's function

[SGFxx, SGFxy, SGFxz, SGFyx, SGFyy, SGFyz, SGFzx, SGFzy, SGFzz] = SGF(Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz, c);


%% Fourier Transform of Current distribution

%depends on the current distribution:

[J_r, J_theta, J_phi] = J_Cyl(D, k0 , th, ph);

%% Far field

c2 = 1j * k0 *zeta * exp(-1j * k0 * r_obs) / (2 * pi * r_obs).*kz;

[E_r, E_theta, E_phi] = FF(c2, SGFxx, SGFxy, SGFxz, SGFyx, SGFyy, SGFyz, SGFzx, SGFzy, SGFzz, J_r, J_theta, J_phi);

E_abs = sqrt((abs(E_r).^2 + abs(E_theta).^2 + abs(E_phi).^2));




%% Current distribution: (In cylinderical coordinates)

Cm = exp(-2j * k0 * focal_length) .* (cos(th/2)).^2 ./ focal_length;

M_rho = - (E_phi./c2) .* exp(-2j * k0 * focal_length);
M_z = 0;
M_ph = E_rho .* exp(-2j * k0 * focal_length);

%% Current distribution: (In cartesian coordinates):

M_cx = cos(th) .* M_rho - sin(th) .* M_ph;
M_cy = sin(th) .* M_rho + cos(th) .* M_ph;
M_cz = M_z;

J_rho = E_rho/zeta;
J_phi = E_phi/zeta;
J_z = 0;

J_cx = cos(th) .* J_rho - sin(th) .* J_phi;
J_cy = sin(th) .* J_rho + cos(th) .* J_phi;
J_cz = J_z;

J_abs = sqrt((J_rho).^2 + (J_phi).^2 + (J_z).^2);





E_cx = E_abs .* cos(th);
E_cy = E_abs .* sin(th);



%% Directivity Calculation
%U = E_abs.^2 ./ (2 * zeta);

%P_rad = Prad(U, dth, dph, th, ph);

%D = (4 * pi) * (U./P_rad);

%% Plots

%plot(th, D);
%polarplot(th, E_abs/max(E_abs))

%surf(kx/k0, ky/k0, E_abs);
figure(1);
%plot(th(1, :)*180/pi, db(abs(J_abs(1, :))/max(abs(J_abs(1, :)))));
plot(rho(1, :), db(abs(E_abs(1, :))/max(abs(E_abs(1, :)))));

figure(2);
surface(rho.*cos(ph), rho.*sin(ph), 20*log10(abs(M_cx)) , 'linestyle' , 'none' );
figure(3);
surface(rho.*cos(ph), rho.*sin(ph), 20*log10(abs(M_cy)) , 'linestyle' , 'none' );



