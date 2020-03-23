%%
f = 60e9;
lambda = 3e8/f;
k0 = 2 * pi / lambda;
eps_r = 1;
zeta = 120*pi;
%% Circular dimensions
focal_length = 1;
d = focal_length/5;
%D = d;
D = 3 * lambda;
a = D/2;


%%

r_obs = 10000 * lambda;
%l = lambda/2;
%w = lambda/20;

%%

%dth = pi/180; dph = dth;
% For entire region of theta and phi
%[th, ph] = meshgrid(1e-7:dth:pi-dth*1e-7, 0:dph:2*pi); 


%[th, ph] = meshgrid(-pi/2:dth:pi/2, eps:dph:2*pi); 

drho = d/500; dphi = pi/ 180;

[rho, ph] = meshgrid(eps:drho:d/2, eps:dphi:2*pi);
th = 2 * atan(rho/(2 * focal_length));

%[rho, ph] = meshgrid(eps:drho:(D/2), eps:dphi:2*pi);

%th = 2 * atan(rho/(2 * focal_length));

% For polar plot at specific phi [edit the phi here based on the question]
%[th, ph] = meshgrid(1e-7:dth:pi-dth*1e-7, pi/2); 

%wavenumbers in all directions (x, y and z)




%% Current Distribution on the aperture:
[Jx, Jy, Jz] = J_Circular(D, k0 , th);
c3 = -2j .* k0 .* zeta .* (exp(-1j * k0 * r_obs))./(2 .* pi * r_obs);
E_r = 0;
E_th = c3 .* Jy .* cos(th) .* sin(ph);
E_ph = c3 .* Jy .* cos(ph);


%[E_r, E_th, E_ph] = Polar_tr(E_abs, th, ph);

E_abs_2 = sqrt((abs(E_r).^2 + abs(E_th).^2 + abs(E_ph).^2));

Constant_term = (r_obs)./(exp(-1j * k0 * r_obs));

E_r_without_phase = Constant_term .* E_r;
E_th_without_phase = Constant_term .* E_th;
E_ph_without_phase = Constant_term .* E_ph;

%% Magnetic current
Cm = exp(-1j .* k0 * 2 * focal_length) .* (cos(th/2)).^2/focal_length;

M_rho = -E_ph_without_phase .* Cm;
M_ph = E_th_without_phase .* Cm;
M_z = 0;

% Cartesian coordinates 

M_cx = cos(ph) .* M_rho - sin(ph) .* M_ph;
M_cy = sin(ph) .* M_rho + cos(ph) .* M_ph;
M_cz = M_z;

% Electric current
J_rho = E_th_without_phase .* Cm ./zeta;
J_ph = E_ph_without_phase .* Cm ./zeta;
J_z = 0;

%cartesian coordinates 
J_cx = cos(ph) .* J_rho - sin(ph) .* J_ph;
J_cy = sin(ph) .* J_rho + cos(ph) .* J_ph;
J_cz = J_z;
J_abs = sqrt(abs(J_cx).^2 + abs(J_cy).^2 + abs(J_cz).^2);


%% Fourier Transform of Current distribution

%depends on the current distribution:
dth = pi/180; dph = dth;
%th_obs = linspace( -5*(lambda/d), 5*(lambda/d), 181) ;

th_obs = linspace( -5*(lambda/d), 5*(lambda/d), 181) ;
[theta_obs, phi_obs] = meshgrid(th_obs, eps:pi/4:2*pi);
%[J_ft_rho, J_ft_phi, J_ft_z] = J_Aperture(J_rho, J_ph, J_z, drho, dphi, k0, rho, ph, theta_obs, phi_obs);


    


[J_ft_x, J_ft_y, J_ft_z] = J_Aperture(J_cx, J_cy, J_cz, drho, dphi, k0, rho, ph, theta_obs, phi_obs);

kx = k0 .* sin(theta_obs) .* cos(phi_obs);
ky = k0 .* sin(theta_obs) .* sin(phi_obs);
kz = k0 .* cos(theta_obs);


c = (-zeta./(2 * k0 .* kz)); %constant term in the equations

[Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz] = Dyad(k0, kx, ky, kz);
[SGFxx, SGFxy, SGFxz, SGFyx, SGFyy, SGFyz, SGFzx, SGFzy, SGFzz] = SGF(Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz, c);

%% Far field

%J_ax = cos(phi_obs) .* J_ft_rho - sin(phi_obs) .* J_ft_phi;
%J_ay = sin(phi_obs) .* J_ft_rho + cos(phi_obs) .* J_ft_phi;
%J_az = J_ft_z;

c2 = 1j * kz * exp(-1j * k0 * r_obs) / (2 * pi * r_obs);

[Ex, Ey, Ez] = FF(c2, SGFxx, SGFxy, SGFxz, SGFyx, SGFyy, SGFyz, SGFzx, SGFzy, SGFzz, J_ft_x, J_ft_y, J_ft_z);

E_abs = sqrt((abs(Ex).^2 + abs(Ey).^2 + abs(Ez).^2));

%% Airy pattern

[J_airy_x, J_airy_y, J_airy_z] = J_Circular(d, k0 , theta_obs);


c2_airy = 1j * kz * exp(-1j * k0 * r_obs) / (2 * pi * r_obs);

[Ex_airy, Ey_airy, Ez_airy] = FF(c2_airy, SGFxx, SGFxy, SGFxz, SGFyx, SGFyy, SGFyz, SGFzx, SGFzy, SGFzz, J_airy_x, J_airy_y, J_airy_z);

E_abs_airy = sqrt((abs(Ex_airy).^2 + abs(Ey_airy).^2 + abs(Ez_airy).^2));



figure(1)
plot(theta_obs(1, :)*180/pi, db(E_abs(1, :)./max(E_abs(1, :))));
hold on;

plot(theta_obs(1, :)*180/pi, db(E_abs_airy(1, :)./max(E_abs_airy(1, :))));

%figure(2)
surface(rho.*cos(ph), rho.*sin(ph), db(abs(J_cx)) , 'linestyle' , 'none' );
%figure(4);
%surface(rho.*cos(ph), rho.*sin(ph), db(abs(J_cy)) , 'linestyle' , 'none' );


