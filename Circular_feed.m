%%
f = 60e9;
lambda = 3e8/f;
k0 = 2 * pi / lambda;
eps_r = 1;
zeta = 120*pi;
%% Circular dimensions

D = 3 * lambda;
a = D/2;
focal_length = 1;
d = focal_length/0.5;

%%

r_obs = 10000 * lambda;
%l = lambda/2;
%w = lambda/20;

%%

dth = pi/180; dph = dth;
% For entire region of theta and phi
% [th, ph] = meshgrid(1e-7:dth:pi-dth*1e-7, 0:dph:2*pi); 


%[th, ph] = meshgrid(-pi/2:dth:pi/2, eps:dph:2*pi); 

drho = d/1000; dphi = pi/ 180;

[rho, ph] = meshgrid(-d/2:drho:d/2, eps:dphi:2*pi);
th = 2 * atan(rho/(2 * focal_length));

%[rho, ph] = meshgrid(eps:drho:(D/2), eps:dphi:2*pi);

%th = 2 * atan(rho/(2 * focal_length));

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

[Jx, Jy, Jz] = J_Circular(D, k0 , th);

%% Far field

c2 = 4j * exp(-1j * k0 * r_obs) / (2 * pi * r_obs).*kz;

[Ex, Ey, Ez] = FF(c2, SGFxx, SGFxy, SGFxz, SGFyx, SGFyy, SGFyz, SGFzx, SGFzy, SGFzz, Jx, Jy, Jz);

E_abs = sqrt((abs(Ex).^2 + abs(Ey).^2 + abs(Ez).^2));

%% Current Distribution on the aperture:

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

%% Far Field from the aperture:




%% Directivity Calculation
%U = E_abs.^2 ./ (2 * zeta);

%P_rad = Prad(U, dth, dph, th, ph);

%D = (4 * pi) * (U./P_rad);

%% Plots

%plot(th, D);
%polarplot(th, E_abs/max(E_abs))

%surf(kx/k0, ky/k0, E_abs);
% figure(1);
% plot(th(1, :)*180/pi, db(abs(E_abs_2(1, :))/max(E_abs_2(1, :))));
%figure(2);
%plot(th(1, :)*180/pi, db(abs(J_abs(1, :))/max(J_abs(1, :))));

%figure(2);
%plot(rho(1, :), db(abs(J_abs(1, :))./max(J_abs(1, :))));

%figure(2);
%plot(rho(1, :), db(abs(J_abs(1, :))));
% figure(3);
% surface(rho.*cos(ph), rho.*sin(ph), 20*log10(abs(Ex/Constant_term)) , 'linestyle' , 'none' );
% figure(4);
% surface(rho.*cos(ph), rho.*sin(ph), 20*log10(abs(Ey/Constant_term)) , 'linestyle' , 'none' );

figure(3);
%surface(rho.*cos(ph), rho.*sin(ph), db(abs(J_cx)) , 'linestyle' , 'none' );
%figure(4);
% %surface(rho.*cos(ph), rho.*sin(ph), db(abs(J_cy)) , 'linestyle' , 'none' );
figure(5);
surface(rho.*cos(ph), rho.*sin(ph), db(abs(M_cx)) , 'linestyle' , 'none' );
figure(6);
surface(rho.*cos(ph), rho.*sin(ph), db(abs(M_cy)) , 'linestyle' , 'none' );
figure(7);
surface(rho.*cos(ph), rho.*sin(ph), db(abs(Ex .* Constant_term .* Cm)) , 'linestyle' , 'none' );
figure(8);
surface(rho.*cos(ph), rho.*sin(ph), db(abs(Ey .* Constant_term .* Cm)) , 'linestyle' , 'none' );


%figure(5);
%plot(rho(1, :), 20*log10(abs(J_abs(1, :))));

%% Far field of the aperture current:

