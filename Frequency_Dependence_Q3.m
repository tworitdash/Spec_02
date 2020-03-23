%%
[th, ph, h] = meshgrid(1e-7:dth:pi-dth*1e-7, 0:dph:2*pi, 1e-3:0.25e-3:20e-3); 
f = 15e9;
he = 1e-3:0.25e-3:20e-3;
%f = 5e9:1e9:15e9;
lambda = 3e8./f;
k0 = 2 * pi ./ lambda;
eps_r = 1;
zeta = 120*pi;

%%

r_obs = 100000 * lambda;
l = 10e-3;
w = 1e-3;

%%

dth = pi/180; dph = dth;
% For entire region of theta and phi
%[th, ph] = meshgrid(1e-7:dth:pi-dth*1e-7, 0:dph:2*pi); 
%th = 0;
%ph = pi/2;
% For polar plot at specific phi [edit the phi here based on the question]
%[th, ph] = meshgrid(1e-7:dth:pi-dth*1e-7, pi/2); 

%wavenumbers in all directions (x, y and z)
%kx = zeros(361,180,length(freq));
%for i = 1:length(freq)
 %   kx(:,:,i) = k0(i) .* sin(th(:, :, i)) .* cos(ph(:, :, i));
  %  ky = k0(i) .* sin(th(:, :, i)) .* sin(ph(:, :, i));
   % kz = k0(i) .* cos(th(:, :, i));
    %c(i) = (-zeta./(2 * k0(i) .* kz));
%end


    kx = k0 .* sin(th) .* cos(ph);
    ky = k0 .* sin(th) .* sin(ph);
    kz = k0 .* cos(th);
    c = (-zeta./(2 * k0 .* kz));


%c = (-zeta./(2 * k0 .* kz)); %constant term in the equations

%% Calculation of the Dyad
[Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz] = Dyad(k0, kx, ky, kz);

%% Calculation of Spectral Green's function

[SGFxx, SGFxy, SGFxz, SGFyx, SGFyy, SGFyz, SGFzx, SGFzy, SGFzz] = SGF(Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz, c);

%% Fourier Transform of Current distribution

%depends on the current distribution:

[Jx, Jy, Jz] = JFT_shifted_PEC(l, w, h, k0, kx, ky, kz);

%% Far field

c2 = 1j .* kz .* exp(-1j .* k0 .* r_obs) ./ (2 .* pi .* r_obs);

[Ex, Ey, Ez] = FF(c2, SGFxx, SGFxy, SGFxz, SGFyx, SGFyy, SGFyz, SGFzx, SGFzy, SGFzz, Jx, Jy, Jz);

E_abs = sqrt(abs(Ex.^2 + Ey.^2 + Ez.^2));

%% Directivity Calculation
U_tot = E_abs.^2 ./ (2 * zeta);
U_3d = U_tot(1, 1, :);
U = squeeze(U_3d(1, :));

P_rad = zeros(length(he), 1);

P_rad = Prad(P_rad, U_tot, dth, dph, th, ph);

D = (4 * pi) * (U./P_rad');

%D_theta_0 = squeeze(D(:, 1, 1:length(freq)));

%% Plots
plot(he/10^(-3), D)

%plot(th, D);
%polarplot(th, E_abs/max(E_abs))

%surf(kx/k0, ky/k0, E_abs);




