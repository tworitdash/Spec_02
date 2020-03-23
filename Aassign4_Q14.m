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

dth = pi/2000; dph = pi/100;
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

%figure(1);
%plot(th(91, :)*(180/pi), db(e_go_abs(91, :)./max(e_go_abs(91, :))), 'LineWidth', 2);
%hold on;
%plot(theta(91, :)*(180/pi), db(E_abs(91, :)./max(E_abs(91, :))), 'LineWidth', 2);

%xlabel('\theta(Deg)', 'FontSize', 12, 'FontWeight', 'bold');
%ylabel('Norm GO field, Norm |E^{ff}_{feed}| [dB]', 'FontSize', 12, 'FontWeight', 'bold');
%title('Absolute GO field incident and E^{ff} of feed at \phi = 90 ', 'FontSize', 12, 'FontWeight', 'bold');
%legend({'Norm GO field','Norm |E^{ff}_{feed}|'},'Location','northeast', 'FontSize', 12, 'FontWeight', 'bold');

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


%% 1.4

[del_theta, phi_i] = meshgrid(-6*pi/180:dth/2:6*pi/180, eps:dph:2*pi);

%del_theta = eps:dth:6*pi/180;

V_g_ang = zeros(size(del_theta, 1), size(del_theta, 2));
P_rx_ang = zeros(size(del_theta, 1), size(del_theta, 2));

%phi_i = eps:dph:2*pi;

kx_1 = k0 .* sin(th) .* cos(ph);
ky_1 = k0 .* sin(th) .* sin(ph);
kz_1 = k0 .* cos(th);

delta_n = (1 - cos(th))./(1 + cos(th));
%e_ang_th = zeros(size(del_theta, 1), size(del_theta, 2));
%e_ang_ph = zeros(size(del_theta, 1), size(del_theta, 2));

for i = 1:size(del_theta, 1)
     for j = 1:size(del_theta, 2)
    
        phase = (focal_length./k0).*(kx_1 .* k0 .* sin(del_theta(i, j)) .* cos(phi_i(i, j)) + ky_1 .* k0 .* sin(del_theta(i, j)) .* sin(phi_i(i, j))); 
            
        e_ang_th = e_go_th .* exp(-1j * phase.*(1 + delta_n));
        e_ang_ph = e_go_ph .* exp(-1j * phase.*(1 + delta_n));
        
        V_go_th = e_ang_th ./ C_go;
        V_go_ph = e_ang_ph ./ C_go;

        V_g_int_th = C_vg .* V_far_th(:,theta(1,:)<=theta_0) .* V_go_th .* sin(th);
        V_g_int_ph = C_vg .* V_far_ph(:,theta(1,:)<=theta_0) .* V_go_ph .* sin(th);

        V_g_int = V_g_int_th + V_g_int_ph;
        
        V_g_ang(i, j) = sum(sum(V_g_int)) * dph * dth;


        P_rx_ang(i, j) = abs(V_g_ang(i, j)).^2 ./ (16 .* P_rad);


        
     end
    
end


figure(2);
plot(del_theta(51, :)*180/pi, 10*log10(P_rx_ang(51, :)./max(P_rx_ang(51, :))));

hold on;
grid on;

%% Assignment 2

%%
f = 100e9;
lambda = 3e8/f;
k0 = 2 * pi / lambda;
eps_r = 1;
zeta = 120*pi;
%% Circular dimensions

d = 50 * lambda;
focal_length = 3 * d;
%D = d;
D = 6 * lambda;
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

drho = d/2000; dphi = pi/100;

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
dth = pi/2000; dph = pi/100;
%th_obs = linspace( -5*(lambda/d), 5*(lambda/d), 181) ;


[theta_obs, phi_obs] = meshgrid(-5*(lambda/d):dth:5*(lambda/d), eps:dph:2*pi);
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




plot(theta_obs(51, :)*180/pi, db(E_abs(51, :)./max(E_abs(51, :))));


%plot(theta_obs(1, :)*180/pi, db(E_abs_airy(1, :)./max(E_abs_airy(1, :))));

%figure(2)
%surface(rho.*cos(ph), rho.*sin(ph), db(abs(J_cx)) , 'linestyle' , 'none' );
%figure(4);
%surface(rho.*cos(ph), rho.*sin(ph), db(abs(J_cy)) , 'linestyle' , 'none' );



xlabel('\theta(Deg)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Received Power(dB), Radiation Pattern', 'FontSize', 12, 'FontWeight', 'bold');
title('Received power at the feed and radiation pattern from tx approach at \phi = 90 ', 'FontSize', 12, 'FontWeight', 'bold');
legend({'Rx -> FO','Tx'},'Location','northeast', 'FontSize', 12, 'FontWeight', 'bold');






