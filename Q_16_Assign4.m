%%
f = 100e9;
lambda = 3e8/f;
k0 = 2 * pi / lambda;
eps_r = 1;
zeta = 120*pi;
r_obs = 10000 * lambda;
%% Circular dimensions
f_number = 3;
D = 6 * lambda;  %diameter of feed


d = 50 * lambda; %diameter of reflector
focal_length = f_number * d;

%% Incident electric field
E_0_pw = 1;

dth = pi/2000; dph = pi/180;
% For entire region of theta and phi
theta_0 = 2 * atan(1/(4 * f_number));

[th, ph] = meshgrid(-theta_0:dth:theta_0, eps:dph:2*pi); 

%[th, ph] = meshgrid(eps:dth:pi-dth*eps, 0:dph:2*pi); 



e_go_th = - (2 ./ (1 + cos(th))) .* E_0_pw .* sin(ph);
e_go_ph = - (2 ./ (1 + cos(th))) .* E_0_pw .* cos(ph);

e_go_abs = sqrt(abs(e_go_th).^2 + abs(e_go_ph).^2);

%% Far field of feed 


%[theta, phi] = meshgrid(eps:dth:pi/2-eps*dth, eps:dph:2*pi); 

[theta, phi] = meshgrid(-theta_0:dth:theta_0, eps:dph:2*pi); 


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

plot(th(1, :)*(180/pi), db(abs(e_go_ph(1, :))./max(abs(e_go_ph(1, :)))), 'LineWidth', 2);
hold on;
plot(theta(1, :)*(180/pi), db(abs(Eph(1, :)./max(abs(Eph(1, :))))), 'LineWidth', 2);
grid on;
xlabel('\theta(Deg)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Norm E^{GO}_{\phi}, Norm E^{ff}_{\phi_{feed}} [dB]', 'FontSize', 12, 'FontWeight', 'bold');
title('\theta component of GO field incident and E^{ff} of feed at \phi = 0 ', 'FontSize', 12, 'FontWeight', 'bold');
legend({'E_{GO}_{\phi} \phi = 0','E^{ff}_{\phi_{feed}} \phi = 0'},'Location','northeast', 'FontSize', 12, 'FontWeight', 'bold');

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

%% Directivity

Denominator_int = P_rx_ang .* sin(del_theta).* (dth/2) .* dph;

Direc = (4 * pi) * P_rx_ang ./ sum(sum(Denominator_int));

Direc_max = max(Direc);

G = 10*log10(eta_ap * max(Direc));




%% Displaced field

dx = D;

Eth_disp = Eth .* exp(1j*kx.*dx);  % displaced field 
Eph_disp = Eph .* exp(1j*kx.*dx);  % displaced field 

V_far_th_disp = Eth_disp * focal_length./(exp(-1j * k0 * focal_length));
V_far_ph_disp = Eph_disp * focal_length./(exp(-1j * k0 * focal_length));

%% efficiency in case of displaced feed

V_far_abs_disp = sqrt(abs(V_far_th_disp).^2 + abs(V_far_ph_disp).^2);

C_rad = 1/(2 * zeta);

P_rad_int_disp = C_rad * (V_far_abs_disp).^2 .* sin(theta) * dth * dph;

P_rad_disp = sum(sum(P_rad_int_disp));



V_g_ang_disp = zeros(size(del_theta, 1), size(del_theta, 2));
P_rx_ang_disp = zeros(size(del_theta, 1), size(del_theta, 2));

for i = 1:size(del_theta, 1)
     for j = 1:size(del_theta, 2)
    
        phase = (focal_length./k0).*(kx_1 .* k0 .* sin(del_theta(i, j)) .* cos(phi_i(i, j)) + ky_1 .* k0 .* sin(del_theta(i, j)) .* sin(phi_i(i, j))); 
            
        e_ang_th = e_go_th .* exp(-1j * phase.*(1 + delta_n));
        e_ang_ph = e_go_ph .* exp(-1j * phase.*(1 + delta_n));
        
        V_go_th = e_ang_th ./ C_go;
        V_go_ph = e_ang_ph ./ C_go;

        V_g_int_th_disp = C_vg .* V_far_th_disp(:,theta(1,:)<=theta_0) .* V_go_th .* sin(th);
        V_g_int_ph_disp = C_vg .* V_far_ph_disp(:,theta(1,:)<=theta_0) .* V_go_ph .* sin(th);

        V_g_int_disp = V_g_int_th_disp + V_g_int_ph_disp;
        
        V_g_ang_disp(i, j) = sum(sum(V_g_int_disp)) * dph * dth;


        P_rx_ang_disp(i, j) = abs(V_g_ang_disp(i, j)).^2 ./ (16 .* P_rad_disp);

        
        
     end
    
end

%% efficiency



P_rx_disp = max(P_rx_ang_disp(1, :));

A_r = pi * (d^2/4);

P_inc = E_0_pw.^2./(2 * zeta) * A_r;

eta_ap_disp = (P_rx_disp)./ P_inc;

figure();
plot(del_theta(1, :)*180/pi, 10*log10(P_rx_ang(1, :)./max(P_rx_ang(1, :))));
hold on;
plot(del_theta(1, :)*180/pi, 10*log10(P_rx_ang_disp(1, :)./max(P_rx_ang_disp(1, :))));



