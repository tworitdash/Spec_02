%%
f = 100e9;
lambda = 3e8/f;
k0 = 2 * pi / lambda;
eps_r = 1;
zeta = 120*pi;
%% Circular dimensions

%d = focal_length/5; %diameter of the reflector
% d = 50*lambda; %make it an array whenever needed

d = linspace(0.2, 2, 100);

% focal_length = 3 * d;


focal_length = 1;

eta_s = zeros(size(d));

n_f = zeros(size(d));
d_f = zeros(size(d));
etp_n_rho = zeros(size(d));
etp_n_phi = zeros(size(d));
etp_n = zeros(size(d));
etp_d = zeros(size(d));
e_tp = zeros(size(d));
D_max = zeros(size(d));
D_a = zeros(size(d));
G = zeros(size(d));

%D = d;
D = 3 * lambda; %Diameter of the feed
a = D/2;

for i = 1:size(d, 2)
%%

    r_obs = 10000 * lambda;
%l = lambda/2;
%w = lambda/20;

%%

dth = pi/1000; dph = pi/180;
% For entire region of theta and phi
[th, ph] = meshgrid(1e-7:dth:pi/2-dth*1e-7, 0:dph:2*pi); 


%[th, ph] = meshgrid(-pi/2:dth:pi/2, eps:dph:2*pi); 

    drho = d(i)/500; dphi = pi/ 180;

    [rho, phi] = meshgrid(eps:drho:d(i)/2, eps:dphi:2*pi);
    theta = 2 * atan(rho/(2 * focal_length));

%[rho, ph] = meshgrid(eps:drho:(D/2), eps:dphi:2*pi);

%th = 2 * atan(rho/(2 * focal_length));

% For polar plot at specific phi [edit the phi here based on the question]
%[th, ph] = meshgrid(1e-7:dth:pi-dth*1e-7, pi/2); 

%wavenumbers in all directions (x, y and z)




%% Current Distribution on the aperture:
    [Jx, Jy, Jz] = J_Circular(D, k0 , th);
    [Jx_rpz, Jy_rpz, Jz_rpz] = J_Circular(D, k0, theta);
    c3 = -2j .* k0 .* zeta .* (exp(-1j * k0 * r_obs))./(2 .* pi * r_obs);
    E_r = 0;
    E_th = c3 .* Jy .* cos(th) .* sin(ph);
    E_ph = c3 .* Jy .* cos(ph);
    
    E_rho = c3 .* Jy_rpz .* sin(theta);
    E_phi = c3 .* Jy_rpz .* cos(theta);
   

%[E_r, E_th, E_ph] = Polar_tr(E_abs, th, ph);

    E_abs_2 = sqrt((abs(E_r).^2 + abs(E_th).^2 + abs(E_ph).^2));

    Constant_term = (r_obs)./(exp(-1j * k0 * r_obs));

    E_r_without_phase = Constant_term .* E_r;
    E_th_without_phase = Constant_term .* E_th;
    E_ph_without_phase = Constant_term .* E_ph;
    
    E_rho_without_phase = Constant_term .* E_rho;
    E_phi_without_phase = Constant_term .* E_phi;
    
    %plot(th(1, :)*180/pi, db(E_abs_2(1, :)./max(E_abs_2(1, :))));

%% Spillover efficiency:
    C_spillover = 1/(2 * zeta);
    
    f_pattern_square = (abs(E_th_without_phase).^2 + abs(E_ph_without_phase).^2);
    U_feed = C_spillover .* f_pattern_square;
    
    f_hash = focal_length/d(i);
    
    theta_0 = 2 * acot(4 * f_hash);
    
    % numerator
    Int_u_n = U_feed(:,th(1,:)<=theta_0) .* sin(th(:,th(1,:)<=theta_0)) .* dth .* dph;
    n_f(i) = sum(sum(Int_u_n));
    
    Int_u_d = U_feed .* sin(th) .* dth .* dph;
    d_f(i) = sum(sum(Int_u_d));
        
        
    eta_s(i) = n_f(i)/d_f(i);
        
    
    %% Tapper efficiency:
    
    Int_etp_n_rho = E_rho_without_phase .* rho .* drho .* dphi;
    etp_n_rho(i) = sum(sum(Int_etp_n_rho));
    
    Int_etp_n_phi = E_phi_without_phase .* rho .* drho .* dphi;
    etp_n_phi(i) = sum(sum(Int_etp_n_phi));
    
    etp_n(i) = abs(etp_n_rho(i)).^2 + abs(etp_n_phi(i)).^2;
    
    E_abs_rpz = abs(E_rho_without_phase).^2 + abs(E_phi_without_phase).^2;
    Int_etp_d = E_abs_rpz .* rho .* drho .* dphi;
    etp_d(i) = sum(sum(Int_etp_d));
    
    
    Area = pi .* (d(i)^2)/4;
    
    e_tp(i) =  (1/Area) * etp_n(i)./etp_d(i);
    
    
    
    D_max(i) = (4 * pi * Area)/(lambda)^2;
    
    D_a(i) = D_max(i) * e_tp(i);
    G(i) = D_max(i) * e_tp(i) * eta_s(i);


   
    
end

e_ap = eta_s .* e_tp;
plot(d, (e_ap));
hold on;
plot(d, eta_s);
hold on;
plot(d, e_tp);

%plot(theta_obs(1, :)*180/pi, db(E_abs_airy(1, :)./max(E_abs_airy(1, :))));

%figure(2)
%surface(rho.*cos(ph), rho.*sin(ph), db(abs(J_cx)) , 'linestyle' , 'none' );
%figure(4);
%surface(rho.*cos(ph), rho.*sin(ph), db(abs(J_cy)) , 'linestyle' , 'none' );


