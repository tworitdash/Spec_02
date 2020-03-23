%% Antenna dimesnions

dx = 20e-3;
dy = 20e-3;

w = 1e-3; % dipole width
l = 15e-3; % dipole length 

dth = pi/180;
dph = pi/180;
Vtm = 1;



th = -pi/2-eps:pi/180:pi/2-eps;

ph = pi/2;

f = 10e9;

AEP_x = zeros(size(th, 1), size(th, 2));
AEP_y = zeros(size(th, 1), size(th, 2));
AEP_z = zeros(size(th, 1), size(th, 2));
AEP_abs = zeros(size(th, 1), size(th, 2));
    

for p = 1:size(th , 2)
    
    lambda = 3e8./f;
    k0 = 2 * pi ./ lambda;
    r_obs = 10000 * lambda;
    zeta = 120*pi;
    
   [mx, my] = meshgrid(-5:1:4, -5:1:4);
    
    %mx = 0;
    %my = 0;
    
    Constant_term = -1/(dx*dy);

    kx = k0 .* sin(th(p)) .* cos(ph);
    ky = k0 .* sin(th(p)) .* sin(ph);
    kz = k0 .* cos(th(p));

    kxm = kx - (2*pi*mx)/dx;
    kym = ky - (2*pi*my)/dy;
    kzm = (-1j)*sqrt(-(k0.^2 - kxm.^2 - kym.^2));

    c = (-zeta./(2 * k0 * kzm)); %constant term in the equations

    %% Calculation of the Dyad
    [Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz] = Dyad(k0, kxm, kym, kzm);

    %% Calculation of Spectral Green's function

    [SGFxx, SGFxy, SGFxz, SGFyx, SGFyy, SGFyz, SGFzx, SGFzy, SGFzz] = SGF(Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz, c);

    %% Fourier transform of the current 
        
        
    I_kxm = 2*k0.*(cos(kxm * l/2) - cos(k0 * l/2))./((k0.^2 - kxm.^2) * sin(k0 * l/2));

    J_kym = sinc(kym*w/2/pi);
            
    zin = Constant_term .* SGFxx .* abs(I_kxm).^2 .* abs(J_kym).^2;
            
    z_in_f = sum(sum(zin));
      
    I0 = 1/(50 + z_in_f);
      
    I_kxm_1 = 2*k0.*(cos(kx * l/2) - cos(k0 * l/2))./((k0.^2 - kx.^2) * sin(k0 * l/2));

    J_kym_1 = sinc(ky*w/2/pi);
      
    AEP_x(p) = 1j .* kz .* SGFxx(6, 6) .* I_kxm_1 .* J_kym_1 .* exp(-1j * k0 .* r_obs)/(2 * pi * r_obs) .* I0;
    AEP_y(p) = 1j .* kz .* SGFyx(6, 6) .* I_kxm_1 .* J_kym_1 .* exp(-1j * k0 .* r_obs)/(2 * pi * r_obs) .* I0;
    AEP_z(p) = 1j .* kz .* SGFzx(6, 6) .* I_kxm_1 .* J_kym_1 .* exp(-1j * k0 .* r_obs)/(2 * pi * r_obs) .* I0;
    
    AEP_abs(p) = sqrt(abs(AEP_x(p)).^2 + abs(AEP_y(p)).^2 + abs(AEP_z(p)).^2);
end
th_deg = th*180/pi;
AEP_norm = db(AEP_abs/max(AEP_abs));

plot(th_deg, AEP_norm, 'LineWidth', 2);


hold on;
%plot(th_deg(61), AEP_norm(61), '*');
%hold on;
%plot(th_deg(121), AEP_norm(121), '*');

 xlabel('Scan angle \theta_0(Deg)', 'FontSize', 12, 'FontWeight', 'bold');
 ylabel('AEP (dB) \phi = \pi/2 rad', 'FontSize', 12, 'FontWeight', 'bold');
 title('AEP on H plane', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
ylim([-20 0]);
%print('AEP_90', '-depsc');
%legend({'AF_{\theta}(\phi = 0)','AF_{\phi}(\phi = \pi/2)'},'Location','northeast', 'FontSize', 12, 'FontWeight', 'bold');
    

    
    
 

