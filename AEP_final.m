%% Antenna dimesnions

dx = 20e-3;
dy = 20e-3;

w = 1e-3; % dipole width
l = 15e-3; % dipole length 

dth = pi/180;
dph = pi/180;
Vtm = 1;



th = [eps 0.5411]; % for 0 and 30 degrees
max_abs = zeros(size(th, 1), size(th, 2));

ph = eps;

f = 10e9;
    


for n = 1:size(th, 2)
    
    lambda = 3e8./f;
    k0 = 2 * pi ./ lambda;
    r_obs = 10000 * lambda;
    zeta = 120*pi;
    
   [mx, my] = meshgrid(-5:1:4, -5:1:4);
    
    %mx = 0;
    %my = 0;
    
    
    
    Constant_term = -1/(dx*dy);

    kx0 = k0 .* sin(th(n)) .* cos(ph);
    ky0 = k0 .* sin(th(n)) .* sin(ph);
    kz0 = k0 .* cos(th(n));

    kxm = kx0 - (2*pi*mx)/dx;
    kym = ky0 - (2*pi*my)/dy;
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
      
    [theta, phi] = meshgrid(-pi/2+eps:pi/180:pi/2+eps, eps:pi/180:2*pi+eps);
    
    kx = k0 * sin(theta) .* cos(phi);
    ky = k0 * sin(theta) .* sin(phi);
      
    I_kxm_1 = 2*k0.*(cos(kx * l/2) - cos(k0 * l/2))./((k0.^2 - kx.^2) * sin(k0 * l/2));

    J_kym_1 = sinc(ky*w/2/pi);
      
    kz = (-1j)*sqrt(-(k0.^2 - kx.^2 - ky.^2));

    c_new = (-zeta./(2 * k0 * kz));
      
    %% Calculation of the Dyad
    [Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz] = Dyad(k0, kx, ky, kz);

    %% Calculation of Spectral Green's function

    [SGFxx, SGFxy, SGFxz, SGFyx, SGFyy, SGFyz, SGFzx, SGFzy, SGFzz] = SGF(Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz, c_new);

    N_x = 10;
    N_y = 10;
      
    AF_x = 0;
    AF_y = 0;
      
    for i = 1:10
        AF_x = AF_x + exp(1j * (kx - kx0) * (i-1) * dx);
      
        AF_y = AF_y + exp(1j * (ky - ky0) * (i-1) * dy);
    end
      
    Kons = 1j * k0 * cos(theta) * exp(-1j * k0 * r_obs)/(2 * pi * r_obs);
    B = I_kxm_1 .* J_kym_1;
    
    ex = e_single_x .* AF_x .* AF_y .* I0;
    ey = e_single_y .* AF_x .* AF_y .* I0;
    ez = e_single_z .* AF_x .* AF_y .* I0;
    
    e_abs = sqrt(ex.^2 + ey.^2 +ez.^2);
    
    e_th = ex .* cos(theta) .* cos(phi) + ey .* cos(theta) .* sin(phi) - ez .* sin(theta);
    e_ph = -sin(phi) .* ex + cos(phi) .* ey;
    
    e_abs_2 = sqrt(e_th.^2 + e_ph.^2);
     
     
    max_abs(n) = abs(e_abs(1, 91));
    
    if n ~=1 
        str = '--';
    else 
        str = '-';
    end
     
     
    figure(2);
    hold on;
    
    plot(theta(1, :)*(180/pi), db(abs(e_abs(91, :))/max_abs(1)), str, 'LineWidth', 2);
    
    
%     plot(theta(1, :)*(180/pi), db(abs(e_ph(91, :))/max_ph), 'LineWidth', 1);

%     Array Factor:

   % AF = AF_x .* AF_y;
    
    %figure(2);
    %plot(theta(1, :)*(180/pi), db(abs(AF(1, :))/max(abs(AF(1, :)))), '--', 'LineWidth', 2);
    %hold on;
    %plot(theta(1, :)*(180/pi), db(abs(AF(91, :))/max(abs(AF(91, :)))), 'LineWidth', 1);
    
    xlabel('\theta(Deg), \theta_0(Deg)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('E_{abs}(\phi = 90) (dB)', 'FontSize', 12, 'FontWeight', 'bold');
    title('Far electric field pattern at scan angle 0 and 30 deg with AEP on H plane', 'FontSize', 12, 'FontWeight', 'bold');
    
    
%     xlabel('\theta(Deg)', 'FontSize', 12, 'FontWeight', 'bold');
%     ylabel('AF_x \times AF_y (dB)', 'FontSize', 12, 'FontWeight', 'bold');
%     title('Array factor', 'FontSize', 12, 'FontWeight', 'bold');
%     legend({'AF_{\theta}(\phi = 0)','AF_{\phi}(\phi = \pi/2)'},'Location','northeast', 'FontSize', 12, 'FontWeight', 'bold');
    
    
   %print('Arr_factor', '-depsc')
end  

legend({'scan \theta_0 = 0','scan \theta_0 = 30', 'AEP'},'Location','northeast', 'FontSize', 12, 'FontWeight', 'bold');

 grid on;
 ylim([-40 0]);

