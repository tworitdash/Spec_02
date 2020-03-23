%% Antenna dimesnions

dx = 15e-3;
dy = 15e-3;

w = 3e-3; % dipole width
l = 14e-3; % dipole length 

dth = pi/180;
dph = pi/180;

f = 10.7e9;

th = [eps, pi/3];
ph = eps; % For E plane; change it to 90 for H plane
max_abs = zeros(size(th, 1), size(th, 2));

for j = 1:size(th, 2)


    
    [mx, my] = meshgrid(-10:1:10, -10:1:10);
    %mx = 0;
    %my = 0;
     
    lambda = 3e8./f;
    k0 = 2 * pi ./ lambda;
    zeta = 120*pi;

    Constant_term = -1/(dx*dy);

            kx0 = k0 .* sin(th(j)) .* cos(ph);
            ky0 = k0 .* sin(th(j)) .* sin(ph);
            kz0 = k0 .* cos(th(j));

            kxm = kx0 - (2*pi*mx)/dx;
            kym = ky0 - (2*pi*my)/dy;
            kzm = (-1j)*sqrt(-(k0.^2 - kxm.^2 - kym.^2));

            c = (-1./(2 * zeta * k0 * kzm)); %constant term in the equations

        %% Calculation of the Dyad
            [Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz] = Dyad(k0, kxm, kym, kzm);

        %% Calculation of Spectral Green's function

            [SGFxx, SGFxy, SGFxz, SGFyx, SGFyy, SGFyz, SGFzx, SGFzy, SGFzz] = SGF(Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz, c);

        %% Fourier transform of the current 

            I_kxm = 2*k0.*(cos(kxm .* l/2) - cos(k0 .* l/2))./((k0.^2 - kxm.^2) .* sin(k0 * l/2));

            J_kym = sinc(kym.*w/2/pi);
        
            %% Z
            yin = 4 * Constant_term .* SGFxx .* abs(I_kxm).^2 .* abs(J_kym).^2;
            
            y_in_f = sum(sum(yin));
            z_in_f = 1./y_in_f;
            
    
            y00 = yin(11, 11);
             
            yin_higher = y_in_f - y00;
            
            z00 = 1./y00;
            zin_higher = 1./yin_higher;
            
            r_obs = 1000 * lambda;
            
            
    I0 = 1/(y_in_f); %I_0 is actually V_0.. I was too lazy to change it everywhere...
                        % So, better don't complain :D
      
    [theta, phi] = meshgrid(-pi/2+eps:pi/180:pi/2+eps, eps:pi/180:2*pi+eps);
    
    kx = k0 * sin(theta) .* cos(phi);
    ky = k0 * sin(theta) .* sin(phi);
      
    I_kxm_1 = 2*k0.*(cos(kx * l/2) - cos(k0 * l/2))./((k0.^2 - kx.^2) * sin(k0 * l/2));

    J_kym_1 = sinc(ky*w/2/pi);
      
    kz = (-1j)*sqrt(-(k0.^2 - kx.^2 - ky.^2));

    c_new = (-1./(2 * zeta * k0 * kz));
      
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
    
    h_single_x = Kons .* SGFxx .* B;
    h_single_y = Kons .* SGFyx .* B;
    h_single_z = Kons .* SGFzx .* B;
    
    hx = h_single_x .* AF_x .* AF_y .* I0;
    hy = h_single_y .* AF_x .* AF_y .* I0;
    hz = h_single_z .* AF_x .* AF_y .* I0;
    
    h_abs = sqrt(hx.^2 + hy.^2 +hz.^2);
    
    max_abs(j) = abs(h_abs(1, 91));
    
    if j ~=1 
        str = '--';
    else 
        str = '-';
    end
     
     
    
    hold on;
    
    plot(theta(1, :)*(180/pi), db(abs(h_abs(1, :))/max_abs(1)), str, 'LineWidth', 2);
    %plot(theta(91, :)*(180/pi), db(abs(h_abs(1, :))), str, 'LineWidth', 2);
    
    
%     plot(theta(1, :)*(180/pi), db(abs(e_ph(91, :))/max_ph), 'LineWidth', 1);

%     Array Factor:

   % AF = AF_x .* AF_y;
    
    %figure(2);
    %plot(theta(1, :)*(180/pi), db(abs(AF(1, :))/max(abs(AF(1, :)))), '--', 'LineWidth', 2);
    %hold on;
    %plot(theta(1, :)*(180/pi), db(abs(AF(91, :))/max(abs(AF(91, :)))), 'LineWidth', 1);
    
    xlabel('\theta(Deg), \theta_0(Deg)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('E_{abs}(\phi = 0) (dB)', 'FontSize', 12, 'FontWeight', 'bold');
    title('Far magnetic field pattern at scan angle 0 and 60 deg with AEP on E plane', 'FontSize', 12, 'FontWeight', 'bold');
            
    grid on;




end

% xlabel('frequency(GHz)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Re Z_{in} , Im Z_{in} (\Omega) ', 'FontSize', 12, 'FontWeight', 'bold');
% title('Fundamental and higher order active input impedance at \theta = 60', 'FontSize', 12, 'FontWeight', 'bold');
% legend({'Real Z_{00}(60, 0)','Imag Z_{00}(60, 0)', 'Real Z_{h0}(60, 0)', 'Imag Z_{h0}(60, 0)'},'Location','northeast', 'FontSize', 12, 'FontWeight', 'bold');
% %ylim([-2000 3000]);
% print('Impedance_60_0h', '-dpng');
% %xlim([-180 180]);
% 


