%% Antenna dimesnions

dx = 15e-3;
dy = 15e-3;

w = 1e-3; % dipole width
l = 14e-3; % dipole length 

dth = pi/180;
dph = pi/180;
Vtm = 1;

e_inc = Vtm / sqrt((dx * dy));

%th = -pi:dth:pi;
th = eps;
ph = eps; % For E plane; change it to 90 for H plane

z_in_f = zeros(size(th, 1), size(th, 2));
i_bf = zeros(size(th, 1), size(th, 2));

e_scat = zeros(size(th, 1), size(th, 2));

Gamma = zeros(size(th, 1), size(th, 2));
Tau = zeros(size(th, 1), size(th, 2));
%% wave properties 

    f = 6e9:0.01e9:14e9;
    

%% Spectral Green's function
    

for p=1:size(f, 2)
    
    lambda = 3e8./f(p);
    k0 = 2 * pi ./ lambda;
    zeta = 120*pi;
    
   [mx, my] = meshgrid(-10:1:10, -10:1:10);
    
    %mx = 0;
    %my = 0;
    
    
    
    Constant_term = -1/(dx*dy);

            kx = k0 .* sin(th) .* cos(ph);
            ky = k0 .* sin(th) .* sin(ph);
            kz = k0 .* cos(th);

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
            
            z_in_f(p) = sum(sum(zin));
            
            I_kxm_1 = 2*k0.*(cos(-kx * l/2) - cos(k0 * l/2))./((k0.^2 - kx.^2) * sin(k0 * l/2));

            J_kym_1 = sinc(-ky*w/2/pi);
            
            B = I_kxm_1 .* J_kym_1;
        
            %% Z
            
            %ref(p) = (real(z_in_f(p)) - 60)/(real(z_in_f(p)) + 60);
            
            
            
   
            
            v = e_inc .* B;
            
            i_bf(p) = v/z_in_f(p);
            
            e_scat(p) = (i_bf(p)./(dx * dy)) * SGFxx(11, 11) .* B;
            
            Gamma(p) = e_scat(p) ./ e_inc;
            Tau(p) = 1 + Gamma(p);

end


plot(f, abs(Gamma), 'LineWidth', 2);
hold on;
plot(f, abs(Tau), 'LineWidth', 2);
grid on;
plot(f, abs(Gamma).^2 + abs(Tau).^2, 'LineWidth', 2);




% plot(th*180/pi, real(z_in_f), 'LineWidth', 3);
% grid on;
% hold on;
% plot(th*180/pi, imag(z_in_f), 'LineWidth', 3);
% 
% 
% xlabel('\theta(Deg)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Re Z_{in} , Im Z_{in} (\Omega) ', 'FontSize', 12, 'FontWeight', 'bold');
% title('higher order Real and Imag part of Z_{in} at \phi = \pi/2 (H plane) dx=dy=2\lambda/3', 'FontSize', 12, 'FontWeight', 'bold');
% legend({'Real Z_{in}','Imag Z_{in}'},'Location','northeast', 'FontSize', 12, 'FontWeight', 'bold');
% 
% xlim([-180 180]);
% ylim([-50 1000]); 



%


%% Grating lobe diagrams for both the configuration:


        

