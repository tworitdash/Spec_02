%% Antenna dimesnions

dx = 15e-3;
dy = 15e-3;

w = 3e-3; % dipole width
l = 14e-3; % dipole length 

dth = pi/180;
dph = pi/180;

f = 5e9:0.05e9:15e9;

th = [pi/3];
ph = eps; % For E plane; change it to 90 for H plane

for j = 1:size(th, 2)

y_in_f = zeros(size(f, 1), size(f, 2));
z_in_f = zeros(size(f, 1), size(f, 2));
yin_higher = zeros(size(f, 1), size(f, 2));
y00 = zeros(size(f, 1), size(f, 2));
zin_higher = zeros(size(f, 1), size(f, 2));
z00 = zeros(size(f, 1), size(f, 2));

%% wave properties 

    %f = 10e9;


   
%% Spectral Green's function
    

for p=1:size(f, 2)
    
    [mx, my] = meshgrid(-30:1:30, -30:1:30);
    %mx = 0;
    %my = 0;
     
    lambda = 3e8./f(p);
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
            
            y_in_f(p) = sum(sum(yin));
            z_in_f(p) = 1./y_in_f(p);
            
    
            y00(p) = yin(11, 11);
             
            yin_higher(p) = y_in_f(p) - y00(p);
            
            z00(p) = 1./y00(p);
            zin_higher(p) = 1./yin_higher(p);
            
            
            
           
            %ref(p) = (real(z_in_f(p)) - 60)/(real(z_in_f(p)) + 60);

end

figure(1);

plot(f*10^-9, real(y_in_f), 'LineWidth', 3);
grid on;
hold on;
plot(f*10^-9, imag(y_in_f), 'LineWidth', 3);
hold on;

% figure(2);
% 
% plot(f*10^-9, real(y00), 'LineWidth', 3);
% grid on;
% hold on;
% plot(f*10^-9, imag(y00), 'LineWidth', 3);
% 
% hold on;
% 
% plot(f*10^-9, real(yin_higher), 'LineWidth', 3);
% grid on;
% hold on;
% plot(f*10^-9, imag(yin_higher), 'LineWidth', 3);
% 
% hold on;

%f_ = 3e8/dx * (1 / (sin(th) + 1));


end

xlabel('frequency(GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Re Y_{in} , Im Y_{in} (S) ', 'FontSize', 12, 'FontWeight', 'bold');
title('Active input admittance at \theta_0 = 0', 'FontSize', 12, 'FontWeight', 'bold');
legend({'Real Y_{in}(0, 0)','Imag Y_{in}(0, 0)', 'Real Z_{h0}(60, 0)', 'Imag Z_{h0}(60, 0)'},'Location','northeast', 'FontSize', 12, 'FontWeight', 'bold');
%ylim([-2000 3000]);
print('Admittance_in', '-dpng');
%xlim([-180 180]);



