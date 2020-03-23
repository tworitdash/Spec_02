%% Antenna dimesnions

% for dipole
clear all;
close all;

c0 = 3e8;
dth = pi/180;
dph = pi/180;

f = 3e9:0.05e9:15e9;

lambda0 =c0./f(end);

w = 0.1 .* lambda0;
deld = 0.1 .* lambda0;
dx = 0.5 .* lambda0;
dy = 0.5 .* lambda0;



th = [eps];
ph = eps; % For E plane; change it to 90 for H plane

for j = 1:size(th, 2)

y_in_f = zeros(size(f, 1), size(f, 2));
z_in_f = zeros(size(f, 1), size(f, 2));
yin_higher = zeros(size(f, 1), size(f, 2));
y00 = zeros(size(f, 1), size(f, 2));
zin = zeros(size(f, 1), size(f, 2));
z00 = zeros(size(f, 1), size(f, 2));

%% wave properties 

    %f = 10e9;


   
%% Spectral Green's function
    

for p=1:size(f, 2)
    
    [mx, my] = meshgrid(-20:1:20, -30:1:30);
    %D_inf = zeros(size(mx));
    
 
     
        lambda = 3e8./f(p);
    
        k0 = 2 * pi ./ lambda;
    
        zeta = 120*pi;
   

  

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

       
            const = 1/dy;
            
            D_int = const .* SGFxx .* besselj(0, (kym .* w)/2);
            
            D_inf = sum(D_int, 1);
            
    
            
            const2 = (-1/dx);
            
            zin_int = const2 .* sinc((kxm(1,:) .* deld)/2/pi).^2./(D_inf);
            
            zin(p) = sum(zin_int);
            
            
            
            
            
end

figure(1);

plot(f/f(end), real(zin), 'LineWidth', 3);
grid on;
hold on;
plot(f/f(end), imag(zin), 'LineWidth', 3);
hold on;




end

xlabel('f/f_0', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Re Z_{in} , Im Z_{in} (\Omega) ', 'FontSize', 12, 'FontWeight', 'bold');
title('Active input impedance at \theta_0 = 0', 'FontSize', 12, 'FontWeight', 'bold');
legend({'Real Z_{in}','Imag Z_{in}'},'Location','northeast', 'FontSize', 12, 'FontWeight', 'bold');
print('A5_Q1_slot_br', '-depsc');

%ylim([-2000 3000]);
%print('Admittance_in', '-dpng');
%xlim([-180 180]);



