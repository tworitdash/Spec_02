%% Antenna dimesnions

% for dipole
clear all;
close all;

c0 = 3e8;
dth = pi/180;
dph = pi/180;

Nx = 7;
Dummy = 2;

f = 2e9:0.05e9:15e9;

lambda0 =c0./10^10;

w = 0.2 .* lambda0;
deld = 0.2 .* lambda0;
dx = 0.5 .* lambda0;
dy = 0.5 .* lambda0;
hd = 0.25 .* lambda0;

% w = 0.1 .* lambda0;
% deld = 0.1 .* lambda0;
% dx = 0.5 .* lambda0;
% dy = 0.5 .* lambda0;
% hd = 0.25 .* lambda0;



th = eps;
ph = eps; % For E plane; change it to 90 for H plane

for j = 1:size(th, 2)

y_in_f = zeros(size(f, 1), size(f, 2));
z_in_f = zeros(size(f, 1), size(f, 2));
yin_higher = zeros(size(f, 1), size(f, 2));
y00 = zeros(size(f, 1), size(f, 2));
zin = zeros(size(f, 1), size(f, 2));
yin = zeros(size(f, 1), size(f, 2));
z00 = zeros(size(f, 1), size(f, 2));
Gamma = zeros(size(f, 1), size(f, 2));

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
    
            c = (-zeta./(2 * k0 * kzm)); %constant term in the equations

        %% Calculation of the Dyad
            [Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz] = Dyad(k0, kxm, kym, kzm);

        %% Calculation of Spectral Green's function

            [SGFxx, SGFxy, SGFxz, SGFyx, SGFyy, SGFyz, SGFzx, SGFzy, SGFzz] = SGF(Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz, c);

       
            const = 1/dy;
            
            D_int = const .* SGFxx .* besselj(0, (kym .* w)/2) .* (1 - exp(-1j .* kzm .* 2 .* hd));
            
            D_inf = sum(D_int, 1);
            
    
            
            const2 = (-1/dx);
            
            yin_int = const2 .* sinc((kxm(1,:) .* deld)/2/pi).^2./(D_inf);
            
           
            
            yin(p) = sum(yin_int);
            
            Size = Nx + 2.*Dummy;
            yin_int_mutual = zeros(size(f, 2), Size);
            del = 0.01 .* k0;
            a = -50.*k0 - 1j.*del;
            b = 50.*k0 + 1j.*del;
            
            for m = 0:(Nx + 2*Dummy-1)
                    int =  @(kx) (-1/(2.*pi)) .* sinc((kx .* deld)/2/pi).^2./(D_inf_func(kx, th(j), ph, w, k0, dy, hd) .* exp(-1j .* kx .* abs(m)./dx));
                   if m == 0
                        yin_int_mutual(p, m+1) = integral(int, a, b, 'Waypoints', [(-1-1j).*del, (1+1j).*del]);
                   else
                       yin_int_mutual(p, m+1) = integral(int, -del-1j.*20.*k0, k0+del-1j*20.*k0, 'Waypoints', [(-1-1j).*del, (1+1j).*del, k0+1j.*del, k0+del]);
                   end
                   
            end
           v0 = 1;
            
           v = [zeros(1, Dummy), v0.*exp(-1j.*kx0.*dx).*(0:Nx-1), zeros(1, Dummy)].';
           z0 = 400;
           Y = toeplitz(real(yin_int_mutual(p, :)))+1j.*toeplitz(imag(yin_int_mutual(p, :)));
           I = eye(Size, Size);
           Z = inv(Y);
           i(p, :) = inv(Z + z0*I)*v;
           
           Z(p, :) = v./i(p, :).' - z0;
           
            %yin_int_mutual = const2 .* sinc((kxm(1,:) .* deld)/2/pi).^2./(D_inf);
            
            zin(p) = 1./yin(p);
            
            
            
            Gamma(p) = (zin(p) - z0) / (zin(p) + z0);
            
            
            
            
            
end

% figure(1);
% 
% plot(f/f(end), real(zin), 'LineWidth', 3);
% grid on;
% hold on;
% plot(f/f(end), imag(zin), 'LineWidth', 3);
% hold on;

% 
% figure(1);
% 
% plot(f*10^(-9), 20*log10(abs(Gamma)), 'LineWidth', 3);
% grid on;
% hold on;

end

% xlabel('f/f_0', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Re Z_{in} , Im Z_{in} (\Omega) ', 'FontSize', 12, 'FontWeight', 'bold');
% title('Active input impedance', 'FontSize', 12, 'FontWeight', 'bold');
% %legend({'Real Z_{in}(0, 0)','Imag Z_{in}(0, 0)', 'Real Z_{in}(45, 0)','Imag Z_{in}(45, 0)'},'Location','northeast', 'FontSize', 12, 'FontWeight', 'bold');
% legend({'Real Z_{in}','Imag Z_{in}'},'Location','northeast', 'FontSize', 12, 'FontWeight', 'bold');
% 
% %ylim([-2000 3000]);
% print('A5_Q1_dipole_br', '-depsc');
% %xlim([-180 180]);


% xlabel('Frequency(GHz)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Active \Gamma (dB)', 'FontSize', 12, 'FontWeight', 'bold');
% title(['Active reflection coefficient at z0 = ', num2str(z0), '\Omega'], 'FontSize', 12, 'FontWeight', 'bold');

%legend({'\Gamma (0, 0)', '\Gamma (45, 0) E Plane', '\Gamma (45, 90) H Plane'},'Location','northeast', 'FontSize', 12, 'FontWeight', 'bold');


