%% Antenna dimesnions

dx = 20e-3;
dy = 20e-3;

w = 1e-3; % dipole width
l = 14e-3; % dipole length 

dth = pi/180;
dph = pi/180;

th = -pi:dth:pi;
ph = 0; % For E plane; change it to 90 for H plane

z_in_f = zeros(size(th, 1), size(th, 2));
%% wave properties 

    f = 10e9;


    
    lambda = 3e8/f;
    k0 = 2 * pi / lambda;
    zeta = 120*pi;

%% Spectral Green's function
    

for p=1:size(th, 2)
    
    [mx, my] = meshgrid(-10:1:10, -10:1:10);
    %my = -5:1:-1;
    %mx = [-10 -9 -8 -7 -6 -5 -4 -3 -2 -1 1 2 3 4 5 6 7 8 8 10];
    %my = [-10 -9 -8 -7 -6 -5 -4 -3 -2 -1 1 2 3 4 5 6 7 8 8 10];
    
    
    
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
        
            %% Z
            zin = Constant_term .* SGFxx .* abs(I_kxm).^2 .* abs(J_kym).^2;
            z_in_f(p) = sum(sum(zin));
    
            %ref(p) = (real(z_in_f(p)) - 60)/(real(z_in_f(p)) + 60);

end



plot(th*180/pi, real(z_in_f), 'LineWidth', 3);
grid on;
hold on;
plot(th*180/pi, imag(z_in_f), 'LineWidth', 3);


xlabel('\theta(Deg)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Re Z_{in} , Im Z_{in} (\Omega) ', 'FontSize', 12, 'FontWeight', 'bold');
title('higher order Real and Imag part of Z_{in} at \phi = \pi/2 (H plane) dx=dy=2\lambda/3', 'FontSize', 12, 'FontWeight', 'bold');
legend({'Real Z_{in}','Imag Z_{in}'},'Location','northeast', 'FontSize', 12, 'FontWeight', 'bold');

xlim([-180 180]);
ylim([-50 1000]); 



%


%% Grating lobe diagrams for both the configuration:

mx1 = -2:1:2;
my1 = -2:1:2;

ang=0:0.01:2*pi;
count = 0;

kx0 = k0 * cos(ph);
ky0 = k0 * sin(ph);

for q = 1:size(mx1, 2)
    for s = 1:size(my1, 2)
        x1 = 2 * pi *mx1(q) / (dx);
        y1 = 2 * pi * my1(s) / (dy);
        
        k_v = sqrt(abs(kx0).^2 + abs(ky0).^2);
         
        xp =  k0 * cos(ang);
        yp =  k0 * sin(ang);
        
        %plot(x1 + xp, y1 + yp, 'LineWidth', 2);
        %hold on;
        %quiver(x1, y1, kx0, ky0, 'LineWidth', 2);
        
    end
end
grid on;
% xlabel('k_{xm}', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('k_{ym}', 'FontSize', 12, 'FontWeight', 'bold');
% title('Grating lobe diagram when dx = dy =2\lambda/3 at \phi = 0 (E Plane)', 'FontSize', 12, 'FontWeight', 'bold');
%legend({'Real Z_{in}','Imag Z_{in}'},'Location','northeast', 'FontSize', 12, 'FontWeight', 'bold');
        

