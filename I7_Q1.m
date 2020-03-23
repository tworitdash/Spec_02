%% Antenna dimesnions

% for dipole
clear all;
close all;

c0 = 3e8;
dth = pi/180;
dph = pi/180;
er = 2.2;

Nx = 7;
Dummy = 2;

%f = 2e9:0.05e9:15e9;
f = 10e9;
%f = linspace(2e9, 15e9, 1000);
f0 = 10e9;

lambda0 = c0./f0;
lambda_ = lambda0 ./ sqrt(er);


% w = 0.2 .* lambda0;
% deld = 0.25 .* lambda0;
% dx = 0.5 .* lambda0;
% dy = 0.5 .* lambda0;
% hd = 0.31 .* lambda0;

w = 0.2 .* lambda0;
deld = 0.25 .* lambda0;
dx =  0.5 * lambda0;
dy =  0.5 * lambda0;
hd = 0.25 .* lambda_;



th = [-pi/2-eps:dth:pi/2-eps];
%th = [eps, pi/6, pi/4, pi/3, pi/2];
%th = eps;

ph = eps; % For E plane; change it to 90 for H plane

Dir = zeros(size(th, 1), size(th, 2));
Sinc_term_j = zeros(size(th, 1), size(th, 2));
AF_x_j = zeros(size(th, 1), size(th, 2));

for j = 1:size(th, 2)

y_in_f = zeros(size(f, 1), size(f, 2));
z_in_f = zeros(size(f, 1), size(f, 2));
yin_higher = zeros(size(f, 1), size(f, 2));
y00 = zeros(size(f, 1), size(f, 2));
zin = zeros(size(f, 1), size(f, 2));
yin = zeros(size(f, 1), size(f, 2));
z00 = zeros(size(f, 1), size(f, 2));
Gamma = zeros(size(f, 1), size(f, 2));
Zn = zeros(Nx + 2.*Dummy, size(f, 2));
Gamma_ = zeros(Nx + 2.*Dummy, size(f, 2));

%% wave properties 

    %f = 10e9;


   
%% Spectral Green's function
    

for p=1:size(f, 2)
    
    [mx, my] = meshgrid(-30:1:30, -30:1:30);
    %D_inf = zeros(size(mx));
    
       
        lambda = 3e8./f(p);
    
        k0 = 2 * pi ./ lambda;
      
            kx0 = k0 .* sin(th(j)) .* cos(ph);
            ky0 = k0 .* sin(th(j)) .* sin(ph);
            %kz0 = k0 .* cos(th(j));
        
        yin_int_mutual = zeros(Nx + 2.*Dummy, 1);
        v0 = 1;
        
        v = [zeros(1,Dummy) v0*exp(-1j*kx0*[0:Nx-1]*dx) zeros(1,Dummy)].';
        %v = (sin([eps:Nx+2*Dummy - 1]./(Nx+2*Dummy - 1).*pi).^2 .* exp(-1j*kx0*[-Dummy:Nx+Dummy-1]*dx)).';
    
            zeta = 120*pi;
   
            kxm = kx0 - (2*pi*mx)/dx;
            kym = ky0 - (2*pi*my)/dy;
            
            krho = sqrt(kxm.^2 + kym.^2);
            %k_rho_2 =  krho.^2;
            
            
            
            kz0 = (-1j)*sqrt(-(k0.^2 - krho.^2));
    
            z = hd+eps;
            
            zeta_s = zeta./sqrt(er); 
            
            
            ksub = k0 .* sqrt(er); 
            
            kzs = (-1j)*sqrt(-(ksub.^2 - krho.^2));

            % For TM
            z0_TM = zeta * kz0 ./ k0;
            zs_TM = zeta_s * kzs ./ ksub;

            % For TE

            z0_TE = zeta * k0 ./ kz0;
            zs_TE = zeta_s * ksub ./ kzs;

            [zup_TE, zdn_TE, zup_TM, zdn_TM] = zupdn(z0_TE, zs_TE, z0_TM, zs_TM, hd, kzs);


            [vtm, vte, itm, ite] = txline(zup_TE, zdn_TE, zup_TM, zdn_TM, hd, z, kz0, kzs, z0_TE, zs_TE, z0_TM, zs_TM);


            [SGFxx, ~, ~] = Green(vtm, vte, itm, ite, kxm, kym, krho.^2, zeta, k0);
            
           % kzm = (-1j)*sqrt(-(k0.^2 - kxm.^2 - kym.^2));
    
%             c = (-zeta./(2 * k0 * kzm)); %constant term in the equations
% 
%         %% Calculation of the Dyad
%             [Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz] = Dyad(k0, kxm, kym, kzm);
% 
%         %% Calculation of Spectral Green's function
% 
%             [SGFxx, SGFxy, SGFxz, SGFyx, SGFyy, SGFyz, SGFzx, SGFzy, SGFzz] = SGF(Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz, c);

       
            const = 1/dy;
            
            D_int = const .* SGFxx .* besselj(0, (kym .* w)/2); %.* (1 - exp(-1j .* kz0 .* 2 .* hd));
            
            D_inf = sum(D_int, 1);
            
    
            
            const2 = (-1/dx);
            
            yin_int = const2 .* sinc((kxm(1,:) .* deld)/2/pi).^2./(D_inf);
            
           z0 = 400;
            
            yin(p) = sum(yin_int);
             
            zin(p) = 1./yin(p);
            
           Gamma(p) = (zin(p) - z0) / (zin(p) + z0);
            
            
            del = 0.01 .* k0;
%             a = -50.*k0 - 1j.*del;
%             b = 50.*k0 + 1j.*del;
            
                 
%             for m = 0:Nx+2*Dummy-1
%                 int =  @(kxy) sinc((kxy .* deld)/2/pi).^2./(D_inf_func(kxy, th(j), ph, w, k0, dy, hd)) .* exp(-1j .* kxy .* abs(m) .* dx);
%                 if m == 0
%                    y = integral(int, -50*k0-1j*del,50*k0+1j*del, 'Waypoints', [(-1-1j).*del, (1+1j).*del]);
%                  else
%                    y = integral(int, -del-1j*20*k0,k0+del-1j*20*k0, 'Waypoints', [-del-del*1j,del+del*1j, k0+1j*del, k0+del]);
%                 end
%                  yin_int_mutual(m+1, 1) = -y/(2 * pi);
%             end
            
             for m = 0:Nx+2*Dummy-1
                int =  @(kxy) sinc((kxy .* deld)/2/pi).^2./(D_inf_func_er(kxy, th(j), ph, w, k0, dy, hd, er)) .* exp(-1j .* kxy .* abs(m) .* dx);
                if m == 0
                   y = integral(int, -50*k0-1j*del,50*k0+1j*del, 'Waypoints', [(-1-1j).*del, (1+1j).*del]);
                 else
                   y = integral(int, -del-1j*20*k0,2.2*k0+del-1j*20*k0, 'Waypoints', [-del-del*1j,del+del*1j, er*k0+1j*del, er*k0+del]);
                end
                 yin_int_mutual(m+1, 1) = -y/(2 * pi);
            end
            %z1 = 200;
            Y = toeplitz(real(yin_int_mutual)) + 1j * toeplitz(imag(yin_int_mutual));
            Z = inv(Y);
            idel = (Z + eye(Nx+2*Dummy)*z0)\v;
            
            Zn(:, p) = v./idel - z0; 
            
            Gamma_(:, p) = (Zn(:, p) - z0)./(Zn(:, p) + z0);
            
end
        [th1, ph1] = meshgrid(-pi/2-eps:dth:pi/2-eps, eps:dph:2*pi-eps);
        k0_ = (2 .* pi .* f0)/c0;
        
        kx = k0_ .* sin(th1) .* cos(ph1);
        ky = k0_ .* sin(th1) .* sin(ph1);
        
        
        
        krho_ = sqrt(kx.^2 + ky.^2);
        ksub_ = k0_ .* sqrt(er);
        
        kz0_ = (-1j)*sqrt(-(k0_.^2 - krho_.^2));
        kzs_ = (-1j)*sqrt(-(ksub_.^2 - krho_.^2));

            % For TM
            z0_TM = zeta * kz0_ ./ k0_;
            zs_TM = zeta_s * kzs_ ./ ksub_;

            % For TE

            z0_TE = zeta * k0_ ./ kz0_;
            zs_TE = zeta_s * ksub_ ./ kzs_;

            [zup_TE, zdn_TE, zup_TM, zdn_TM] = zupdn(z0_TE, zs_TE, z0_TM, zs_TM, hd, kzs_);


            [vtm, vte, itm, ite] = txline(zup_TE, zdn_TE, zup_TM, zdn_TM, hd, z, kz0_, kzs_, z0_TE, zs_TE, z0_TM, zs_TM);


            [SGFxx, SGFyx, SGFzx] = Green(vtm, vte, itm, ite, kx, ky, krho_.^2, zeta, k0_);
            
        
        AF_x_1 = zeros(size(kx));
        AF_x_2 = zeros(size(kx));
        AF_y = zeros(size(ky));
        AF_x = zeros(size(kx));
        
        for nx = 0:Nx-1
            
            AF_x_1 = AF_x_1 + (-v0 .* exp(1j .* (kx - kx0).*nx.*dx)); 
        end
        
        for nx = 0:Nx+1+Dummy
            AF_x_2 = AF_x_2 + z0 .* idel(nx + 1) .* exp(1j .* kx .* (nx - 2) .* dx); 
        end
        
       % AF_x = AF_x_1 + AF_x_2;
        
        for ny = 1:Nx
          AF_x = AF_x - exp(1j * (kx - kx0) * (ny-1) * dx);
          AF_y = AF_y + exp(1j * (ky - ky0) * (ny-1) * dy);
        end
      
      Sinc_term = sinc(kx0.*deld/2/pi)./(D_inf_func_er(kx0, th(j), ph, w, k0_, dy, hd, er));
      
      I_spectrum = Sinc_term .* besselj(0, (ky0*w/2)) .* AF_x .* AF_y;
      DD = 1./I_spectrum;
      
      Sinc_term_j(j) = Sinc_term;
      
      if j == 126
        figure(4);
      
         %plot(th1(1, :).*180/pi, abs(I_spectrum(1, :)), 'LineWidth', 3);
        hold on;
        %plot(th1(1, :).*180/pi, abs(DD(1, :)), 'LineWidth', 3);
        hold on;
        plot(th1(1, :).*180/pi, abs(AF_x(1, :)), 'LineWidth', 3);
      end
    
      
      
      r_obs = 1000 .* lambda0;
      
      Kons = 1j * kz0_ * exp(-1j * k0_ * r_obs)/(2 * pi * r_obs);
      
      ex = Kons .* SGFxx .* I_spectrum .* exp(-1j.*kzs_.*abs(z - hd));
      ey = Kons .* SGFyx .* I_spectrum .* exp(-1j.*kzs_.*abs(z - hd));
      ez = Kons .* SGFzx .* I_spectrum .* exp(-1j.*kzs_.*abs(z - hd));
      
      e_abs = sqrt(ex.^2 + ey.^2 + ez.^2);
      
      eth = ex .* cos(th1) .* cos(ph1) + ey .* cos(th1) .* sin(ph1) - ez .* sin(th1);
      eph = -sin(ph1) .* ex + cos(ph1) .* ey;
      
      V_far_th = eth * r_obs./(exp(-1j * k0_ * r_obs));
      V_far_ph = eph * r_obs./(exp(-1j * k0_ * r_obs));
    
        C_rad = 1/(2 * zeta);
        V_far_abs = sqrt(abs(V_far_th).^2 + abs(V_far_ph).^2);
        U = C_rad * (V_far_abs).^2;
    if th(j) >= 0
        P_rad_int = C_rad * (V_far_abs(:, 91:end)).^2 .* sin(th1(:, 91:end)) * dth * dph;
    else
        P_rad_int = C_rad * (V_far_abs(:, 1:90)).^2 .* sin(th1(:, 1:90)) * dth * dph;
    end
        P_rad = nansum(nansum(P_rad_int));
        
       if th(j) > 0 
        Dir(j) = 4 * pi * max(max(U(:, 90:end))) ./ P_rad;
       elseif th(j) == eps
           Dir(j) = 4 * pi * max(max(U(1, 91))) ./ P_rad;
       else
        Dir(j) = 4 * pi * max(max(U(:, 1:90))) ./ P_rad;
       end
           
      Dir_e = 4 * pi * U ./ P_rad;
        %if th(j) == eps || th(j) == pi/6 || th(j) == pi/4 || th(j) == pi/3 || th(j) == pi/2
        %if j == 31 || j == 31+12 || j == 31+7 || j == 31+20
        if j == 115 || j == 116 || j == 117 || j == 118 || j == 119 || j == 120 || j == 121 ...
                || j == 122 || j == 123 || j == 124 || j == 125 || j == 126 || j == 127 || ...
                j == 128 || j == 129
            hold on;
            figure(1);
            hold on;
            plot(th1(1, :).*180/pi, 10*log10(abs(Dir_e(1, :))), 'LineWidth', 3);
        end
        if j == 126
            hi = Dir_e;
            max(max(hi))
        end
%         
%         figure(1);
%         
%         plot(th1(1, :).*180/pi, db(e_abs(1, :)/max(e_abs(1, :))), 'LineWidth', 3)
%        
%         grid on;
%         hold on;
%         
%         plot(th1(91, :).*180/pi, db(e_abs(91, :)/max(e_abs(91, :))), 'LineWidth', 3)
        
       
      
% figure(1);
% 
% for i = Dummy+1:length(v) - Dummy
%     plot(f*10^(-9), real(Zn(i, :)), 'r', 'LineWidth', 1);
%     grid on;
%     hold on;
%     plot(f*10^(-9), imag(Zn(i, :)), '--r', 'LineWidth', 1);
%     hold on;
% end
% plot(f*10^(-9), real(zin), 'k', 'LineWidth', 3);
% hold on;
% plot(f*10^(-9), imag(zin), '--k', 'LineWidth', 3);
% 
% 
% xlabel('frequency(GHz)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Re Z_{a, dipole} , Im Z_{a, dipole} (\Omega) ', 'FontSize', 12, 'FontWeight', 'bold');
% title('Active input impedance of all finite array (Red), Black Curve(Infinite Array)', 'FontSize', 12, 'FontWeight', 'bold');
% %legend({'Real Z_{in}(0, 0)','Imag Z_{in}(0, 0)', 'Real Z_{in}(45, 0)','Imag Z_{in}(45, 0)'},'Location','northeast', 'FontSize', 12, 'FontWeight', 'bold');
% %legend({'Real Z_{in}','Imag Z_{in}'},'Location','northeast', 'FontSize', 12, 'FontWeight', 'bold');
% %print(['Final_Q1_Z_dipole_th_', num2str(th(j)*180/pi), '_ph_', num2str(ph*180/pi), '_', num2str(Nx), '_er_', num2str(er)], '-depsc');
% 
% 
% %% Gamma 
% 
% 
% 
% figure(2);
% 
% 
% for i = Dummy+1:length(v) - Dummy
%     plot(f*10^(-9), db(abs(Gamma_(i, :))), 'r', 'LineWidth', 1);
%     hold on;
% end
% plot(f*10^(-9), db(abs(Gamma)), 'k', 'LineWidth', 3);
% hold on;
% grid on;
% %ylim([-30 0]);
% xlabel('frequency(GHz)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Active |\Gamma| dB', 'FontSize', 12, 'FontWeight', 'bold');
% title('Active reflection coefficients of all elements finite array (Red), Black Curve(Infinite Array)', 'FontSize', 12, 'FontWeight', 'bold');
% %legend({'Real Z_{in}(0, 0)','Imag Z_{in}(0, 0)', 'Real Z_{in}(45, 0)','Imag Z_{in}(45, 0)'},'Location','northeast', 'FontSize', 12, 'FontWeight', 'bold');
% %legend({'Real Z_{in}','Imag Z_{in}'},'Location','northeast', 'FontSize', 12, 'FontWeight', 'bold');
% %print(['Final_Q1_Gamma_dipole_th_', num2str(th(j)*180/pi), '_ph_', num2str(ph*180/pi), '_', num2str(Nx), '_er_', num2str(er)], '-depsc');



% 
% figure(1);
% 
% plot(f*10^(-9), 20*log10(abs(Gamma)), 'LineWidth', 3);
% grid on;
% hold on;

end



% plot(th.*180/pi, abs(Dir), 'LineWidth', 3, 'Color', [0.6350, 0.0780, 0.1840]);
% grid on;
figure(1);
plot(th.*180/pi, 10*log10(abs(Dir)), 'LineWidth', 3);
hold on;
grid on;

% grid on;
xlabel('Scan angle and Observation Angle (\theta (deg))', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Directivity and Envelope (dB) at \phi = 0', 'FontSize', 12, 'FontWeight', 'bold');
title('Directivity and Envelope', 'FontSize', 12, 'FontWeight', 'bold');
%legend({'Real Z_{in}(0, 0)','Imag Z_{in}(0, 0)', 'Real Z_{in}(45, 0)','Imag Z_{in}(45, 0)'},'Location','northeast', 'FontSize', 12, 'FontWeight', 'bold');
%legend({'Real Z_{in}','Imag Z_{in}'},'Location','northeast', 'FontSize', 12, 'FontWeight', 'bold');
%ylim([-40 0]);
% %ylim([-2000 3000]);
% print('A5_Q1_dipole_br', '-depsc');
% %xlim([-180 180]);

  figure(3);
      plot(th(1, :).*180/pi, Sinc_term_j(1, :), 'LineWidth', 3);
      hold on;
      
      


% xlabel('Frequency(GHz)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Active \Gamma (dB)', 'FontSize', 12, 'FontWeight', 'bold');
% title(['Active reflection coefficient at z0 = ', num2str(z0), '\Omega'], 'FontSize', 12, 'FontWeight', 'bold');

%legend({'\Gamma (0, 0)', '\Gamma (45, 0) E Plane', '\Gamma (45, 90) H Plane'},'Location','northeast', 'FontSize', 12, 'FontWeight', 'bold');

