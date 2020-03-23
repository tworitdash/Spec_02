close all;
clear all;


f = 2e9:0.1e9:10e9;
%f = 2e9;
zeta = 120 .* pi;
e0 = 8.85418782e-12;

Gamma_TE = zeros(size(f));
Gamma_TM = zeros(size(f));
S11_TE = zeros(size(f));
S12_TE = zeros(size(f));
S11_TM = zeros(size(f));
S12_TM = zeros(size(f));

c = 3e8;

for i = 1:size(f, 2)

    lambda = c./f(end);
    dx = 0.2 .*  lambda;
    dy = 0.2 .* lambda;
    

    w = 0.01 .* lambda;
    dz = 0.01 .* lambda;
    
    m = [-30:-1,1:30];
    
    Const = (2 .* f(i) .* dy .* e0);
    
    B_ = Const .* abs(sinc(pi .* m .* w ./ dy./ pi)).^2 ./ abs(m);
    
    B_inf = B_ .* 1j .* tan((-1j .* pi .* abs(m) .* dz)./dy);
    B_semiInf = B_ .* (1/2 + 1j./2 .* tan((-1j .* pi .* abs(m) .* dz)./dy));
    
    Binf = sum(B_inf);
    BsemiInf = sum(B_semiInf);  
    
    th = pi/3;
    
    k0 = 2.*pi.*f(i)./c;
    
    kz = k0 .* cos(th);
    
    Z0_TM = zeta .* kz./k0;
    Z0_TE = zeta .* k0./kz;
    
    ZTE_in = (-1j./Binf) .* (1./(1 - (sin(th).^2./2)));
    
    ZTM_in = (-1j./Binf);
    
    ZTE_out = (-1j./BsemiInf) .* (1./(1 - (sin(th).^2./2)));
    ZTM_out = (-1j./BsemiInf);
    
    ABCD_tr_TE = [cos(kz.*dz) 1j.*Z0_TE.*sin(kz.*dz); 1j.*(1./Z0_TE).*sin(kz.*dz) cos(kz.*dz)];
    ABCD_out_TE = [1 0; 1./ZTE_out 1];
    ABCD_in_TE = [1 0; 1./ZTE_in 1];
    
    ABCD_tr_TM = [cos(kz.*dz) 1j.*Z0_TM.*sin(kz.*dz); 1j.*(1./Z0_TM).*sin(kz.*dz) cos(kz.*dz)];
    ABCD_out_TM = [1 0; 1./ZTM_out 1];
    ABCD_in_TM = [1 0; 1./ZTM_in 1];
    
    N = 5; %number of layers
    
    ABCD_TE = [1 0; 0 1];
    ABCD_TM = [1 0; 0 1];
    
    for k = 1:N-2
        % for TE
        ABCD_TE = ABCD_TE * ABCD_in_TE * ABCD_tr_TE;
        
        % for TM
        ABCD_TM = ABCD_TM * ABCD_in_TM * ABCD_tr_TM;
    end
    
    ABCD_tot_TE = ABCD_out_TE * ABCD_tr_TE * ABCD_TE * ABCD_tr_TE * ABCD_out_TE;
    
    ABCD_tot_TM = ABCD_out_TM * ABCD_tr_TM * ABCD_TM * ABCD_tr_TM * ABCD_out_TM;
    
    A_TE = ABCD_tot_TE(1, 1);
    A_TM = ABCD_tot_TM(1, 1);
    B_TE = ABCD_tot_TE(1, 2);
    B_TM = ABCD_tot_TM(1, 2);
    C_TE = ABCD_tot_TE(2, 1);
    C_TM = ABCD_tot_TM(2, 1);
    D_TE = ABCD_tot_TE(2, 2);
    D_TM = ABCD_tot_TM(2, 2);
    
    S11_TE(i) = (A_TE + B_TE./Z0_TE - C_TE.*Z0_TE - D_TE)./(A_TE + B_TE./Z0_TE + C_TE.*Z0_TE + D_TE);
    S12_TE(i) = 2.*(A_TE .* D_TE - B_TE .* C_TE)./(A_TE + B_TE./Z0_TE + C_TE.*Z0_TE + D_TE);
    
    S11_TM(i) = (A_TM + B_TM./Z0_TM - C_TM.*Z0_TM - D_TM)./(A_TM + B_TM./Z0_TM + C_TM.*Z0_TM + D_TM);
    S12_TM(i) = 2.*(A_TM .* D_TM - B_TM .* C_TM)./(A_TM + B_TM./Z0_TM + C_TM.*Z0_TM + D_TM);
    
    
%     
%     ZTE_p = (ZTE .* Z0_TE) ./ (ZTE + Z0_TE);
%     ZTM_p = (ZTM .* Z0_TM) ./ (ZTM + Z0_TM);
%     
%     Gamma_TE(i) = (ZTE_p - Z0_TE)./(ZTE_p + Z0_TE);
%     Gamma_TM(i) = (ZTM_p - Z0_TM)./(ZTM_p + Z0_TM);
%     
end

plot(f .* 10^(-9), abs(S11_TE).^2, 'LineWidth', 2);
hold on;
plot(f .* 10^(-9), abs(S11_TM).^2, 'LineWidth', 2);

grid on;

plot(f .* 10^(-9), abs(S12_TE).^2, '--', 'LineWidth', 2);
hold on;
plot(f .* 10^(-9), abs(S12_TM).^2, '--', 'LineWidth', 2);

xlabel('frequency(GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('|S_{11}|^2, |S_{12}|^2', 'FontSize', 12, 'FontWeight', 'bold');
title(['Reflection and Transmission Coefficients at \theta = ', num2str(th*180/pi)], 'FontSize', 12, 'FontWeight', 'bold');
legend({'|S_{11}|^2 TE', '|S_{11}|^2 TM', '|S_{12}|^2 TE', '|S_{12}|^2 TM'},'Location','south', 'FontSize', 12, 'FontWeight', 'bold');

