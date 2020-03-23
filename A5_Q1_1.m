close all;
clear all;


f = 2e9:0.1e9:10e9;
zeta = 120 .* pi;
e0 = 8.85418782e-12;

Gamma_TE = zeros(size(f));
Gamma_TM = zeros(size(f));


c = 3e8;

for i = 1:size(f, 2)

    lambda = c./f(end);
    dx = 0.2 .*  lambda;
    dy = 0.2 .* lambda;

    w = 0.01 .* lambda;
    
    m = [-10:-1,1:10];
    
    Const = (2 .* f(i) .* dy .* e0);
    
    B_ = Const .* abs(sinc(pi .* m .* w ./ dy./ pi)).^2 ./ abs(m);
    
    B = sum(B_);
    
    th = pi/3;
    
    k0 = 2.*pi.*f(i)./c;
    
    kz = k0 .* cos(th);
    
    Z0_TM = zeta .* kz./k0;
    Z0_TE = zeta .* k0./kz;
    
    ZTE = (-1j./B) .* (1./(1 - (sin(th).^2./2)));
    
    ZTM = (-1j./B);
    
    ZTE_p = (ZTE .* Z0_TE) ./ (ZTE + Z0_TE);
    ZTM_p = (ZTM .* Z0_TM) ./ (ZTM + Z0_TM);
    
    Gamma_TE(i) = (ZTE_p - Z0_TE)./(ZTE_p + Z0_TE);
    Gamma_TM(i) = (ZTM_p - Z0_TM)./(ZTM_p + Z0_TM);
    
end

plot(f .* 10^(-9), abs(Gamma_TE).^2, 'LineWidth', 2);
hold on;
plot(f .* 10^(-9), abs(Gamma_TM).^2, 'LineWidth', 2);

grid on;

plot(f .* 10^(-9), 1 - abs(Gamma_TE).^2, '--', 'LineWidth', 2);
hold on;
plot(f .* 10^(-9), 1 - abs(Gamma_TM).^2, '--', 'LineWidth', 2);

xlabel('frequency(GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('|S_{11}|^2, |S_{12}|^2', 'FontSize', 12, 'FontWeight', 'bold');
title(['Reflection and Transmission Coefficients at \theta =', num2str(th)], 'FontSize', 12, 'FontWeight', 'bold');
legend({'|S_{11}|^2 TE', '|S_{11}|^2 TM', '|S_{12}|^2 TE', '|S_{12}|^2 TM'},'Location','west', 'FontSize', 12, 'FontWeight', 'bold');

print('A5_Q1_Gamma_T', '-depsc');
