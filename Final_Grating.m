%% Grating lobe diagrams for both the configuration:

dx = 15e-3;
dy = 15e-3;

f = 10.7e9;
ph = eps;

lambda = 3e8./f;
k0 = 2 * pi ./ lambda;
zeta = 120*pi;

mx1 = -1:1:1;
my1 = -1:1:1;

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
        
        plot(x1 + xp, y1 + yp, 'LineWidth', 2);
        hold on;
        quiver(x1, y1, kx0, ky0, 'LineWidth', 2);
        
    end
end
grid on;
 xlabel('k_{xm}', 'FontSize', 12, 'FontWeight', 'bold');
 ylabel('k_{ym}', 'FontSize', 12, 'FontWeight', 'bold');
 title('Grating lobe diagram when dx = dy =0.53\lambda at \phi = 0 (E Plane)', 'FontSize', 12, 'FontWeight', 'bold');
%legend({'Real Z_{in}','Imag Z_{in}'},'Location','northeast', 'FontSize', 12, 'FontWeight', 'bold');
print('Final_Grating', '-dpng')
