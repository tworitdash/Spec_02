
%% PDF of 2 hypothesis; One without the target and other with the target

x = -10:0.001:10; % window for the normal distributon without the target 
nu = 4;           % degrees of freedom for the chi squared distribution 
                  %in the case of target present
                  
x_c = 0:0.001:30; % window for chi squared distribution 
                  %(because this distribution always starts from 0)

pd = makedist('Normal','mu',-1,'sigma',1); %Normal distribution object

pdf_norm = pdf(pd, x);                 %PDF of Normal distribution                         
pdf_chi = chi2pdf(x_c, nu);            %PDF of Chi Squared distribution

figure(1);
plot(x, pdf_norm, 'LineWidth', 2);      %Plot of Normal distribution
grid on;
hold on;

plot(x_c, pdf_chi, 'LineWidth', 2);    %Plot of Chi squared distribution




xlabel('x', 'FontSize', 12, 'FontWeight', 'bold'); 
ylabel('PDF', 'FontSize', 12, 'FontWeight', 'bold');
title('PDF of received signal', 'FontSize', 12, 'FontWeight', 'bold');

legend({'No target','With 1 target'},'Location','northeast', 'FontSize', 12, 'FontWeight', 'bold');

%% ROC Plot

Pfa = eps:0.01:1-eps;    %Range of false alarm probabilities 
                         % for which detection  
                         %probablities need to be evaluated
mu = -1;                 % mean of normal distribution
sigma = 1;               %standard deviation of normal distribution
gamma = icdf('norm', 1-Pfa, mu, sigma); 

%Inverse CDF (Area under the PDF) to calculate the threshold


Pmd = chi2cdf(gamma, 4);     % Probablility of missed detection

Pd = 1 - Pmd;                % Probability of detection

figure(2);

plot(Pfa, Pd, 'LineWidth', 2) %ROC plot (Pd vs Pfa)
grid on;
xlabel('P_{fa}', 'FontSize', 12, 'FontWeight', 'bold'); 
ylabel('P_d', 'FontSize', 12, 'FontWeight', 'bold');
title('ROC', 'FontSize', 12, 'FontWeight', 'bold');