%% ROC Plot

Pfa = eps:0.01:1-eps;               %Range of false alarm probabilities for which detection  
                                    %probablities need to be evaluated
mu = -1;                            % mean of normal distribution
sigma = 1;                          %standard deviation of normal distribution
gamma = icdf('norm', 1-Pfa, mu, sigma); %Inverse CDF (Area under the PDF) to calculate the threshold


Pmd = chi2cdf(gamma, 4);            % Probablility of missed detection

Pd = 1 - Pmd;                       % Probability of detection

plot(Pfa, Pd, 'LineWidth', 2)       %ROC plot (Pd vs Pfa)
grid on;
xlabel('P_{fa}', 'FontSize', 12, 'FontWeight', 'bold'); 
ylabel('P_d', 'FontSize', 12, 'FontWeight', 'bold');
title('ROC', 'FontSize', 12, 'FontWeight', 'bold');
