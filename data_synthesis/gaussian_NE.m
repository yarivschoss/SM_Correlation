function [E_i, n_i, samples] = gaussian_NE(E_mean, sigma2, nSamples)
% gaussian_NE
%   Generates the Gaussian distribution N(E) and returns it as indexed series
%
% Outputs:
%   E_i     - energy axis (E(i))
%   n_i     - N(E(i)) Gaussian PDF values
%   samples - random samples drawn from N(E_mean, sigma2)
%
% Inputs:
%   E_mean  - mean energy E
%   sigma2  - variance
%   nSamples - number of random samples (default = 10000)

    if nargin < 3
        nSamples = 10000;
    end

    sigma = sqrt(sigma2);

    % Energy axis (indexed)
    E_i = linspace(E_mean - 4*sigma, E_mean + 4*sigma, 1000);

    % N(E): Gaussian probability density function
    n_i = (1/(sigma*sqrt(2*pi))) * exp(-0.5*((E_i - E_mean)/sigma).^2);

    % Random samples from N(E)
    samples = E_mean + sigma * randn(nSamples,1);

    % Plot using the indexed series
    figure('Color','w');
    histogram(samples, 60, 'Normalization','pdf');
    hold on;
    plot(E_i, n_i, 'r', 'LineWidth', 2);
    grid on;
    xlabel("Energy E(i)");
    ylabel("N(E(i))");
    title(sprintf("Gaussian Distribution  N(E = %.0f, \\sigma = %.0f)", E_mean, sigma));
    legend("Samples","N(E(i))");
end
