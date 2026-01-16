function out = synthesize_customer_end2end(E, startTime, outXlsx)
% synthesize_customer_end2end
%   1) Builds Gaussian N(E) with sigma = E/4 using gaussian_NE
%   2) Generates a continuous 18000-sample quarter-hourly Excel log using generate_from_indexed_NE
%   3) Plots the generated Excel log using plot_excel_energy_log
%
% Inputs:
%   E         - mean energy value
%   startTime - datetime scalar
%   outXlsx   - output Excel file name/path (e.g., "synthetic_customer.xlsx")
%
% Output:
%   out - struct with fields: E, sigma, sigma2, outXlsx

    if nargin < 3
        outXlsx = "synthetic_customer.xlsx";
    end

    sigma  = E/4;
    sigma2 = sigma^2;

    % 1) Gaussian indexed distribution: [E_i, n_i]
    [E_i, n_i] = gaussian_NE(E, sigma2);

    % 2) Generate continuous log + write Excel
    %    (Assumes your generate_from_indexed_NE signature is:
    %     generate_from_indexed_NE(E_i, n_i, startTime, outXlsx))
    generate_from_indexed_NE(E_i, n_i, startTime, outXlsx);

    % 3) Plot from Excel
    plot_excel_energy_log(outXlsx);

    % Return useful info
    out.E = E;
    out.sigma = sigma;
    out.sigma2 = sigma2;
    out.outXlsx = outXlsx;
end
