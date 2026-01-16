function generate_100_customers()
% generate_100_customers
%   Generates 100 synthetic customer energy logs and saves them into:
%       SM_Correlation/data/virtual_customers
%
%   Each customer:
%       E ~ Uniform[5500,8000]
%       sigma = E/4
%       18000 quarter-hourly samples (~6 months)
%
%   Uses:
%       gaussian_NE.m
%       generate_from_indexed_NE.m

    rng(1);                        % Reproducibility
    nCustomers = 100;
    startTime = datetime(2022,11,23,1,0,0);

    % --- Resolve project root ---
    thisFile = mfilename('fullpath');
    dataSynthesisDir = fileparts(thisFile);
    projectRoot = fileparts(dataSynthesisDir);

    % --- Output directory: data/virtual_customers ---
    dataDir = fullfile(projectRoot, "data");
    outDir  = fullfile(dataDir, "virtual_customers");

    if ~exist(outDir,'dir')
        mkdir(outDir);
    end

    % --- Metadata table ---
    meta = table('Size',[nCustomers 4], ...
        'VariableTypes',{'double','double','double','string'}, ...
        'VariableNames',{'CustomerID','E_mean','sigma','filename'});

    for k = 1:nCustomers

        % --- Random E and sigma ---
        E = 5500 + (8000-5500)*rand();
        sigma = E/4;

        % --- Gaussian indexed N(E) ---
        [E_i, n_i] = gaussian_NE(E, sigma^2);

        % --- Output filename ---
        fname = sprintf("customer_%03d.xlsx",k);
        fpath = fullfile(outDir, fname);

        % --- Generate Excel log ---
        generate_from_indexed_NE(E_i, n_i, startTime, fpath);

        % --- Save ground truth ---
        meta.CustomerID(k) = k;
        meta.E_mean(k) = E;
        meta.sigma(k) = sigma;
        meta.filename(k) = fname;

        fprintf("Generated %s   E = %.1f   sigma = %.1f\n", fname, E, sigma);
    end

    % --- Save metadata ---
    writetable(meta, fullfile(outDir, "customers_metadata.xlsx"));

    disp("----------------------------------------------------");
    disp("100 virtual premises generated in:");
    disp(outDir);
    disp("Ground truth saved to customers_metadata.xlsx");
end
