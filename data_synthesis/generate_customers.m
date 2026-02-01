function generate_customers(nCustomers, nSons, varargin)
% generate_customers
% Generates nCustomers synthetic customer energy logs into:
%   SM_Correlation/data/virtual_customers
%
% nSons customers will be marked as sons (is_son=1), rest as orphans (0).
%
% Each customer:
%   E ~ Uniform[5500,8000]
%   sigma = E/4
%   N quarter-hourly samples 
%
% Uses:
%   gaussian_NE.m
%   generate_from_indexed_NE.m

p = inputParser;
p.addParameter('rngSeed', [], @(x) isempty(x) || (isnumeric(x)&&isscalar(x)));
p.addParameter('periods', [], @(x) isempty(x) || (isnumeric(x)&&isscalar(x)));
p.parse(varargin{:});
opt = p.Results;

if nargin < 2
    nSons = nCustomers; % default: all are sons
end
assert(nSons >= 0 && nSons <= nCustomers, 'nSons must be in [0, nCustomers].');

if isempty(opt.rngSeed)
    rng('shuffle');   % new seed for every run
else
    rng(opt.rngSeed); % reproduce seed
end

startTime = datetime(2022,11,23,1,0,0);

% --- Resolve project root ---
thisFile = mfilename('fullpath');
dataSynthesisDir = fileparts(thisFile);
projectRoot = fileparts(dataSynthesisDir);

% --- Output directory: data/virtual_customers ---
outDir = fullfile(projectRoot, "data", "virtual_customers");
if ~exist(outDir,'dir')
    mkdir(outDir);
end

% --- Decide ground truth sons/orphans ---
is_son = false(nCustomers,1);
is_son(1:nSons) = true;                 % simplest: first nSons are sons
is_son = is_son(randperm(nCustomers));  % shuffle so it's not trivial

% --- Metadata table ---
meta = table('Size',[nCustomers 5], ...
    'VariableTypes',{'double','double','double','string','logical'}, ...
    'VariableNames',{'CustomerID','E_mean','sigma','filename','is_son'});

for k = 1:nCustomers
    % --- Random E and sigma ---
    E = 5500 + (8000-5500)*rand();
    sigma = E/4;

    % --- Gaussian indexed N(E) ---
    [E_i, n_i] = gaussian_NE(E, sigma^2);

    % --- Output filename ---
    fname = sprintf("customer_%03d.xlsx", k);
    fpath = fullfile(outDir, fname);


    if isempty(opt.periods)
        nPeriods = 1500;
    else
        nPeriods = opt.periods;
    end

    % --- Generate Excel log ---
    generate_from_indexed_NE(E_i, n_i, startTime, fpath, nPeriods);

    % --- Save metadata ---
    meta.CustomerID(k) = k;
    meta.E_mean(k) = E;
    meta.sigma(k) = sigma;
    meta.filename(k) = fname;
    meta.is_son(k) = is_son(k);

    fprintf("Generated %s | E=%.1f | sigma=%.1f | is_son=%d\n", fname, E, sigma, meta.is_son(k));
end

% --- Save metadata ---
writetable(meta, fullfile(outDir, "customers_metadata.xlsx"));

fprintf("----------------------------------------------------\n");
fprintf("%d virtual premises generated in:\n%s\n", nCustomers, outDir);
fprintf("Ground truth saved to customers_metadata.xlsx (column: is_son)\n");
end
