function out = preprocess_synthetic_data(varargin)
% preprocess_synthetic_data
% Creates a corrupted copy of the synthetic dataset by injecting:
%   - measurement noise
%   - missing data
%   - anomalies (spikes / level shifts / dropout bursts)
%
% Default input dataset:
%   data/virtual_customers/customer_*.xlsx
%   data/virtual_customers/customers_metadata.xlsx (copied if exists)
%   data/transformer_profile.xlsx
%
% Output dataset (defaults):
%   data/virtual_customers_pp/customer_*.xlsx
%   data/virtual_customers_pp/customers_metadata.xlsx
%   data/transformer_profile_pp.xlsx
%
% Usage examples:
%   preprocess_synthetic_data(); % default mild corruption
%   preprocess_synthetic_data('noiseStdCustomer',5,'noiseStdTransformer',10);
%   preprocess_synthetic_data('missingProbCustomer',0.02,'missingProbTransformer',0.01);
%   preprocess_synthetic_data('spikeProb',0.002,'spikeAmp',200);
%   preprocess_synthetic_data('seed',7);
%
% Notes:
% - Corruptions are applied to EC column only.
% - Missing values are written as NaN in the Excel files.
% - UKF should handle NaNs only if your loader/interpolator supports it.
%   (If needed, add an interpolation step in your loader.)

% ---------------- Parse inputs ----------------
p = inputParser;

p.addParameter('projectRoot', "", @(x)isstring(x)||ischar(x));

% input locations (relative to projectRoot)
p.addParameter('customersSubdirIn',  fullfile("data","virtual_customers"), @(x)isstring(x)||ischar(x));
p.addParameter('transformerFileIn',  fullfile("data","transformer_profile.xlsx"), @(x)isstring(x)||ischar(x));

% output locations
p.addParameter('customersSubdirOut', fullfile("data","virtual_customers_pp"), @(x)isstring(x)||ischar(x));
p.addParameter('transformerFileOut', fullfile("data","transformer_profile_pp.xlsx"), @(x)isstring(x)||ischar(x));

% RNG control (optional)
p.addParameter('seed', [], @(x) isempty(x) || (isnumeric(x)&&isscalar(x)));

% Measurement noise (Gaussian)
p.addParameter('noiseStdCustomer',    0, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('noiseStdTransformer', 0, @(x)isnumeric(x)&&isscalar(x)&&x>=0);

% Missingness (MCAR)
p.addParameter('missingProbCustomer',    0, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
p.addParameter('missingProbTransformer', 0, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);

% Anomalies
p.addParameter('enableSpikes', true, @(x)islogical(x)||ismember(x,[0 1]));
p.addParameter('spikeProb', 0, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1); % per-sample probability
p.addParameter('spikeAmp',  0, @(x)isnumeric(x)&&isscalar(x)&&x>=0);       % absolute added magnitude (same units as EC)

p.addParameter('enableLevelShift', false, @(x)islogical(x)||ismember(x,[0 1]));
p.addParameter('levelShiftProb', 0, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1); % per-series probability
p.addParameter('levelShiftAmp',  0, @(x)isnumeric(x)&&isscalar(x));             % added offset magnitude

p.addParameter('enableDropoutBursts', false, @(x)islogical(x)||ismember(x,[0 1]));
p.addParameter('dropoutBurstProb', 0, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1); % per-series probability
p.addParameter('dropoutBurstLen',  0, @(x)isnumeric(x)&&isscalar(x)&&x>=0);       % length in samples

p.parse(varargin{:});
opt = p.Results;

% ---------------- Resolve project root ----------------
if strlength(string(opt.projectRoot))==0
    thisFile = mfilename('fullpath');
    thisDir  = fileparts(thisFile);
    projectRoot = string(fileparts(thisDir));
else
    projectRoot = string(opt.projectRoot);
end

customersDirIn   = fullfile(projectRoot, string(opt.customersSubdirIn));
customersDirOut  = fullfile(projectRoot, string(opt.customersSubdirOut));
transformerIn    = fullfile(projectRoot, string(opt.transformerFileIn));
transformerOut   = fullfile(projectRoot, string(opt.transformerFileOut));

if ~isfolder(customersDirIn)
    error("Input customers folder not found: %s", customersDirIn);
end
if ~isfile(transformerIn)
    error("Input transformer file not found: %s", transformerIn);
end

% RNG
if isempty(opt.seed)
    rng('shuffle');
    seedUsed = NaN;
else
    rng(opt.seed);
    seedUsed = opt.seed;
end

% ---------- Clean previous PP data ----------
if exist(customersDirOut,'dir')
    fprintf("Cleaning previous preprocessed data...\n");
    delete(fullfile(customersDirOut, '*'));   % delete all files
else
    mkdir(customersDirOut);
end

% Delete old transformer pp file if exists
if exist(transformerOut,'file')
    delete(transformerOut);
end

fprintf("\n=== preprocess_synthetic_data ===\n");
fprintf("Input customers   : %s\n", customersDirIn);
fprintf("Output customers  : %s\n", customersDirOut);
fprintf("Input transformer : %s\n", transformerIn);
fprintf("Output transformer: %s\n", transformerOut);
if isnan(seedUsed), fprintf("RNG seed          : shuffle\n"); else, fprintf("RNG seed          : %d\n", seedUsed); end

% ---------------- Copy & corrupt customer files ----------------
files = dir(fullfile(customersDirIn, "customer_*.xlsx"));
if isempty(files)
    error("No customer_*.xlsx found in: %s", customersDirIn);
end
[~,ix] = sort({files.name}); files = files(ix);
Nfiles = numel(files);

stats = struct();
stats.nCustomers = Nfiles;
stats.customers = repmat(struct( ...
    'file',"", 'nSamples',0, 'nMissing',0, 'nSpikes',0, 'levelShift',0, 'dropoutBurst',0), Nfiles, 1);

for i = 1:Nfiles
    inPath  = fullfile(files(i).folder, files(i).name);
    outPath = fullfile(customersDirOut, files(i).name);

    T = readtable(inPath);
    assert(all(ismember({'STARTOFINTERVAL_TIMES','EC'}, T.Properties.VariableNames)), ...
        "Customer file missing required columns: %s", files(i).name);

    x = double(T.EC(:));
    Tn = numel(x);

    % Noise
    if opt.noiseStdCustomer > 0
        x = x + opt.noiseStdCustomer * randn(Tn,1);
    end

    % Spikes (per-sample)
    nSp = 0;
    if opt.enableSpikes && opt.spikeProb > 0 && opt.spikeAmp > 0
        isSpike = rand(Tn,1) < opt.spikeProb;
        sgn = sign(randn(Tn,1)); sgn(sgn==0)=1;
        x(isSpike) = x(isSpike) + sgn(isSpike)*opt.spikeAmp;
        nSp = sum(isSpike);
    end

    % Level shift (per-series)
    lvl = 0;
    if opt.enableLevelShift && opt.levelShiftProb > 0
        if rand() < opt.levelShiftProb
            lvl = opt.levelShiftAmp;
            x = x + lvl;
        end
    end

    % Dropout burst (per-series): consecutive NaNs
    doBurst = 0;
    if opt.enableDropoutBursts && opt.dropoutBurstProb > 0 && opt.dropoutBurstLen > 0
        if rand() < opt.dropoutBurstProb
            doBurst = 1;
            L = min(Tn, round(opt.dropoutBurstLen));
            startIdx = randi(max(1, Tn-L+1));
            x(startIdx:startIdx+L-1) = NaN;
        end
    end

    % Missingness MCAR (per-sample)
    nMiss = 0;
    if opt.missingProbCustomer > 0
        isMiss = rand(Tn,1) < opt.missingProbCustomer;
        x(isMiss) = NaN;
        nMiss = sum(isMiss) + sum(isnan(x)) - sum(isnan(T.EC(:))); % include burst NaNs too
    else
        nMiss = sum(isnan(x));
    end

    % Physical constraint (optional)
    x = max(x, 0); % keeps NaN as NaN

    T.EC = x;
    writetable(T, outPath, 'FileType','spreadsheet');

    stats.customers(i).file = string(files(i).name);
    stats.customers(i).nSamples = Tn;
    stats.customers(i).nMissing = sum(isnan(T.EC));
    stats.customers(i).nSpikes = nSp;
    stats.customers(i).levelShift = lvl;
    stats.customers(i).dropoutBurst = doBurst;
end

% Copy metadata if exists
metaIn  = fullfile(customersDirIn,  "customers_metadata.xlsx");
metaOut = fullfile(customersDirOut, "customers_metadata.xlsx");
if isfile(metaIn)
    copyfile(metaIn, metaOut);
end

% ---------------- Corrupt transformer profile ----------------
Ttr = readtable(transformerIn);
assert(all(ismember({'STARTOFINTERVAL_TIMES','EC'}, Ttr.Properties.VariableNames)), ...
    "Transformer file must include STARTOFINTERVAL_TIMES and EC");

y = double(Ttr.EC(:));
Ty = numel(y);

% noise
if opt.noiseStdTransformer > 0
    y = y + opt.noiseStdTransformer * randn(Ty,1);
end

% spikes
if opt.enableSpikes && opt.spikeProb > 0 && opt.spikeAmp > 0
    isSpikeY = rand(Ty,1) < opt.spikeProb;
    sgn = sign(randn(Ty,1)); sgn(sgn==0)=1;
    y(isSpikeY) = y(isSpikeY) + sgn(isSpikeY)*opt.spikeAmp;
end

% missingness
if opt.missingProbTransformer > 0
    isMissY = rand(Ty,1) < opt.missingProbTransformer;
    y(isMissY) = NaN;
end

y = max(y, 0); % keeps NaN
Ttr.EC = y;
writetable(Ttr, transformerOut, 'FileType','spreadsheet');

% ---------------- Output struct ----------------
out = struct();
out.projectRoot = projectRoot;
out.customersDirIn = customersDirIn;
out.customersDirOut = customersDirOut;
out.transformerIn = transformerIn;
out.transformerOut = transformerOut;
out.seedUsed = seedUsed;
out.options = opt;
out.stats = stats;

fprintf("Done. Corrupted dataset ready.\n\n");
end
