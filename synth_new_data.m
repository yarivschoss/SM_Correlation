function synth_new_data(varargin)
% synth_new_data
% Wrapper to clear + synthesize customers + build transformer profile.
%
% Usage examples:
%   synth_new_data('meters',120,'sons',100)
%   synth_new_data('meters',120,'sons',100,'periods',200)
%   synth_new_data('meters',120,'sons',100,'seed',42)
%   synth_new_data('meters',120,'sons',100,'periods',200,'seed',42)
%
% Params (Name-Value):
%   'meters'  : number of customers (default 120)
%   'sons'    : number of sons (default 100)
%   'periods' : number of periods (optional; if omitted -> not passed)
%   'seed'    : rngSeed (optional; if omitted -> not passed; generate_customers does shuffle)

% ---------- Parse Name-Value ----------
p = inputParser;
p.addParameter('meters', 120, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('sons',   100, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('periods', [], @(x) isempty(x) || (isnumeric(x)&&isscalar(x)&&x>0));
p.addParameter('seed',    [], @(x) isempty(x) || (isnumeric(x)&&isscalar(x)));
p.parse(varargin{:});
opt = p.Results;

nMeters  = opt.meters;
nSons    = opt.sons;

assert(nSons <= nMeters, "'sons' must be in [0, meters].");

% ---------- Paths ----------
projectRoot  = pwd;
customersDir = fullfile(projectRoot, "data", "virtual_customers");
outXlsx      = fullfile(projectRoot, "data", "transformer_profile.xlsx");

fprintf("\n=== synth_new_data ===\n");
fprintf("meters   : %d\n", nMeters);
fprintf("sons     : %d\n", nSons);
if isempty(opt.periods)
    fprintf("periods  : (not passed) -> generate_customers default\n");
else
    fprintf("periods  : %d\n", opt.periods);
end
if isempty(opt.seed)
    fprintf("seed     : (not passed) -> generate_customers rng('shuffle')\n");
else
    fprintf("seed     : %d (passed as rngSeed)\n", opt.seed);
end

% ---------- Clear old data ----------
clear_synthetic_data('projectRoot', projectRoot);

% ---------- Build argument list for generate_customers ----------
gcArgs = {nMeters, nSons};

if ~isempty(opt.periods)
    gcArgs = [gcArgs, {'periods', opt.periods}];
end

if ~isempty(opt.seed)
    gcArgs = [gcArgs, {'rngSeed', opt.seed}];
end

% ---------- Generate customers ----------
generate_customers(gcArgs{:});

% ---------- Sanity check ----------
d = dir(fullfile(customersDir, "*.xlsx"));
assert(~isempty(d), "No customer XLSX files generated.");

% ---------- Generate transformer ----------
generate_transformer_profile(customersDir, outXlsx);

fprintf("✔ Synthetic data generated successfully\n\n");
end
