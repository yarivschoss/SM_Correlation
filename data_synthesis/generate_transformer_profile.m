function generate_transformer_profile(customersDir, outXlsx, varargin)
% generate_transformer_profile
% Aggregates customer energy logs to create transformer load profile.
%
% Inputs:
%   customersDir : folder containing customer_*.xlsx and optionally customers_metadata.xlsx
%   outXlsx      : output xlsx path for transformer profile
%
% Optional Name-Value:
%   'useMetadata' (logical) default true   - if true, use customers_metadata.xlsx when exists
%   'sonsOnly'    (logical) default true   - if true and metadata has is_son, aggregate only sons
%   'measNoiseStd'(double)  default 0      - additive Gaussian noise std on transformer EC (optional)
%   'lossAlpha'   (double)  default 0      - multiplicative losses: y = (1+alpha)*sum(...)
%   'timeColumn'  (string)  default "STARTOFINTERVAL_TIMES"
%   'ecColumn'    (string)  default "EC"

p = inputParser;
p.addParameter('useMetadata', true, @(x)islogical(x)||ismember(x,[0 1]));
p.addParameter('sonsOnly', true, @(x)islogical(x)||ismember(x,[0 1]));
p.addParameter('measNoiseStd', 0, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('lossAlpha', 0, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('timeColumn', "STARTOFINTERVAL_TIMES", @(x)isstring(x)||ischar(x));
p.addParameter('ecColumn', "EC", @(x)isstring(x)||ischar(x));
p.parse(varargin{:});
opt = p.Results;

customersDir = string(customersDir);
outXlsx = string(outXlsx);
timeCol = string(opt.timeColumn);
ecCol   = string(opt.ecColumn);

files = dir(fullfile(customersDir, "customer_*.xlsx"));
if isempty(files)
    error("No customer_*.xlsx files found in: %s", customersDir);
end
[~,ix] = sort({files.name});
files = files(ix);

N = numel(files);

% --- Decide which customers to include ---
includeMask = true(N,1);

metaPath = fullfile(customersDir, "customers_metadata.xlsx");
if opt.useMetadata && isfile(metaPath)
    M = readtable(metaPath);

    % find is_son column (case-insensitive)
    varNames = string(M.Properties.VariableNames);
    idx_is_son = find(lower(varNames) == "is_son" | lower(varNames) == "isson", 1);

    if ~isempty(idx_is_son) && opt.sonsOnly
        is_son = M{:, idx_is_son};
        is_son = logical(is_son);

        if numel(is_son) == N
            includeMask = is_son(:);
        else
            warning("customers_metadata.xlsx has is_son but length (%d) != num customers (%d). Using all customers.", ...
                numel(is_son), N);
        end
    end
end

useFiles = files(includeMask);
if isempty(useFiles)
    error("No customers selected for aggregation. (Maybe all is_son are false?)");
end

% --- Read first selected customer as time master ---
T0 = readtable(fullfile(useFiles(1).folder, useFiles(1).name));
assert(any(strcmpi(T0.Properties.VariableNames, timeCol)), "Missing time column in %s", useFiles(1).name);
assert(any(strcmpi(T0.Properties.VariableNames, ecCol)),   "Missing EC column in %s",   useFiles(1).name);

t_master = toDatetimeCol(T0.(findVar(T0, timeCol)));
y_sum    = double(T0.(findVar(T0, ecCol)));

tt_master = timetable(y_sum, 'RowTimes', t_master, 'VariableNames', {'SUM'});

% --- Aggregate others on intersection timestamps ---
for k = 2:numel(useFiles)
    Tk = readtable(fullfile(useFiles(k).folder, useFiles(k).name));
    assert(any(strcmpi(Tk.Properties.VariableNames, timeCol)), "Missing time column in %s", useFiles(k).name);
    assert(any(strcmpi(Tk.Properties.VariableNames, ecCol)),   "Missing EC column in %s",   useFiles(k).name);

    t_k = toDatetimeCol(Tk.(findVar(Tk, timeCol)));
    x_k = double(Tk.(findVar(Tk, ecCol)));

    tt_k      = timetable(x_k, 'RowTimes', t_k, 'VariableNames', {'X'});

    ttJ = synchronize(tt_master, tt_k, 'intersection');

    % if intersection shrinks, update master vectors accordingly
    t_master = ttJ.Properties.RowTimes;
    y_sum    = ttJ.SUM + ttJ.X;

    tt_master = timetable(y_sum, 'RowTimes', t_master, 'VariableNames', {'SUM'});
    
end

% --- Apply losses and optional measurement noise ---
y = (1 + opt.lossAlpha) * y_sum;

if opt.measNoiseStd > 0
    y = y + opt.measNoiseStd * randn(size(y));
end

% --- Write transformer profile ---
Tout = table();
Tout.(timeCol) = t_master;
Tout.(ecCol)   = y;

outDir = fileparts(outXlsx);
if ~isempty(outDir) && ~exist(outDir,'dir')
    mkdir(outDir);
end

writetable(Tout, outXlsx);

fprintf("Transformer profile saved: %s\n", outXlsx);
fprintf("Aggregated customers: %d / %d (sonsOnly=%d)\n", numel(useFiles), N, opt.sonsOnly);

end

% ---------- helpers ----------
function dt = toDatetimeCol(col)
    if isdatetime(col)
        dt = col(:);
    else
        dt = datetime(col);
        dt = dt(:);
    end
end

function vn = findVar(T, nameWanted)
% returns actual variable name (case-insensitive match)
    vars = string(T.Properties.VariableNames);
    idx = find(lower(vars) == lower(string(nameWanted)), 1);
    if isempty(idx)
        error("Column '%s' not found.", nameWanted);
    end
    vn = T.Properties.VariableNames{idx};
end
