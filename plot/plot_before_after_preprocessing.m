function fig = plot_before_after_preprocessing(varargin)
% plot_before_after_preprocessing
% Compare ORIGINAL synthetic data vs PREPROCESSED (corrupted) data:
%   - Transformer: overlay y(t) original vs preprocessed
%   - One customer: overlay x_i(t) original vs preprocessed
%   - Marks NaNs (missing) on the plots
%
% Usage:
%   plot_before_after_preprocessing()
%   plot_before_after_preprocessing('customerIdx', 7)
%   plot_before_after_preprocessing('customersSubdirOut', fullfile("data","virtual_customers_pp"), ...
%                                   'transformerFileOut', fullfile("data","transformer_profile_pp.xlsx"))

p = inputParser;
p.addParameter('projectRoot', "", @(x)isstring(x)||ischar(x));

% ORIGINAL
p.addParameter('customersSubdirIn', fullfile("data","virtual_customers"), @(x)isstring(x)||ischar(x));
p.addParameter('transformerFileIn', fullfile("data","transformer_profile.xlsx"), @(x)isstring(x)||ischar(x));

% PREPROCESSED (output of preprocess_synthetic_data)
p.addParameter('customersSubdirOut', fullfile("data","virtual_customers_pp"), @(x)isstring(x)||ischar(x));
p.addParameter('transformerFileOut', fullfile("data","transformer_profile_pp.xlsx"), @(x)isstring(x)||ischar(x));

% Which customer to show
p.addParameter('customerIdx', 1, @(x)isnumeric(x)&&isscalar(x)&&x>=1);

% Plot options
p.addParameter('showDiff', true, @(x)islogical(x)||ismember(x,[0 1]));
p.addParameter('title', "", @(x)isstring(x)||ischar(x));
p.parse(varargin{:});
opt = p.Results;

% Resolve projectRoot like your other funcs
if strlength(string(opt.projectRoot))==0
    thisFile = mfilename('fullpath');
    thisDir  = fileparts(thisFile);
    projectRoot = string(fileparts(thisDir));
else
    projectRoot = string(opt.projectRoot);
end

custInDir  = fullfile(projectRoot, string(opt.customersSubdirIn));
custOutDir = fullfile(projectRoot, string(opt.customersSubdirOut));

trInFile   = fullfile(projectRoot, string(opt.transformerFileIn));
trOutFile  = fullfile(projectRoot, string(opt.transformerFileOut));

assert(isfolder(custInDir),  "Missing original customers folder: %s", custInDir);
assert(isfolder(custOutDir), "Missing preprocessed customers folder: %s", custOutDir);
assert(isfile(trInFile),     "Missing original transformer file: %s", trInFile);
assert(isfile(trOutFile),    "Missing preprocessed transformer file: %s", trOutFile);

% Pick customer file by sorted order customer_*.xlsx
filesIn  = dir(fullfile(custInDir,  "customer_*.xlsx"));
filesOut = dir(fullfile(custOutDir, "customer_*.xlsx"));
assert(~isempty(filesIn) && ~isempty(filesOut), "customer_*.xlsx not found in input/output folders.");

[~,ix] = sort({filesIn.name});  filesIn  = filesIn(ix);
[~,ox] = sort({filesOut.name}); filesOut = filesOut(ox);

N = min(numel(filesIn), numel(filesOut));
idx = min(opt.customerIdx, N);

custFileIn  = fullfile(filesIn(idx).folder,  filesIn(idx).name);
custFileOut = fullfile(filesOut(idx).folder, filesOut(idx).name);

% ----- Load transformer -----
Ttr0 = readtable(trInFile);
Ttr1 = readtable(trOutFile);
assert(all(ismember({'STARTOFINTERVAL_TIMES','EC'}, Ttr0.Properties.VariableNames)), "Bad transformerFileIn columns");
assert(all(ismember({'STARTOFINTERVAL_TIMES','EC'}, Ttr1.Properties.VariableNames)), "Bad transformerFileOut columns");

t0 = toDatetimeCol(Ttr0.STARTOFINTERVAL_TIMES);
y0 = double(Ttr0.EC(:));
t1 = toDatetimeCol(Ttr1.STARTOFINTERVAL_TIMES);
y1 = double(Ttr1.EC(:));

% Align transformer on intersection times
tt0 = timetable(t0(:), y0, 'VariableNames', {'Y0'});
tt1 = timetable(t1(:), y1, 'VariableNames', {'Y1'});
ttT = synchronize(tt0, tt1, 'intersection');

tT = ttT.Time;
y0a = ttT.Y0;
y1a = ttT.Y1;

% ----- Load customer -----
Tc0 = readtable(custFileIn);
Tc1 = readtable(custFileOut);
assert(all(ismember({'STARTOFINTERVAL_TIMES','EC'}, Tc0.Properties.VariableNames)), "Bad customer columns (in)");
assert(all(ismember({'STARTOFINTERVAL_TIMES','EC'}, Tc1.Properties.VariableNames)), "Bad customer columns (out)");

tc0 = toDatetimeCol(Tc0.STARTOFINTERVAL_TIMES);
x0  = double(Tc0.EC(:));
tc1 = toDatetimeCol(Tc1.STARTOFINTERVAL_TIMES);
x1  = double(Tc1.EC(:));

ttc0 = timetable(tc0(:), x0, 'VariableNames', {'X0'});
ttc1 = timetable(tc1(:), x1, 'VariableNames', {'X1'});
ttC  = synchronize(ttc0, ttc1, 'intersection');

tC = ttC.Time;
x0a = ttC.X0;
x1a = ttC.X1;

% ----- Plot -----
fig = figure('Color','w','Name','Before vs After Preprocessing');
mainTitle = string(opt.title);
if strlength(mainTitle)==0
    mainTitle = sprintf("Before vs After Preprocessing | Customer %d (%s)", idx, filesIn(idx).name);
end
sgtitle(mainTitle, 'FontWeight','bold');

% Transformer overlay
subplot(2,2,1);
plot(tT, y0a, 'k', 'DisplayName','Transformer (original)'); hold on;
plot(tT, y1a, 'r', 'DisplayName','Transformer (preprocessed)');
markNaNs(tT, y1a, [1 0 0]); % mark NaNs in preprocessed
grid on; legend('Location','best'); ylabel('EC'); xlabel('Time');
title('Transformer: original vs preprocessed');

% Customer overlay
subplot(2,2,2);
plot(tC, x0a, 'k', 'DisplayName','Customer (original)'); hold on;
plot(tC, x1a, 'r', 'DisplayName','Customer (preprocessed)');
markNaNs(tC, x1a, [1 0 0]);
grid on; legend('Location','best'); ylabel('EC'); xlabel('Time');
title(sprintf('Customer %d: original vs preprocessed', idx));

% Optional difference plots
if opt.showDiff
    subplot(2,2,3);
    dY = y1a - y0a;
    plot(tT, dY, 'b'); hold on;
    yline(0,'k-');
    grid on; ylabel('\Delta EC'); xlabel('Time');
    title('Transformer difference (pre - orig)');

    subplot(2,2,4);
    dX = x1a - x0a;
    plot(tC, dX, 'b'); hold on;
    yline(0,'k-');
    grid on; ylabel('\Delta EC'); xlabel('Time');
    title(sprintf('Customer %d difference (pre - orig)', idx));
else
    subplot(2,2,3);
    histogram(y1a(~isnan(y1a)) - y0a(~isnan(y0a)));
    grid on; title('Transformer diff histogram'); xlabel('\Delta EC'); ylabel('count');

    subplot(2,2,4);
    histogram(x1a(~isnan(x1a)) - x0a(~isnan(x0a)));
    grid on; title(sprintf('Customer %d diff histogram', idx)); xlabel('\Delta EC'); ylabel('count');
end

end

% -------- helpers --------
function dt = toDatetimeCol(col)
    if isdatetime(col), dt = col(:); else, dt = datetime(col); dt = dt(:); end
end

function markNaNs(t, y, color)
    idx = isnan(y);
    if any(idx)
        yl = ylim;
        scatter(t(idx), yl(1)*ones(sum(idx),1), 10, 'filled', ...
            'MarkerFaceColor', color, 'MarkerEdgeColor', color, ...
            'DisplayName','Missing (NaN)');
        ylim(yl);
    end
end
