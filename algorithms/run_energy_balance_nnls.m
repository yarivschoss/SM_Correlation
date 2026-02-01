function results = run_energy_balance_nnls(varargin)
% run_energy_balance_nnls
% Baseline NNLS for energy balance:
%   min ||X w - y||^2  s.t.  w >= 0
%
% Default project structure:
%   data/transformer_profile.xlsx
%   data/virtual_customers/customer_*.xlsx
%   data/virtual_customers/customers_metadata.xlsx (optional; is_son)
%
% Usage:
%   results = run_energy_balance_nnls();
%   results = run_energy_balance_nnls('projectRoot', pwd, 'doPlots', true);

% -------------------- Parse inputs --------------------
p = inputParser;
p.addParameter('projectRoot', "", @(x)isstring(x)||ischar(x));
p.addParameter('customersSubdir', fullfile("data","virtual_customers"), @(x)isstring(x)||ischar(x));
p.addParameter('transformerFile', fullfile("data","transformer_profile.xlsx"), @(x)isstring(x)||ischar(x));
p.addParameter('doPlots', true, @(x)islogical(x)||ismember(x,[0 1]));

% Classification options (optional)
p.addParameter('fixedTau', [], @(x)isnumeric(x)&&isscalar(x) || isempty(x));      % if empty -> can auto-pick (if GT exists)
p.addParameter('autoTau', true, @(x)islogical(x)||ismember(x,[0 1]));        % only used when tau is empty and GT exists

p.parse(varargin{:});
opt = p.Results;

% -------------------- Resolve project root --------------------
if strlength(string(opt.projectRoot))==0
    thisFile = mfilename('fullpath');     % .../algorithms/run_energy_balance_nnls.m
    thisDir  = fileparts(thisFile);       % .../algorithms
    opt.projectRoot = string(fileparts(thisDir)); % project root = parent
end
projectRoot = string(opt.projectRoot);

customersDir    = fullfile(projectRoot, string(opt.customersSubdir));
transformerXls  = fullfile(projectRoot, string(opt.transformerFile));

% -------------------- Load transformer --------------------
Ttr = readtable(transformerXls);
assert(all(ismember({'STARTOFINTERVAL_TIMES','EC'}, Ttr.Properties.VariableNames)), ...
    'Transformer file must include STARTOFINTERVAL_TIMES and EC');

t_tr = toDatetimeCol(Ttr.STARTOFINTERVAL_TIMES);
y    = double(Ttr.EC(:));

tt_tr = timetable(t_tr(:), y, 'VariableNames', {'Y'});

% -------------------- Load customers -> build X aligned --------------------
files = dir(fullfile(customersDir, "customer_*.xlsx"));
if isempty(files)
    error("No customer_*.xlsx files found in: %s", customersDir);
end
[~,ix] = sort({files.name});
files = files(ix);
N = numel(files);

X = zeros(height(tt_tr), N);

for i = 1:N
    Tc = readtable(fullfile(files(i).folder, files(i).name));
    assert(all(ismember({'STARTOFINTERVAL_TIMES','EC'}, Tc.Properties.VariableNames)), ...
        'Customer file %s must include STARTOFINTERVAL_TIMES and EC', files(i).name);

    t_c = toDatetimeCol(Tc.STARTOFINTERVAL_TIMES);
    x_c = double(Tc.EC(:));

    tt_c = timetable(t_c(:), x_c, 'VariableNames', {'X'});
    ttJ  = synchronize(tt_tr, tt_c, 'intersection');

    if i == 1
        % Update master time to the intersection
        tt_tr = ttJ(:, 'Y');
        X = zeros(height(tt_tr), N);
    else
        % Keep aligning to current master
        tt_master = tt_tr;
        ttJ = synchronize(tt_master, tt_c, 'intersection');
        if height(ttJ) ~= height(tt_tr)
            % shrink master + X if needed
            keepTimes = ttJ.Time;
            X = X(ismember(tt_tr.Time, keepTimes), :);
            tt_tr = tt_tr(ismember(tt_tr.Time, keepTimes), :);
        end
    end

    X(:,i) = ttJ.X;
end

y = tt_tr.Y;

% -------------------- NNLS solve --------------------
w_hat = lsqnonneg(X, y);

% -------------------- Reconstruction --------------------
y_hat = X * w_hat;
rmse  = sqrt(mean((y - y_hat).^2, 'omitnan'));

% -------------------- Optional GT + classification --------------------
cls = struct();
metaPath = fullfile(customersDir, "customers_metadata.xlsx");
if isfile(metaPath)
    M = readtable(metaPath);
    vars = lower(string(M.Properties.VariableNames));
    idx_is_son = find(vars=="is_son" | vars=="isson", 1);

    if ~isempty(idx_is_son)
        gt_isSon = logical(M{:, idx_is_son});
        gt_isSon = gt_isSon(:);

        if numel(gt_isSon) ~= N
            warning("GT is_son length (%d) != N customers (%d). Skipping classification.", numel(gt_isSon), N);
        else
            cls.gt_isSon = gt_isSon;

            if ~isempty(opt.fixedTau)
                tau_use = opt.fixedTau;
                cls.tau_source = "fixed";
            else
                if opt.autoTau
                    % simple mid-gap (works well when separated)
                    if any(~gt_isSon) && any(gt_isSon)
                        tau_use = 0.5 * (max(w_hat(~gt_isSon)) + min(w_hat(gt_isSon)));
                    else
                        tau_use = 0.5;
                    end
                    cls.tau_source = "auto_midgap";
                else
                    tau_use = 0.5;
                    cls.tau_source = "default_0.5";
                end
            end

            pred_isSon = (w_hat >= tau_use);

            cm.TP = sum(pred_isSon & gt_isSon);
            cm.FP = sum(pred_isSon & ~gt_isSon);
            cm.TN = sum(~pred_isSon & ~gt_isSon);
            cm.FN = sum(~pred_isSon & gt_isSon);

            Prec = cm.TP / max(1, cm.TP + cm.FP);
            Rec  = cm.TP / max(1, cm.TP + cm.FN);
            F1   = 2*(Prec*Rec) / max(1e-12, Prec+Rec);

            cls.tau = tau_use;
            cls.pred_isSon = pred_isSon;
            cls.cm = cm;
            cls.Precision = Prec;
            cls.Recall = Rec;
            cls.F1 = F1;
        end
    end
end

% -------------------- Pack results --------------------
results = struct();
results.projectRoot = projectRoot;
results.customersDir = customersDir;
results.transformerFile = transformerXls;
results.customerFiles = string({files.name}).';

results.time = tt_tr.Time;
results.y = y;
results.X = X;

results.w_hat = w_hat;
results.y_hat = y_hat;
results.rmse  = rmse;
results.isSon = pred_isSon;

results.classification = cls;

% -------------------- Plots --------------------
if opt.doPlots
    figure('Color','w','Name','NNLS: Transformer Reconstruction');
    plot(results.time, results.y, 'k', 'DisplayName','y (measured)'); hold on;
    plot(results.time, results.y_hat, 'r', 'DisplayName','yhat (NNLS)');
    grid on; xlabel('Time'); ylabel('Energy');
    title(sprintf('NNLS reconstruction (RMSE=%.4g)', results.rmse));
    legend('Location','best');

    figure('Color','w','Name','NNLS: Weights');
    stem(results.w_hat, 'filled');
    grid on; xlabel('Customer index'); ylabel('w');
    title('NNLS weights');

    if isfield(cls,'F1')
        fprintf("\nNNLS classification (tau=%g, source=%s): Prec=%.3f Rec=%.3f F1=%.3f\n", ...
            cls.tau, cls.tau_source, cls.Precision, cls.Recall, cls.F1);
    end
end

end

% ---------- helper ----------
function dt = toDatetimeCol(col)
if isdatetime(col)
    dt = col(:);
else
    dt = datetime(col);
    dt = dt(:);
end
end
