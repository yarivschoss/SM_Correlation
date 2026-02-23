function results = run_energy_balance_ukf(varargin)
% run_energy_balance_ukf
% End-to-end UKF energy-balance estimator for synthetic transformer/customer data.
%
% Data layout (projectRoot):
%   data/
%     transformer_profile.xlsx                 (columns: STARTOFINTERVAL_TIMES, EC)
%     virtual_customers/
%        customer_001.xlsx ... customer_100.xlsx (columns: STARTOFINTERVAL_TIMES, EC)
%        customers_metadata.xlsx               (optional; may include is_son column)
%
% Model:
%   y(t) = (1 + alpha) * sum_i w_i * x_i(t) + v(t)
%   state s = [w_1..w_N, alpha]^T
%   random-walk: s_{k+1} = s_k + q_k
%
% Outputs:
%   results struct with fields:
%     w_hat, alpha_hat, y_hat, rmse, time, customers, etc.
%
% Optional classification (if GT available):
%   threshold selection on calibration set, evaluation on test set.

% -------------------- Parse inputs --------------------
p = inputParser;
p.addParameter('projectRoot', "", @(x)isstring(x)||ischar(x));
p.addParameter('customersSubdir', fullfile("data","virtual_customers"), @(x)isstring(x)||ischar(x));
p.addParameter('transformerFile', fullfile("data","transformer_profile.xlsx"), @(x)isstring(x)||ischar(x));
p.addParameter('doPlots', true, @(x)islogical(x)||ismember(x,[0 1]));
p.addParameter('fixedTau', [], @(x)isnumeric(x)&&isscalar(x) || isempty(x));
p.addParameter('rngSeed', [], @(x) isempty(x) || (isnumeric(x)&&isscalar(x)));

% UKF params
p.addParameter('Qw', 1e-6, @(x)isnumeric(x)&&isscalar(x));      % process noise for weights
p.addParameter('Qalpha', 1e-8, @(x)isnumeric(x)&&isscalar(x));  % process noise for alpha
p.addParameter('R', 1e-4, @(x)isnumeric(x)&&isscalar(x));       % meas noise variance
p.addParameter('alpha0', 0.01, @(x)isnumeric(x)&&isscalar(x));  % initial alpha
p.addParameter('w0', 0.5, @(x)isnumeric(x)&&isscalar(x));       % initial w for all
p.addParameter('wBounds', [0 1], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('alphaBounds', [-0.2 0.2], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('missingPolicy',"skip", @(x) any(strcmpi(string(x),["skip","hold"])));

% UKF sigma-point tuning
p.addParameter('ukf_alpha', 1e-3, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('ukf_beta', 2, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('ukf_kappa', 0, @(x)isnumeric(x)&&isscalar(x));

% Optional ground truth:
% - either provide 'isSon' logical vector length N (customers order),
% - or provide 'sonsIdx' indices
p.addParameter('isSon', [], @(x)islogical(x)||isempty(x));
p.addParameter('sonsIdx', [], @(x)isnumeric(x)||isempty(x));

% Split for calibration/test (by time index)
p.addParameter('calibFrac', 0.7, @(x)isnumeric(x)&&isscalar(x)&&x>0&&x<1);

p.parse(varargin{:});
opt = p.Results;

if isempty(opt.rngSeed)
    rng('shuffle');   % different seed every run
else
    rng(opt.rngSeed); % reproducible
end

% resolve project root
if strlength(string(opt.projectRoot))==0
    % default: assume this file is in project somewhere; project root = parent of this file's folder
    thisFile = mfilename('fullpath');
    thisDir  = fileparts(thisFile);
    opt.projectRoot = string(fileparts(thisDir));
end
projectRoot = string(opt.projectRoot);

customersDir   = fullfile(projectRoot, string(opt.customersSubdir));
transformerXls = fullfile(projectRoot, string(opt.transformerFile));

% -------------------- Load transformer --------------------
Ttr = readtable(transformerXls);
assert(all(ismember({'STARTOFINTERVAL_TIMES','EC'}, Ttr.Properties.VariableNames)), ...
    'Transformer file must include STARTOFINTERVAL_TIMES and EC');

t_tr = toDatetimeCol(Ttr.STARTOFINTERVAL_TIMES);
y    = double(Ttr.EC(:));

% -------------------- Load customers (sequential) --------------------
files = dir(fullfile(customersDir, "customer_*.xlsx"));
if isempty(files)
    error("No customer_*.xlsx files found in: %s", customersDir);
end
% sort by name
[~,ix] = sort({files.name});
files = files(ix);

Ncust = numel(files); % Ncust = number of customers

% Build X matrix aligned to transformer timestamps
X = zeros(numel(t_tr), Ncust); % t_tr rows with Ncust columns

for i = 1:Ncust
    fpath = fullfile(files(i).folder, files(i).name);
    Tc = readtable(fpath);
    assert(all(ismember({'STARTOFINTERVAL_TIMES','EC'}, Tc.Properties.VariableNames)), ...
        'Customer file %s must include STARTOFINTERVAL_TIMES and EC', files(i).name);

    t_c = toDatetimeCol(Tc.STARTOFINTERVAL_TIMES);
    x_c = double(Tc.EC(:));

    % align to transformer time axis
    % use timetable synchronize for robust alignment
    tt_tr = timetable(t_tr(:), y, 'VariableNames', {'Y'});
    tt_c  = timetable(t_c(:),  x_c, 'VariableNames', {'X'});

    ttJ = synchronize(tt_tr, tt_c, 'intersection'); % keep only common timestamps
    if i == 1
        % update master time/y to intersection (so all customers align identically)
        t = ttJ.Time;
        y_al = ttJ.Y;
        % rebuild transformer vectors
        t_tr = t;
        y    = y_al;
        X    = zeros(numel(t_tr), Ncust);
    end

    % for customers after i=1, force intersection with current master time
    tt_master = timetable(t_tr(:), y, 'VariableNames', {'Y'});
    ttJ2 = synchronize(tt_master, tt_c, 'intersection');
    if numel(ttJ2.Time) ~= numel(t_tr)
        % if mismatch, shrink all to the new intersection
        t_tr = ttJ2.Time;
        y    = ttJ2.Y;
        X    = X(ismember(tt_master.Time, t_tr), :);
    end

    X(:,i) = ttJ2.X;
end

T = numel(t_tr);

% -------------------- Ground Truth (optional) --------------------
isSon = [];
if ~isempty(opt.isSon)
    isSon = opt.isSon(:);
elseif ~isempty(opt.sonsIdx)
    isSon = false(Ncust,1);
    isSon(opt.sonsIdx) = true;
else
    % try read from customers_metadata.xlsx if exists
    metaPath = fullfile(customersDir, "customers_metadata.xlsx");
    if isfile(metaPath)
        M = readtable(metaPath);
        % Accept either 'is_son' or 'isSon'
        if any(strcmpi(M.Properties.VariableNames, 'is_son'))
            isSon = logical(M{:, strcmpi(M.Properties.VariableNames,'is_son')});
        elseif any(strcmpi(M.Properties.VariableNames, 'isSon'))
            isSon = logical(M{:, strcmpi(M.Properties.VariableNames,'isSon')});
        end
        if ~isempty(isSon) && numel(isSon) ~= Ncust
            warning("GT isSon exists but length != N customers; ignoring GT.");
            isSon = [];
        end
    end
end

% -------------------- UKF init --------------------
% state: [w(1..N), alpha]
n = Ncust + 1;

w0 = opt.w0 * ones(Ncust,1);
alpha0 = opt.alpha0;

xk = [w0; alpha0];

% covariance init
P = diag([1e-1*ones(Ncust,1); 1e-3]);  % fairly loose

% process noise
Q = diag([opt.Qw*ones(Ncust,1); opt.Qalpha]);

% measurement noise
R = opt.R;

% sigma-point params
ukf_a = opt.ukf_alpha;
ukf_b = opt.ukf_beta;
ukf_k = opt.ukf_kappa;

% storage
w_hist = zeros(Ncust, T);
alpha_hist = zeros(1, T);
y_hat = zeros(T,1);

% -------------------- Run UKF over time --------------------
for k = 1:T

    x_in = X(k,:).';     % customers at time k  (Ncust x 1)
    yk   = y(k);         % transformer measurement

    % ---------- Handle missing values ----------
    hasMissing = isnan(yk) || any(isnan(x_in));

    if hasMissing
        switch lower(string(opt.missingPolicy))

            case "skip"
                % Prediction only (no measurement update)
                [x_pred, P_pred] = ukf_predict(xk, P, Q, ukf_a, ukf_b, ukf_k);

                xk = x_pred;
                P  = P_pred;

                y_hat(k) = NaN;

                w_hist(:,k)   = xk(1:Ncust);
                alpha_hist(k) = xk(end);
                continue;

            case "hold"
                % Replace missing with previous available sample
                if k > 1
                    if isnan(yk)
                        yk = y(k-1);
                    end

                    miss = isnan(x_in);
                    if any(miss)
                        x_in(miss) = X(k-1, miss);
                    end
                else
                    % First sample missing → skip
                    [x_pred, P_pred] = ukf_predict(xk, P, Q, ukf_a, ukf_b, ukf_k);

                    xk = x_pred;
                    P  = P_pred;

                    y_hat(k) = NaN;

                    w_hist(:,k)   = xk(1:Ncust);
                    alpha_hist(k) = xk(end);
                    continue;
                end

            otherwise
                error("Unknown missingPolicy: %s", opt.missingPolicy);
        end
    end

    % ---------- UKF Predict ----------
    [x_pred, P_pred] = ukf_predict(xk, P, Q, ukf_a, ukf_b, ukf_k);

    % ---------- UKF Update ----------
    h = @(s) meas_model(s, x_in);  % scalar measurement model
    [x_upd, P_upd, yk_hat] = ukf_update(x_pred, P_pred, yk, R, h, ukf_a, ukf_b, ukf_k);

    % ---------- Clamp to physical bounds ----------
    x_upd(1:Ncust) = min(max(x_upd(1:Ncust), opt.wBounds(1)), opt.wBounds(2));
    x_upd(end)     = min(max(x_upd(end), opt.alphaBounds(1)), opt.alphaBounds(2));

    % ---------- Save state ----------
    xk = x_upd;
    P  = P_upd;

    w_hist(:,k)   = xk(1:Ncust);
    alpha_hist(k) = xk(end);
    y_hat(k)      = yk_hat;

end


% final estimates (e.g., average over last chunk for stability)
w_hat = mean(w_hist(:, max(1,T-200):T), 2);
alpha_hat = mean(alpha_hist(max(1,T-200):T));

rmse = sqrt(mean((y - y_hat).^2, 'omitnan'));

% -------------------- Optional: threshold selection + classification --------------------
cls = struct();
if ~isempty(isSon)
    calibT = floor(opt.calibFrac*T);
    idxCal = 1:calibT;
    idxTes = (calibT+1):T;

    % Use summary score per customer: average w over interval
    w_cal = mean(w_hist(:, idxCal), 2);
    w_tes = mean(w_hist(:, idxTes), 2);

    if isempty(opt.fixedTau)
    % choose tau on calibration to maximize F1
    taus = linspace(0,1,201);
    best = struct('tau',0.5,'F1',-inf,'Prec',0,'Rec',0,'cm',[]);
    for tt = taus
        pred = w_cal >= tt;
        cm = confusionCounts(isSon, pred);
        [Prec, Rec, F1] = prfFromCM(cm);
        if F1 > best.F1
            best.tau = tt; best.F1 = F1; best.Prec = Prec; best.Rec = Rec; best.cm = cm;
        end
    end

    tau_use = best.tau;
    cls.calibration = best;

    else
        tau_use = opt.fixedTau;   % from last process
        cls.calibration = [];
    end

    % evaluate on test
    predTest = w_tes >= tau_use;
    cmTest = confusionCounts(isSon, predTest);
    [PrecT, RecT, F1T] = prfFromCM(cmTest);

    cls.tau = tau_use;
    cls.test.cm = cmTest;
    cls.test.Precision = PrecT;
    cls.test.Recall = RecT;
    cls.test.F1 = F1T;
    cls.w_cal = w_cal;
    cls.w_test = w_tes;
end

% -------------------- Pack results --------------------
results = struct();
results.projectRoot = projectRoot;
results.customersDir = customersDir;
results.transformerFile = transformerXls;

results.time = t_tr;
results.y = y;
results.X = X;
results.customerFiles = string({files.name}).';

results.w_hist = w_hist;
results.alpha_hist = alpha_hist;
results.w_hat = w_hat;
results.alpha_hat = alpha_hat;

results.y_hat = y_hat;
results.rmse = rmse;

results.classification = cls;

if isfield(cls,'tau') && ~isempty(cls.tau)
    results.isSon = (w_hat >= cls.tau);
else
    results.isSon = [];
end

% -------------------- Plots --------------------
if opt.doPlots
    figure('Color','w','Name','Transformer: Measured vs Reconstructed');
    plot(t_tr, y, 'DisplayName','y (measured)'); hold on;
    plot(t_tr, y_hat, 'DisplayName','yhat (UKF)');
    grid on; xlabel('Time'); ylabel('Energy');
    title(sprintf('Transformer profile (RMSE = %.4g)', rmse));
    legend('Location','best');

    figure('Color','w','Name','Estimated weights (w_i)');
    plot(w_hist.');
    grid on; xlabel('Time index'); ylabel('w_i');
    title('UKF estimated weights over time');

    figure('Color','w','Name','Estimated alpha');
    plot(alpha_hist);
    grid on; xlabel('Time index'); ylabel('alpha');
    title('UKF estimated loss factor \alpha(t)');

    if ~isempty(isSon)
        figure('Color','w','Name','w score with GT');
        stem(w_hat, 'DisplayName','GT orphans'); hold on;
        stem(find(isSon), w_hat(isSon), 'DisplayName','GT sons');
        grid on; xlabel('Customer index'); ylabel('w\_hat');
        title('Final Weights');
        legend('Location','best');

        if isfield(cls,'tau') && ~isempty(cls.tau)
           yline(cls.tau,'k--','DisplayName','\tau');
        end

        fprintf('\nClassification (if GT available):\n');
        fprintf('  tau* = %.3f\n', cls.tau);
        fprintf('  Test Precision=%.3f Recall=%.3f F1=%.3f\n', cls.test.Precision, cls.test.Recall, cls.test.F1);
    end
end

end

% ========================= Helper functions =========================

function yk = meas_model(s, x_in)
% s = [w; alpha], x_in = customers at time k
w = s(1:end-1);
alpha = s(end);
yk = (1 + alpha) * sum(w .* x_in);
end

function dt = toDatetimeCol(col)
if isdatetime(col)
    dt = col(:);
else
    dt = datetime(col);
    dt = dt(:);
end
end

function [x_pred, P_pred] = ukf_predict(x, P, Q, a, b, kappa)
% random walk: x_{k+1} = x_k
% so predicted mean same, covariance adds Q
x_pred = x;
P_pred = P + Q;
% keep symmetric
P_pred = (P_pred + P_pred')/2;
end

function [x_upd, P_upd, yhat] = ukf_update(x_pred, P_pred, y, R, h, a, b, kappa)
% UKF measurement update for scalar measurement
n = numel(x_pred);

[Xsig, Wm, Wc] = sigmaPoints(x_pred, P_pred, a, b, kappa);

% propagate through measurement
m = size(Xsig,2);
Ysig = zeros(1,m);
for i = 1:m
    Ysig(i) = h(Xsig(:,i));
end

yhat = sum(Wm .* Ysig);

% innovation covariance S and cross covariance Pxy
S = 0;
Pxy = zeros(n,1);
for i = 1:m
    dy = Ysig(i) - yhat;
    dx = Xsig(:,i) - x_pred;
    S = S + Wc(i) * (dy*dy);
    Pxy = Pxy + Wc(i) * (dx * dy);
end
S = S + R;

% Kalman gain
K = Pxy / S;

% update
x_upd = x_pred + K * (y - yhat);
P_upd = P_pred - (K * S * K');

% stabilize
P_upd = (P_upd + P_upd')/2;
P_upd = nearestSPD(P_upd);
end

function [Xsig, Wm, Wc] = sigmaPoints(x, P, a, b, kappa)
n = numel(x);
lambda = a^2 * (n + kappa) - n;

% Cholesky (robust)
S = chol(nearestSPD(P), 'lower');

Xsig = zeros(n, 2*n+1);
Xsig(:,1) = x;

gamma = sqrt(n + lambda);
for i = 1:n
    Xsig(:, 1+i)     = x + gamma * S(:,i);
    Xsig(:, 1+n+i)   = x - gamma * S(:,i);
end

Wm = zeros(1,2*n+1);
Wc = zeros(1,2*n+1);

Wm(1) = lambda / (n + lambda);
Wc(1) = Wm(1) + (1 - a^2 + b);

for i = 2:(2*n+1)
    Wm(i) = 1 / (2*(n + lambda));
    Wc(i) = Wm(i);
end
end

function A = nearestSPD(A)
% Make matrix symmetric positive definite (simple approach)
A = (A + A')/2;
% jitter if needed
[~,p] = chol(A);
if p ~= 0
    % add diagonal until SPD
    diagAdd = 1e-12;
    for k = 1:20
        A = A + diagAdd*eye(size(A));
        [~,p] = chol(A);
        if p == 0, break; end
        diagAdd = diagAdd*10;
    end
end
end

function cm = confusionCounts(gt, pred)
% gt, pred are logical vectors
gt = logical(gt(:));
pred = logical(pred(:));
cm.TP = sum(pred & gt);
cm.FP = sum(pred & ~gt);
cm.TN = sum(~pred & ~gt);
cm.FN = sum(~pred & gt);
end

function [Prec, Rec, F1] = prfFromCM(cm)
Prec = cm.TP / max(1, (cm.TP + cm.FP));
Rec  = cm.TP / max(1, (cm.TP + cm.FN));
F1   = 2 * (Prec*Rec) / max(1e-12, (Prec + Rec));
end
