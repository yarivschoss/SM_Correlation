function T = generate_from_indexed_NE(E_i, n_i, startTime, outXlsx)
% generate_excel_log_from_indexed_NE
%   Generates a continuous quarter-hourly synthetic energy log (18000 samples)
%   whose empirical distribution matches the given indexed N(E), and saves it
%   to an Excel file with columns:
%       STARTOFINTERVAL_TIMES, EC
%
% Inputs:
%   E_i      - energy axis vector E(i)
%   n_i      - PDF values vector N(E(i)) (nonnegative)
%   startTime- datetime scalar (first timestamp)
%   outXlsx  - output Excel file path, e.g. "synthetic_customer.xlsx"
%
% Output:
%   T - table written to Excel

    N = 18000;
    dt = minutes(15);

    % --- Validation ---
    E_i = E_i(:);
    n_i = n_i(:);

    if numel(E_i) ~= numel(n_i)
        error("E_i and n_i must have the same length.");
    end
    if any(n_i < 0)
        error("n_i must be nonnegative.");
    end
    if ~isdatetime(startTime) || ~isscalar(startTime)
        error("startTime must be a datetime scalar.");
    end
    if ~(isstring(outXlsx) || ischar(outXlsx))
        error("outXlsx must be a string/char path.");
    end

    % --- Normalize to discrete PMF ---
    if sum(n_i) == 0
        error("n_i sums to zero. Cannot sample from an empty distribution.");
    end
    p = n_i / sum(n_i);

    % --- Inverse transform sampling on the discrete grid ---
    cdf = cumsum(p);
    cdf(end) = 1;

    u = rand(N,1);
    idx = arrayfun(@(x) find(cdf >= x, 1, 'first'), u);
    X = E_i(idx);

    % --- Physical constraint (optional) ---
    X = max(X, 0);

    % --- Continuous timestamps ---
    t = startTime + (0:N-1)' * dt;

    % --- Build table with template column names ---
    T = table(t, X, 'VariableNames', {'STARTOFINTERVAL_TIMES','EC'});

    % --- Save to Excel ---
    [~,~,ext] = fileparts(outXlsx);
    if ~strcmpi(ext, ".xlsx") && ~strcmpi(ext, ".xls")
        error("outXlsx must be an Excel file path ending with .xlsx or .xls");
    end
    writetable(T, outXlsx, 'FileType', 'spreadsheet');
end
