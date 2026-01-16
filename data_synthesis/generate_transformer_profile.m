function Ttr = generate_transformer_profile(customersDir, outXlsx)
% generate_transformer_profile
%   Reads customer Excel logs sequentially and aggregates (sums) EC per timestamp
%   to create a transformer load profile.
%
% Expected customer file format:
%   - STARTOFINTERVAL_TIMES (datetime or text convertible to datetime)
%   - EC (numeric)
%
% Inputs:
%   customersDir - folder containing customer_###.xlsx files
%   outXlsx      - output transformer Excel path (e.g. ".../data/transformer_001.xlsx")
%
% Output:
%   Ttr          - table with columns: STARTOFINTERVAL_TIMES, EC

    if nargin < 2
        outXlsx = fullfile(customersDir, "transformer_profile.xlsx");
    end

    files = dir(fullfile(customersDir, "*.xlsx"));
    % Ignore metadata file if exists
    files = files(~strcmpi({files.name}, "customers_metadata.xlsx"));

    if isempty(files)
        error("No customer .xlsx files found in: %s", customersDir);
    end

    aggTT = [];  % will become a timetable with variable 'EC'

    for k = 1:numel(files)
        fpath = fullfile(files(k).folder, files(k).name);

        T = readtable(fpath);

        % Validate columns
        req = {'STARTOFINTERVAL_TIMES','EC'};
        if ~all(ismember(req, T.Properties.VariableNames))
            error("File %s missing required columns: STARTOFINTERVAL_TIMES, EC", files(k).name);
        end

        % Parse time and values
        t = T.STARTOFINTERVAL_TIMES;
        if ~isdatetime(t)
            t = datetime(t);
        end
        ec = double(T.EC);

        % Build timetable for this customer
        tt = timetable(t(:), ec(:), 'VariableNames', {'EC'});

        % Aggregate
        if isempty(aggTT)
            aggTT = tt;
        else
            % Align by timestamps (outer join), then sum (missing => 0)
            J = synchronize(aggTT, tt, 'union', 'fillwithmissing');
            a = J.EC_aggTT;  a(isnan(a)) = 0;
            b = J.EC_tt;     b(isnan(b)) = 0;

            aggTT = timetable(J.Time, a + b, 'VariableNames', {'EC'});
        end

        fprintf("Aggregated %d/%d: %s\n", k, numel(files), files(k).name);
    end

    % Output table in the same template format
    Ttr = table(aggTT.Time, aggTT.EC, 'VariableNames', {'STARTOFINTERVAL_TIMES','EC'});

    % Save
    writetable(Ttr, outXlsx, 'FileType','spreadsheet');

    fprintf("Transformer profile saved to: %s\n", outXlsx);
end
