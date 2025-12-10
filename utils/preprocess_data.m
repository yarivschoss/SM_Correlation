function data = preprocess_data(rawData, cfg)
% preprocess_data
%   Cleans and preprocesses smart-meter time series for all fields in rawData.
%   For each series (e.g., transformer, yard), the function:
%       1) aligns data on a regular time grid
%       2) detects obvious bad data and outliers
%       3) handles missing values and anomalies via interpolation
%       4) optionally smooths the signal
%
%   Input:
%       rawData - struct; each field is a struct with:
%                   .timeVec    (datetime vector)
%                   .voltVec    (numeric vector)
%                   .timeToVolt (containers.Map)
%       cfg     - global configuration struct
%
%   Output:
%       data    - struct; same fields as rawData, each containing:
%                   .timeVec
%                   .value          (cleaned + smoothed)
%                   .valueRaw       (aligned but uncleaned)
%                   .maskMissingOrig
%                   .maskAnomaly
%                   .maskLongGap
%                   .maskValid

    % ---------- Preprocessing defaults (if not present in cfg) ----------
    if ~isfield(cfg, "preproc")
        cfg.preproc = struct();
    end

    if ~isfield(cfg.preproc, "maxInterpGapSamples")
        cfg.preproc.maxInterpGapSamples = 4;      % e.g. up to 1 hour for 15-min data
    end
    if ~isfield(cfg.preproc, "iqrMultiplier")
        cfg.preproc.iqrMultiplier = 3;            % IQR multiplier for outlier detection
    end
    if ~isfield(cfg.preproc, "zeroRunThresholdSamples")
        cfg.preproc.zeroRunThresholdSamples = 96; % e.g. 1 day for 15-min data
    end
    if ~isfield(cfg.preproc, "smoothWindowSamples")
        cfg.preproc.smoothWindowSamples = 3;      % moving median window
    end

    dt = seconds(cfg.sampleTime_sec);

    seriesNames = fieldnames(rawData);
    data = struct();

    for k = 1:numel(seriesNames)
        name = seriesNames{k};
        s = rawData.(name);

        % Defensive check: require at least timeVec and voltVec
        if ~isfield(s, "timeVec") || ~isfield(s, "voltVec")
            warning("preprocess_data:MissingFields", ...
                    "Series '%s' does not contain timeVec/voltVec. Skipped.", name);
            continue;
        end

        [timeVec, valueClean, valueRaw, ...
         maskMissingOrig, maskAnomaly, maskLongGap, maskValid] = ...
            clean_single_series(s.timeVec, s.voltVec, dt, cfg.preproc);

        % Store results
        data.(name).timeVec        = timeVec;
        data.(name).value          = valueClean;
        data.(name).valueRaw       = valueRaw;
        data.(name).maskMissingOrig = maskMissingOrig;
        data.(name).maskAnomaly    = maskAnomaly;
        data.(name).maskLongGap    = maskLongGap;
        data.(name).maskValid      = maskValid;
    end
end

% ======================================================================
% Local helper: clean a single time series
% ======================================================================
function [timeVec, valueClean, valueRaw, ...
          maskMissingOrig, maskAnomaly, maskLongGap, maskValid] = ...
          clean_single_series(timeIn, valuesIn, dtSeconds, preproc)
% clean_single_series
%   Clean one smart-meter time series using:
%       - time alignment to regular grid
%       - outlier detection (IQR and long zero runs)
%       - interpolation of short gaps
%       - optional smoothing

    % Ensure column vectors
    timeIn   = timeIn(:);
    valuesIn = double(valuesIn(:));

    % ---------- 1) Align on regular time grid ----------
    tt = timetable(timeIn, valuesIn, 'VariableNames', {'value'});

    % Create regular timetable with NaN for missing timestamps
    ttReg = retime(tt, 'regular', 'fillwithmissing', ...
                   'TimeStep', seconds(dtSeconds));

    timeVec   = ttReg.Time;
    valueRaw  = ttReg.value;           % may contain NaN where samples are missing

    % Original missing mask (before any cleaning)
    maskMissingOrig = isnan(valueRaw);

    % ---------- 2) Detect outliers and bad data ----------
    vals = valueRaw;
    maskValidValues = ~isnan(vals);

    % 2a) Negative values (not expected for energy/power)
    negMask = vals < 0 & maskValidValues;

    % 2b) High spikes based on IQR (robust to outliers)
    if any(maskValidValues)
        v = vals(maskValidValues);
        q1 = quantile(v, 0.25);
        q3 = quantile(v, 0.75);
        iqrVal = q3 - q1;

        highThr = q3 + preproc.iqrMultiplier * iqrVal;
        % Lower threshold: typically zero for energy
        lowThr = 0;

        spikeMask = (vals > highThr | vals < lowThr) & maskValidValues;
    else
        spikeMask = false(size(vals));
    end

    % 2c) Long runs of zeros (meter stuck / disconnected)
    zeroMask = (vals == 0) & maskValidValues;
    zeroRunThreshold = preproc.zeroRunThresholdSamples;

    zeroRunMask = false(size(vals));
    if any(zeroMask)
        d = diff([false; zeroMask(:); false]);
        runStarts = find(d == 1);
        runEnds   = find(d == -1) - 1;
        runLengths = runEnds - runStarts + 1;

        for i = 1:numel(runStarts)
            if runLengths(i) >= zeroRunThreshold
                zeroRunMask(runStarts(i):runEnds(i)) = true;
            end
        end
    end

    % Combine anomaly masks
    maskAnomaly = negMask | spikeMask | zeroRunMask;

    % Treat anomalies as missing for interpolation
    valsClean = vals;
    valsClean(maskAnomaly) = NaN;

    % ---------- 3) Interpolate short gaps ----------
    maxGap = preproc.maxInterpGapSamples;

    % Interpolate NaN segments up to maxGap samples
    valsInterp = fillmissing(valsClean, 'linear', 'MaxGap', maxGap);

    % Identify remaining NaNs as long gaps (too long to be trusted)
    maskLongGap = isnan(valsInterp);

    % For long gaps we keep NaN so that algorithms can ignore these times
    % when necessary. maskValid is everything that is not NaN after interp.
    maskValid = ~maskLongGap;

    % ---------- 4) Optional smoothing ----------
    win = preproc.smoothWindowSamples;
    if win > 1
        valsSmoothed = movmedian(valsInterp, win, ...
                                 'omitnan', 'Endpoints', 'shrink');
    else
        valsSmoothed = valsInterp;
    end

    valueClean = valsSmoothed;
end
