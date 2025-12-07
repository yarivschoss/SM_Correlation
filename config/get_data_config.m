function cfg = get_data_config()
    % get_data_config – central place for all global parameters related
    % to data collection and simulation

    % Time parameters
    cfg.sampleTime_sec      = 900;     % Sampling interval: 15 minutes
    cfg.totalDuration_days  = 30;      % Total data duration in days

    % Network structure
    cfg.numTransformers     = 3;
    cfg.numCustomersTotal   = 120;
    cfg.maxCustomersPerTr   = 60;

    % Noise, missing data and measurement errors
    cfg.missingDataRate     = 0.02;    % 2% of samples are missing
    cfg.measurementNoiseStd = 0.01;    % Relative standard deviation of noise

    % Energy-related parameters
    cfg.nominalTransformerPower_kVA = 400;
    cfg.lossFactorApprox            = 0.03;   % Approximate technical losses (3%)

    % Energy-balance algorithm parameters
    cfg.energyBalanceTolerance      = 0.05;   % 5% allowed mismatch
    cfg.maxBalanceIterations        = 500;

    % Optimization parameters
    cfg.optim.learningRate          = 0.01;
    cfg.optim.maxEpochs             = 200;
    cfg.optim.stochasticJumpProb    = 0.1;    % Probability for random jump in SGD

    % E-CPC parameters
    cfg.ecpc.fftWindowSamples       = 96;     % Window size (e.g. one day of 15-min samples)
    cfg.ecpc.minEventAmplitudeRatio = 1e-3;   % Minimum event amplitude relative to total signal
    cfg.ecpc.correlationThreshold   = 0.8;

    % Paths for data and results
    cfg.paths.dataDir               = "data";
    cfg.paths.resultsDir            = "results";
end
