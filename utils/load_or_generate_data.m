function rawData = load_or_generate_data(cfg)
    % load_or_generate_data – load measurement files if available,
    % otherwise generate synthetic data.

    % Example skeleton:
    % 1) Try loading an existing .mat file
    matFile = fullfile(cfg.paths.dataDir, "rawData.mat"); %% rawData is a temp name - need to change!
    if isfile(matFile)
        s = load(matFile);
        rawData = s.rawData;
        return;
    end

    % 2) If file does not exist – create simple synthetic data as a placeholder
    numSamples    = cfg.totalDuration_days * 24 * 3600 / cfg.sampleTime_sec;
    T             = numSamples;
    Tr            = cfg.numTransformers;
    C             = cfg.numCustomersTotal;

    % TODO: replace random data with a more realistic consumption model
    rawData.trPower   = rand(T, Tr);    % Transformer power/energy time-series
    rawData.custPower = rand(T, C);     % Customer power/energy time-series
end
