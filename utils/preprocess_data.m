function data = preprocess_data(rawData, cfg)
    % preprocess_data – clean and transform raw data:
    % add noise, drop samples, normalize, align time axis, etc.

    data = rawData;

    % TODO:
    %   - Add measurement noise according to cfg.measurementNoiseStd
    %   - Remove a fraction of samples according to cfg.missingDataRate
    %   - Apply basic filters / time alignment if needed
end
