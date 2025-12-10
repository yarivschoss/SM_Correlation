function rawData = load_or_generate_data(cfg)
% load_or_generate_data
%   If Excel files exist in the data folder, read them using read_data()
%   and build a struct of structs.
%   Otherwise, generate synthetic data as a fallback.

    dataDir = cfg.paths.dataDir;
    xlsFiles = dir(fullfile(dataDir, '*.xlsx'));

    if ~isempty(xlsFiles)
        % ----- Real data path -----
        rawData = struct();

        for k = 1:numel(xlsFiles)
            filePath = fullfile(xlsFiles(k).folder, xlsFiles(k).name);

            % Use filename (without extension) as struct field
            [~, baseName, ~] = fileparts(xlsFiles(k).name);

            % Example: rawData.Transformer1, rawData.YardA, etc.
            rawData.(baseName) = read_data(filePath);
        end

        return;
    end

    % ----- Synthetic data fallback  -----
    numSamples = cfg.totalDuration_days * 24 * 3600 / cfg.sampleTime_sec;
    T          = numSamples;
    Tr         = cfg.numTransformers;
    C          = cfg.numCustomersTotal;

    rawData.trPower   = rand(T, Tr);    % placeholder transformer series
    rawData.custPower = rand(T, C);     % placeholder customer series
end
