function setup_paths()
    % setup_paths – add all relevant subfolders to the MATLAB path
    rootDir = fileparts(mfilename('fullpath'));
    addpath(fullfile(rootDir, "config"));
    addpath(fullfile(rootDir, "algorithms"));
    addpath(fullfile(rootDir, "utils"));
    addpath(fullfile(rootDir, "data"));
end
