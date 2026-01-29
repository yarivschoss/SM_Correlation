function setup_paths()
    % setup_paths – add all relevant folders to the MATLAB path

    % thisFileDir  = folder where setup_paths.m is located (project_root)
    thisFileDir = fileparts(mfilename('fullpath'));

    % repoRoot = parent folder of project_root (SM_Correlation)
    repoRoot = fileparts(thisFileDir);

    % Add project_root itself (main.m, setup_paths.m)
    addpath(thisFileDir);

    % Add subfolders that sit next to project_root
    addpath(fullfile(repoRoot, "config"));
    addpath(fullfile(repoRoot, "algorithms"));
    addpath(fullfile(repoRoot, "utils"));
    addpath(fullfile(repoRoot, "data"));
    addpath(fullfile(repoRoot, "data_synthesis"));
    addpath(fullfile(repoRoot, "plot"));

    disp("Paths added successfully.");
end
