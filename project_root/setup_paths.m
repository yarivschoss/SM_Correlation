function setup_paths()
    % setup_paths – add all relevant subfolders to the MATLAB path

    projectRoot = fileparts(mfilename('fullpath'));
    
    addpath(fullfile(projectRoot, "project_root"));
    addpath(fullfile(projectRoot, "config"));
    addpath(fullfile(projectRoot, "algorithms"));
    addpath(fullfile(projectRoot, "utils"));
    addpath(fullfile(projectRoot, "data"));

    disp("Paths added successfully.");
end
