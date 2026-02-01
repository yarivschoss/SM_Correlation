function clear_synthetic_data(varargin)
% clear_synthetic_data
% Deletes all synthetic customer files and the transformer profile.
%
% Default project structure:
%   data/virtual_customers/*
%   data/transformer_profile.xlsx
%
% Usage:
%   clear_synthetic_data();
%   clear_synthetic_data('projectRoot', pwd);

% -------------------- Parse inputs --------------------
p = inputParser;
p.addParameter('projectRoot', "", @(x)isstring(x)||ischar(x));
p.parse(varargin{:});
opt = p.Results;

% -------------------- Resolve project root --------------------
if strlength(string(opt.projectRoot)) == 0
    thisFile = mfilename('fullpath');      % .../clear_synthetic_data.m
    thisDir  = fileparts(thisFile);
    projectRoot = string(fileparts(thisDir));
else
    projectRoot = string(opt.projectRoot);
end

% -------------------- Paths --------------------
customersDir   = fullfile(projectRoot, "data", "virtual_customers");
transformerXls = fullfile(projectRoot, "data", "transformer_profile.xlsx");

fprintf("Clearing synthetic data under:\n%s\n\n", projectRoot);

% -------------------- Delete customer files --------------------
if isfolder(customersDir)
    files = dir(customersDir);
    files = files(~[files.isdir]);  % only files

    for k = 1:numel(files)
        fpath = fullfile(files(k).folder, files(k).name);
        delete(fpath);
        fprintf("Deleted: %s\n", fpath);
    end
else
    warning("Folder not found: %s", customersDir);
end

% -------------------- Delete transformer profile --------------------
if isfile(transformerXls)
    delete(transformerXls);
    fprintf("Deleted: %s\n", transformerXls);
else
    fprintf("Transformer profile not found: %s\n", transformerXls);
end

fprintf("\n✔ Synthetic data cleared successfully.\n");
end
