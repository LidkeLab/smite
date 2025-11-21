%uiopen to shadow MATLAB's built-in functionality, allowing dragging and
% dropping of SMITE .h5 files to the command window to work naturally, while
% retaining MATLAB's functionality for other file types.  See also the
% documentation for openh5.m
%
% Make sure this uiopen.m is on your path before MATLAB's built-in version.
% Check with:
%
%   matlabwhich uiopen -all

function uiopen(type, direct)
% Custom uiopen that intercepts H5 files only

% Check if this is an H5 file being dragged
if nargin >= 2 && direct == 1
    [~, ~, ext] = fileparts(type);
    if ~(strcmpi(ext, '.h5') || strcmpi(ext, '.hdf5'))
        % Not an H5 file - pass to original MATLAB uiopen
        callOriginalUiopen(type, direct);
        return;
    end
    % It's H5 - handle it
    autoOpen(type);
    return;
end

% For other cases (no arguments or type specified), handle H5 dialog
if nargin == 0
    [file, path] = uigetfile('*.h5;*.hdf5', 'Select H5 file');
    if isequal(file, 0), return; end
    autoOpen(fullfile(path, file));
else
    callOriginalUiopen(type);
end
end

function callOriginalUiopen(varargin)
    % Temporarily remove our uiopen from path and call MATLAB's
    thisFile = which('uiopen');
    thisPath = fileparts(thisFile);
    rmpath(thisPath);
    try
        uiopen(varargin{:});
    catch ME
        addpath(thisPath);
        rethrow(ME);
    end
    addpath(thisPath);
end

function autoOpen(filepath)
    [~, ~, ext] = fileparts(filepath);
    ext = ext(2:end);
    
    funcName = ['open' ext];
    
    if exist(funcName, 'file')
        fprintf('Opening with %s: %s\n', funcName, filepath);
        feval(funcName, filepath);
    else
        fprintf('Loading H5 file: %s\n', filepath);
        try
            info = h5info(filepath);
            [~, name, ~] = fileparts(filepath);
            varname = matlab.lang.makeValidName(name);
            assignin('base', varname, info);
            fprintf('Assigned to workspace variable: %s\n', varname);
        catch ME
            fprintf('Error loading H5: %s\n', ME.message);
        end
    end
end
