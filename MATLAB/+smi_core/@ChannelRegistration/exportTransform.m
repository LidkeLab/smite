function exportTransform(obj, FileDir, FileName)
%exportTransform exports transform information into a .mat file.
% This method will save a bunch of relevant class fields into a .mat file
% in the specified location.
%
% INPUTS:
%   FileDir: Directory in which we'll save the transform information.
%            (char array/string)(Default = pwd())
%   FileName: Name of the .mat file we'll save. 
%             (char array/string)
%             (Default = 'RegistrationTransform' plus a string containing 
%             time)

% Created by:
%   David J. Schodt (Lidke Lab, 2020)


% Set defaults if needed.
if (~exist('FileDir', 'var') || isempty(FileDir))
    FileDir = pwd();
end
if (~exist('FileName', 'var') || isempty(FileName))
    CurrentTime = cellstr(num2str(round(clock().')));
    DateString = erase(strjoin(CurrentTime, '_'), ' ');
    FileName = ['RegistrationTransform', '_', DateString, '.mat'];
end

% Save the registration information.
RegistrationTransform = obj.RegistrationTransform;
Coordinates = obj.Coordinates;
NNeighborPoints = obj.NNeighborPoints;
PolynomialDegree = obj.PolynomialDegree;
save(fullfile(FileDir, FileName), 'RegistrationTransform', ...
    'Coordinates', 'NNeighborPoints', 'PolynomialDegree')


end