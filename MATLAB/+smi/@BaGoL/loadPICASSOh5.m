function SMD=loadPICASSOh5(DataDir,FileName)
%loadH5 Loads a dataset saved in H5-format by the PICASSO software.
% SMD=BaGoL.loadPICASSOh5(DataDir,FileName)
%
% Loads data from PICASSO h5 file and returns a SMD structure, which is
% readable by BaGoL.
% Note: PICASSO data may not be frame connected. 
%
% See Schnitzbauer, et al, "Super-resolution microscopy with DNA-PAINT",
% Nature Protocol, 2017.
%
% INPUTS:
%    DataDir:  The data directory.
%    FileName: Name of the file to be loaded.
%
% OUTPUT:
%    SMD:      SMD structure (see BaGoL properties)
%

% Created by:
%    Mohamadreza Fazel (LidkeLab 2019)

if ~strcmp('h5',FileName(end-1:end)) && ~strcmp('f5',FileName(end-1:end))
    error('File type has to be given at the end of file name.'); 
end
SM=h5read(fullfile(DataDir,FileName),'/locs');
SMD.X = SM.x;
SMD.Y = SM.y;
SMD.X_SE = SM.lpx;
SMD.Y_SE = SM.lpy;
if isfield(SM,'z')
    SMD.Z = SM.z;
    SMD.Z_SE = SM.lpz;
else
    SMD.Z = [];
    SMD.Z_SE = [];
end
SMD.FrameNum = SM.frame;

end
