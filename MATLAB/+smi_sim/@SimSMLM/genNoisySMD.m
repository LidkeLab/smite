function [SMD_Data] = genNoisySMD(obj,SMD_Model)

% This function takes the SMD_Model struct as an input, generates the noisy 
% SMD struct and return SMD_Data with the following fields: 

% OUTPUTS:

% SMD_Data.X
% SMD_Data.Y
% SMD_Data.X_SE
% SMD_Data.Y_SE
% SMD_Data.PSFSigma
% SMD_Data.FrameNum
% SMD_Data.Photons
% SMD_Data.Bg

% Copy all fields from SMD_Model into SMD_Data.
if nargin<2
    SMD_Model=obj.SMD_Model;
end

SMD_Data = SMD_Model;

SMD_Data.X_SE=(obj.PSFSigma)./sqrt(SMD_Model.Photons);
SMD_Data.X=SMD_Data.X+randn(size(SMD_Data.X_SE)).*SMD_Data.X_SE;

SMD_Data.Y_SE=(obj.PSFSigma)./sqrt(SMD_Model.Photons);
SMD_Data.Y=SMD_Data.Y+randn(size(SMD_Data.Y_SE)).*SMD_Data.Y_SE;

%SMD_Data.Photons=SMD_Model.Photons;

SMD_Data.PSFSigma=SMD_Model.PSFSigma;

%SMD_Data.FrameNum=SMD_Model.FrameNum;

SMD_Data.Bg=obj.Bg.*ones([length(SMD_Model.Photons),1]);

end


