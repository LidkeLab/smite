classdef Threshold < handle

% Threshold localizations based on various properties of the localizations.
% This is done by creating a ThreshFlag field for the SMD structure using the
% method setThreshFlag and the bounds in the SMF.MinMax (only those bounds
% defined will be applied).  applyThresh then applies the ThreshFlag to the
% localizations.  rejectedLocalizations will visually display what
% localizations were eliminated and for what reasons.  translateThreshFlag can
% also be used to produce human readable interpretations of the ThreshFlag.

% =============================================================================
properties

%   % MinMax - structure with each field that reqires thresholding
%   %          (same or subset of fields in SMD)
%   %          format:[Minvalue Maxvalue] for each field
%   %          { E.g.: MinMax.Bg = [BgMin BgMax] }
%   % NOTE: Only those fields provided will be used in thresholding.  The
%   %       default is no thresholding.  See Fields for possible thresholding
%   %       fields.  Typical fields used are:
%   %          X, Y, Z, X_SE, Y_SE, Z_SE, Photons, Bg, PSFSigma, PValue
%   MinMax=[]
    Verbose = 1;   % verbosity level

end % properties
% =============================================================================

% =============================================================================
properties(Constant = true)

   % Defining the possible fields for SMD thresholding.
   Fields={'X';'Y';'Z';'Photons';'Bg';'PSFSigma';'X_SE';'Y_SE';'Z_SE';'PValue'};
   %Fields={'X';'Y';'Z';'Photons';'Bg';'PSFSigma';'PSFSigmaX';'PSFSigmaY';   ...
   %        'X_SE';'Y_SE';'Z_SE';'Photons_SE';'Bg_SE';'PSFSigma_SE';         ...
   %        'PSFSigmaX_SE';'PSFSigmaY_SE';'ZOffset';'DatasetNum';'FrameNum'; ...
   %        'PValue';'LogLikelihood';'ConnectID'};

end % properties(Constant = true)
% =============================================================================

% =============================================================================
methods

   rejectedLocalizations(obj, SMD, options, SaveDir)
   [SMD, TFlag] = setThreshFlag(obj, SMD, MinMax)
   [ThreshFlagReadable, HotBits] = translateThreshFlagNS(obj, ThreshFlag)

end % methods
% =============================================================================

% =============================================================================
methods(Static)

   [SMR] = applyThresh(SMD, Verbose)
   MinMax = setMinMax(SMF)
   [ThreshFlagReadable, HotBits] = translateThreshFlag(ThreshFlag);
   success = unitTest()

end % methods(Static)
% =============================================================================

end % classdef Threshold
