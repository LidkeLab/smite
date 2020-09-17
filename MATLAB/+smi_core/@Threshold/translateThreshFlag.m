function [ThreshFlagReadable, HotBits] = translateThreshFlag(ThreshFlag)
%translateThreshFlag is a static "wrapper" for translateThreshFlagNS
% This method is a "wrapper" for the method
% smi_core.Threshold.translateThreshFlagNS(). The purpose of this static
% method is to be a more user-friendly version of the non-static method
% translateThreshFlagNS(). This method will create an instance of the
% Threshold class solely for the purpose of calling
% translateThreshFlagNS().
%
% NOTE: This method assumes that ThreshFlag is 32 bits, i.e.,
%       max possible ThreshFlag = 2^32 - 1, but it does not directly
%       attempt to typecast/validate the input type of ThreshFlag.
%
% EXAMPLE:
%   Given an SMD structure with SMD.ThreshFlag populated, this method can
%   be used as
%   ThreshFlagReadable = smi_core.Threshold.translateThreshFlag(...
%       SMD.ThreshFlag);
%
% INPUTS:
%   ThreshFlag: Array containing the decimal integer threshold flags as
%               described in setThreshFlag().
%               (unsigned integer, up to uint32 array)
%
% OUTPUTS:
%   ThreshFlagReadable: Cell array containing the translated
%                       (human-readable) identifier corresponding to
%                       ThreshFlag.  There is a direct mapping between the
%                       index of ThreshFlag and ThreshFlagReadable, i.e.
%                       ThreshFlagReadable{ii} <-> ThreshFlag(ii)
%                       (numel(ThreshFlag) x 1 cell array of strings)
%   HotBits: Cell array containing the "hot" bits corresponding to each
%            ThreshFlag (e.g., instead of the string "PValue", HotBits
%            would have the integer 20, corresponding to bit 20 being hot).
%            (numel(ThreshFlag) x 1 cell array of uint64)
%
% REQUIRES:
%   MATLAB R2016b to use the datatype 'string'

% Created by:
%   David James Schodt (Lidkelab, 2020)


% Create an instance of the smi_core.Threshold class.
Threshold = smi_core.Threshold();

% Call the non-static class method translateThreshFlagNS().
[ThreshFlagReadable, HotBits] = ...
    Threshold.translateThreshFlagNS(ThreshFlag);


end % translateThreshFlag