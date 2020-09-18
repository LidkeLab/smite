function [ThreshFlagReadable, HotBits] = ...
    translateThreshFlagNS(obj, ThreshFlag)
%translateThreshFlagNS translates ThreshFlag to a human readable string.
% This method will take the numeric array ThreshFlag (decimal integers) and
% output the human-readable threshold identifier given by the mapping in
% setThreshFlag().  The 'NS' suffix stands for non-static, as this is the
% non-static version of this method (there is a static version
% translateThreshFlag() which itself calls this method).
%
% NOTE: This method assumes that ThreshFlag is 32 bits, i.e.,
%       max possible ThreshFlag = 2^32 - 1, but it does not directly
%       attempt to typecast/validate the input type of ThreshFlag.
%
% EXAMPLE:
%   Given an SMD structure with SMD.ThreshFlag populated, this method can
%   be used as
%   ThreshFlagReadable = ...
%      smi_core.Threshold.translateThreshFlag(SMD.ThreshFlag);
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

% Define some parameters and initialize/pre-allocate arrays as needed.
NFlags = numel(ThreshFlag);
ThreshFlagReadable = cell(NFlags, 1);
HotBits = ThreshFlagReadable; % initialize
NBits = 32; % defined to emphasize this assumption

% Define our bit-to-string map such that ThreshFlagMap(bb) <-> threshold
% flag string for bit bb (note that this method will use the uncommon
% convention that bits start at 1, as is used in setThreshFlag()).
ThreshFlagMap = repelem("Unknown flag", NBits);
ThreshFlagFields = convertCharsToStrings(obj.Fields).';
ThreshFlagMap(1:numel(ThreshFlagFields)) = ThreshFlagFields;

% Loop through the input elements of ThreshFlag and translate.
for ii = 1:NFlags
    % Determine which bits of ThreshFlag(ii) are "hot", i.e., set to 1.
    HotBitBool = logical(bitget(ThreshFlag(ii), 1:NBits));
    HotBits{ii} = uint32(find(HotBitBool));

    % Populate the ii-th element of the ThreshFlagReadable array.
    ThreshFlagReadable{ii} = ThreshFlagMap(HotBitBool);
end

end % translateThreshFlagNS
