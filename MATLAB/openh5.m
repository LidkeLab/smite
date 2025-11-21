function [H5Structure] = openh5(FilePath)
%Extracts contents of an h5 file into H5Structure.
% This method will extract the Data and Attributes from a .h5 file specified by
% FilePath.  It is called when open or uiopen is invoked on a .h5 file.  For
% example, uiopen is invoked when a file is dragged and dropped into the MATLAB
% command window.  See REQUIRES for MATLAB version oddities.
%
% Examples:
%   H5Structure = readH5File('C:\file.h5') will extract all
%       contents of file.h5 and store them in H5Structure.
%
% INPUTS:
%   FilePath: String containing the path to the .h5 file of interest.
%
% OUTPUTS:
%   H5Structure: Structured array containing the information extracted from
%                the .h5 file at FilePath.
%
% REQUIRES:
%   MATLAB 2016b or later (this function no longer works alone in 2024a, but
%                          requires a modified uiopen.m that shadows the
%                          built-in version as well---see uiopen.m in this
%                          directory; note that the modified uiopen.m also
%                          works in MATLABs earlier than 2024a also, so this
%                          remedy is general, but if anything weird happens
%                          when dragging and dropping, this could easily be
%                          due to putting the supplied uiopen.m on MATLAB's
%                          path)
 
% Created by;
%     Michael J. Wester and Claude (2025)

fprintf('Opening h5 with custom loader: %s\n', FilePath);
%  GroupName: (optional) Name of a specific group in the .h5 file to be
%             extracted.
GroupName = inputdlg('Specific group in the .h5 file to be extracted.');
if isempty(GroupName)
   disp('Canceled!');
   return
end
GroupName = GroupName{1};
if ~isempty(GroupName)
   [H5Structure] = smi_core.LoadData.readH5File(FilePath, GroupName);
else
   [H5Structure] = smi_core.LoadData.readH5File(FilePath);
end
assignin('base', 'h5Struct', H5Structure);
disp('Loaded to variable: h5Struct');

end
