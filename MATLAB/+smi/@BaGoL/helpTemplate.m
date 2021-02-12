function [Output1,Output2] = helpTemplate(Input1,Input2,Input3)
%helpTemplate One line summary that appears in class method summary 
%   [Output1,Output2] = helpTemplate(Input1,Input2,Input3)
%
%   Detailed explanation goes here. Yes, detailed. All help 
%   must have the same format. This tempate describes what is required in
%   the help and also defines a formatting convention.  Check the
%   formatting of your function using 'doc SMA_Core', 
%   'help SMA_Core.helpTemplate' and 'doc SMA_Core.helpTemplate'
%
%   Use blank lines to make help more readable. Help should go as a 
%   contigous block below the 'function' command. Help text in the editor
%   must not go beyond the right margin line ->
%   
%   When needed, Inputs should have addition information in () in the order
%   of: (Units) (Size) (Default = ...). Inputs without (Default = ..) are \
%   required. Inputs with (Default = ..) must be set to this value by the 
%   function if input is not given. 
%   
%   If inputs/outputs can't be adequately described one line, more
%   information should be included in this section of the help. 
%
% INPUTS:
%   Input1:     Short Description (Pixels)
%   Input2:     Short Description (Default=100)
%   Input3:     Short Description (microns) (Nx1)
%
% OUTPUTS:
%   Output1:    Short Description
%   OutPut2:    Short Description. A structure with fields:
%       A:      Short Description
%       B:      Short Description
%
% REQUIRES:
%   Statistics Toolbox
%   Parallel Procesing Toolbox
%   NVidia GPU
%
% CITATION:
%   Full citatation with all authors. (Citation for a paper if applicable)        

%created by 
% Full Name (Lidkelab, 2018)

doc SMA_Core
help SMA_Core.helpTemplate
doc SMA_Core.helpTemplate

end

