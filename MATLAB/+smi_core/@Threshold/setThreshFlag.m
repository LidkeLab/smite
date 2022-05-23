function [SMD, TFlag] = setThreshFlag(obj, SMD, MinMax)
% Creates ThreshFlag field for SMD, the same size as SMD.X .
% Values in the ThreshFlag are decimal equivalent of a 32 bit binary number.
% Each bit indicates whether the value for a field falls within the range
% provided by MinMax (success=0) or not (failure=1). A localization which
% succeeds to pass through all the threshold conditions will have a
% ThreshFlag value of '0'.  Note that if a field is not present in the MinMax
% structure, the corresponding bit in ThreshFlag will be set to '0'.
%
% Format for the bits (flags) in the binary number (starting from left) is:
%     bit 11   'LogL'
%     bit 10   'PValue'
%     bit  9   'Z_SE'
%     bit  8   'Y_SE'
%     bit  7   'X_SE'
%     bit  6   'PSFSigma'
%     bit  5   'Bg'
%     bit  4   'Photons'
%     bit  3   'Z'
%     bit  2   'Y'
%     bit  1   'X'
% x marks fields that do not seem to be used elsewhere.  Corresponding tests
% for these fields below has been removed in order to simplify the code.
%
% Example: Value of SMD.Threshflag '536' indicates a binary value of
% '00000000000000000000001000011000' (from dec2bin(536,32)) - showing the
% particular localisation failed the threshold conditions for fields
% 'Photons' (bit 4), 'Bg' (bit 5) and 'PValue' (bit 10).
%
% INPUTS:
%    SMD      with fields as defined (or a subset)
%    MinMax   structure with each field that reqires thresholding
%             (same or subset of fields in SMD)
%             format:[Minvalue Maxvalue] for each field
%             { E.g.: MinMax.Bg = [BgMin BgMax] }
%    NOTE: Only those fields provided and non-empty will be used in
%          thresholding.  The default is no thresholding.  See Fields for
%          possible thresholding fields.  Typical fields used are:
%             X, Y, Z, X_SE, Y_SE, Z_SE, Photons, Bg, PSFSigma, PValue
%
% OUTPUTS:
%    SMD      the same as input, but with the field 'ThreshFlag' added/appended
%    TFlag    a copy of the SMD.ThreshFlag
%
% REQUIRES:
%    Statistics Toolbox

%Created by
%   Sandeep Pallikkuth, Lidke Lab 2017
%Revised by Hanieh Mazloom-Farsibaf (Lidke Lab, 2018)
%Revised by Michael Wester, 2020

    if nargin<1
        error('Please input structure SMD\n')
    end

    % Checking MinMax has fields that are subset of those of SMD
    SMDFN=fieldnames(SMD);
    MMFN=fieldnames(MinMax);
    if isempty(SMDFN)>0
        error('SMD structure empty. Please add fields.\n');
    end
    if isempty(MMFN)>0
        error('MinMax structure empty. Please add fields to threshold.\n')
    end
    Check1=ismember(MMFN,obj.Fields);
    if any(Check1(:,:)==0)
        error('Fields in MinMax should be a subset of those in SMD');
    end

    % Case of no data!
    if length(SMD.X) == 0
       warning('No localizations!');
       TFlag = [];
       SMD.ThreshFlag=TFlag;
       return;
    end
     
    % checking if SMD.ThreshFlag exists
    if sum(strcmp(fieldnames(SMD), 'ThreshFlag')) == 1
        if isempty(SMD.ThreshFlag)
            TFlag(length(SMD.X),1)=uint32(0);
        else
            TFlag=SMD.ThreshFlag;
            % if ThreshFlag exists, then initialize the bits corresponding to
            % fields in MinMax to 0
            for i=1:length(MMFN)
                fn1=MMFN(i);
                xx=strmatch(fn1,obj.Fields,'exact');
                TFlag(:,:)=bitset(TFlag(:,:),xx,0);
            end
        end
    else
        TFlag(length(SMD.X),1)=uint32(0);
    end

    % checking for MinMax field values and assigning failure flags (value 1)
    % for corresponding SMD fields
    if sum(strcmp(fieldnames(MinMax), 'X')) == 1
        if ~isempty(MinMax.X)==1 && ~isempty(SMD.X)
            tflagX=SMD.X<MinMax.X(1)|SMD.X>MinMax.X(2)|~isreal(SMD.X)| isnan(SMD.X);
            k=find(tflagX); % find the position of failures in the array
            xx=strmatch('X',obj.Fields,'exact');
            TFlag(k)=bitset(TFlag(k),xx);
        end
    end
    if sum(strcmp(fieldnames(MinMax), 'Y')) == 1
        if ~isempty(MinMax.Y)==1 && ~isempty(SMD.Y)
            tflagY=SMD.Y<MinMax.Y(1)|SMD.Y>MinMax.Y(2)|~isreal(SMD.Y)|isnan(SMD.Y);
            k=find(tflagY); % find the position of failures in the array
            xx=strmatch('Y',obj.Fields,'exact');
            TFlag(k)=bitset(TFlag(k),xx);
        end
    end
    if sum(strcmp(fieldnames(MinMax), 'Z')) == 1
        if ~isempty(MinMax.Z)==1 && ~isempty(SMD.Z)
            tflagZ=SMD.Z<MinMax.Z(1)|SMD.Z>MinMax.Z(2)|~isreal(SMD.Z)|isnan(SMD.Z);
            k=find(tflagZ); % find the position of failures in the array
            xx=strmatch('Z',obj.Fields,'exact');
            TFlag(k)=bitset(TFlag(k),xx);
        end
    end
    if sum(strcmp(fieldnames(MinMax), 'Photons')) == 1
        if ~isempty(MinMax.Photons)==1 && ~isempty(SMD.Photons)
            tflagI=SMD.Photons<MinMax.Photons(1)|SMD.Photons>MinMax.Photons(2)...
                |~isreal(SMD.Photons)|isnan(SMD.Photons);
            k=find(tflagI); % find the position of failures in the array
            xx=strmatch('Photons',obj.Fields,'exact');
            TFlag(k)=bitset(TFlag(k),xx);
        end
    end
    if sum(strcmp(fieldnames(MinMax), 'Bg')) == 1
        if ~isempty(MinMax.Bg)==1 && ~isempty(SMD.Bg)
            tflagBg=SMD.Bg<MinMax.Bg(1)|SMD.Bg>MinMax.Bg(2)|~isreal(SMD.Bg)...
                |isnan(SMD.Bg);
            k=find(tflagBg); % find the position of failures in the array
            xx=strmatch('Bg',obj.Fields,'exact');
            TFlag(k)=bitset(TFlag(k),xx);
        end
    end
    if sum(strcmp(fieldnames(MinMax), 'PSFSigma')) == 1
        if ~isempty(MinMax.PSFSigma)==1 && ~isempty(SMD.PSFSigma)
            tflagPSFSigma=SMD.PSFSigma<MinMax.PSFSigma(1)|...
                SMD.PSFSigma>MinMax.PSFSigma(2)|~isreal(SMD.PSFSigma)...
                |isnan(SMD.PSFSigma);
            k=find(tflagPSFSigma); % find the position of failures in the array
            xx=strmatch('PSFSigma',obj.Fields,'exact');
            TFlag(k)=bitset(TFlag(k),xx);
        end
    end
    if sum(strcmp(fieldnames(MinMax), 'X_SE')) == 1
        if ~isempty(MinMax.X_SE)==1 && ~isempty(SMD.X_SE)
            tflagXSE=SMD.X_SE<MinMax.X_SE(1)|SMD.X_SE>MinMax.X_SE(2)|...
                SMD.X_SE~=real(SMD.X_SE)|isnan(SMD.X_SE);

            k=find(tflagXSE); % find the position of failures in the array
            xx=strmatch('X_SE',obj.Fields,'exact');
            TFlag(k)=bitset(TFlag(k),xx);
        end
    end
    if sum(strcmp(fieldnames(MinMax), 'Y_SE')) == 1
        if ~isempty(MinMax.Y_SE)==1 && ~isempty(SMD.Y_SE)
            tflagYSE=SMD.Y_SE<MinMax.Y_SE(1)|SMD.Y_SE>MinMax.Y_SE(2)|...
                SMD.Y_SE~=real(SMD.Y_SE)|isnan(SMD.Y_SE);
            k=find(tflagYSE); % find the position of failures in the array
            xx=strmatch('Y_SE',obj.Fields,'exact');
            TFlag(k)=bitset(TFlag(k),xx);
        end
    end
    if sum(strcmp(fieldnames(MinMax), 'Z_SE')) == 1
        if ~isempty(MinMax.Z_SE)==1 && ~isempty(SMD.Z_SE)
            tflagZSE=SMD.Z_SE<MinMax.Z_SE(1)|SMD.Z_SE>MinMax.Z_SE(2)|...
                SMD.Z_SE~=real(SMD.Z_SE)|isnan(SMD.Z_SE);
            k=find(tflagZSE); % find the position of failures in the array
            xx=strmatch('Z_SE',obj.Fields,'exact');
            TFlag(k)=bitset(TFlag(k),xx);
        end
    end
    if sum(strcmp(fieldnames(MinMax), 'PValue')) == 1
        if ~isempty(MinMax.PValue)==1 && ~isempty(SMD.PValue)
            tflagPValue=SMD.PValue<MinMax.PValue(1)|...
                SMD.PValue>MinMax.PValue(2)|~isreal(SMD.PValue)...
                |isnan(SMD.PValue);
            k=find(tflagPValue); % find the position of failures in the array
            xx=strmatch('PValue',obj.Fields,'exact');
            TFlag(k)=bitset(TFlag(k),xx);
        end
    end
    if sum(strcmp(fieldnames(MinMax), 'LogLikelihood')) == 1
        if ~isempty(MinMax.LogLikelihood)==1 && ~isempty(SMD.LogLikelihood)
            tflagLogL=SMD.LogLikelihood<MinMax.LogLikelihood(1)|...
                SMD.LogLikelihood>MinMax.LogLikelihood(2)|~isreal(SMD.LogLikelihood)...
                |isnan(SMD.LogLikelihood);
            k=find(tflagLogL); % find the position of failures in the array
            xx=strmatch('LogLikelihood',obj.Fields,'exact');
            TFlag(k)=bitset(TFlag(k),xx);
        end
    end

    % adding/appending ThresholdFlag field of SMD
    SMD.ThreshFlag=TFlag;

end % setThreshFlag
