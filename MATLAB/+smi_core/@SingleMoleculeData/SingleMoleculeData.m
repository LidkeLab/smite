classdef SingleMoleculeData
    % SingleMoleculeData A class defining the Single Molecule Data structure
    %
    % This datatype is one of the primary results structures in the smite
    % enviroment. The SMD structure is an input and output of many smi
    % methods. It intended to be extensible.
    % The SMD class implements tools for working with SMD structures,
    % but the data strcuture itself is not an object of the class.
    %
    % The structure has the following Properties:
    %
    % SMD:
    %   NDims:          Number of dimensions in localization information (2,3)
    %   NFrames:        Number of image frames in raw data sets
    %   NDatasets:      Number of 3D image stacks
    %   XSize:          Number of pixels in X dimension of raw data
    %   YSize:          Number of pixels in Y dimension of raw data
    %   XBoxCorner:     X coordinate of top right box corner
    %   YBoxCorner:     Y coordinate of top right box corner
    %   ZOffset:        Z position of focal plane of sequence
    %   X:              Estimated X position
    %   Y:              Estimated Y position
    %   Z:              Estimated Z position
    %   Photons:        Estimated Photons  (Integrated collected photons)
    %   Bg:             Estimated Background (Photons/Pixel)
    %   PSFSigma:       Estimated or Fixed Sigma of 2D Gaussian PSF Model (symmetric PSF)
    %   PSFSigmaX:      Estimated or FixedX Sigma of 2D Gaussian PSF Model (asymmetric PSF)
    %   PSFSigmaY:      Estimated or FixedY Sigma of 2D Gaussian PSF Model (asymmetric PSF)
    %   X_SE:           Standard Error of X
    %   Y_SE:           Standard Error of Y
    %   Z_SE:           Standard Error of Z
    %   Photons_SE:     Standard Error of Photons
    %   Bg_SE:          Standard Error of Bg
    %   PSFSigma_SE:    Standard Error of PSFSigma
    %   PSFSigmaX_SE:   Standard Error of PSFSigmaX
    %   PSFSigmaY_SE:   Standard Error of PSFSigmaY
    %   DatasetNum:     File number from which localization originates
    %   FrameNum:       Frame number from which localization originates
    %   PValue:         p-value of fit
    %   LogLikelihood:  Log likelihood of fit
    %   ConnectID:      Identifies the same emitter accross multiple frames
    %   ThreshFlag:     Indicates a valid fit.  0=valid.  See SMA_Core.ThresholdSM
    %   DriftX:         X drift relative to first frame (Pixels) (NFrames x NDatasets)
    %   DriftY:         Y drift relative to first frame (Pixels) (NFrames x NDatasets)
    %   DriftZ:         Z drift relative to first frame (Pixels) (NFrames x NDatasets)
    %
    %
    % SEE ALSO:
    %   smi_core.SMF, smi_core.TR
    %
    
    properties
        
    end
    
    methods (Static)
        function [SMD] = createSMD()
            %createSMD Creates an empty Single Molecule Data (SMD) structure
            SMD.NDims=[];
            SMD.NFrames=[];
            SMD.NDatasets=[];
            SMD.XSize=[];
            SMD.YSize=[];
            SMD.XBoxCorner=[];
            SMD.YBoxCorner=[];
            SMD.ZOffset=[];
            SMD.X=[];
            SMD.Y=[];
            SMD.Z=[];
            SMD.Photons=[];
            SMD.Bg=[];
            SMD.PSFSigma=[];
            SMD.PSFSigmaX=[];
            SMD.PSFSigmaY=[];
            SMD.X_SE=[];
            SMD.Y_SE=[];
            SMD.Z_SE=[];
            SMD.Photons_SE=[];
            SMD.Bg_SE=[];
            SMD.PSFSigma_SE=[];
            SMD.PSFSigmaX_SE=[];
            SMD.PSFSigmaY_SE=[];
            SMD.DatasetNum=[];
            SMD.FrameNum=[];
            SMD.PValue=[];
            SMD.LogLikelihood=[];
            SMD.ConnectID=[];
            SMD.ThreshFlag=[];
            SMD.DriftX=[];
            SMD.DriftY=[];
            SMD.DriftZ=[];
        end
        
        function [SMD] = catSMD(SMD1, SMD2)
            %This method concatenates two SMD structures into one.
            
            % Define which SMD fields are vector fields vs. scalar fields.
            % NOTE: We could just check which fields are numeric vectors
            %       with (numel() > 1), but that convenience limits us to
            %       SMD structures with more than one localization.  Fields
            %       that aren't present in this list will still be checked
            %       in that manner just to be safe.
            VectorFields = {'XBoxCorner', 'YBoxCorner', ...
                'X', 'Y', 'Z', 'X_SE', 'Y_SE', 'Z_SE', ...
                'Photons', 'Photons_SE', 'Bg', 'Bg_SE', ...
                'PSFSigma', 'PSFSigmaX', 'PSFSigmaY', ...
                'PSFSigma_SE', 'PSFSigmaX_SE', 'PSFSigmaY_SE', ...
                'PValue', 'LogL', 'ThreshFlag', ...
                'DatasetNum', 'FrameNum', 'ConnectID', ...
                'DriftX', 'DriftY', 'DriftZ'};
            
            % Create a list of fields present in each of SMD1 and SMD2.
            FieldNames1 = fieldnames(SMD1);
            FieldNames2 = fieldnames(SMD2);
            
            % Check if SMD1 and SMD2 have the same fields, warning the user
            % if this isn't the case.
            UniqueFields1 = FieldNames1(...
                ~ismember(FieldNames1, FieldNames2));
            UniqueFields2 = FieldNames2(...
                ~ismember(FieldNames2, FieldNames1));
            if ~isempty(UniqueFields1)
                PrintFriendlyFields1 = cell2mat(cellfun(@(X) [X, ' '], ...
                    UniqueFields1, 'UniformOutput', false).');
                warning(['smi_core.SingleMoleculeData.catSMD() - ', ...
                    'The following SMD1 fields aren''t ', ...
                    'present in SMD2: ', PrintFriendlyFields1])
            end
            if ~isempty(UniqueFields2)
                PrintFriendlyFields2 = cell2mat(cellfun(@(X) [X, ' '], ...
                    UniqueFields2, 'UniformOutput', false).');
                warning(['smi_core.SingleMoleculeData.catSMD() - ', ...
                    'The following SMD2 fields aren''t ', ...
                    'present in SMD1: ', PrintFriendlyFields2])
            end
            
            % Place the unique fields (fields that are only in one of SMD1 
            % or SMD2) in the SMD structure directly.
            SMD = smi_core.SingleMoleculeData.createSMD();
            for uu = 1:numel(UniqueFields1)
                SMD.(UniqueFields1{uu}) = SMD1.(UniqueFields1{uu});
            end
            for uu = 1:numel(UniqueFields2)
                SMD.(UniqueFields2{uu}) = SMD2.(UniqueFields2{uu});
            end
            
            % Loop through fields present in SMD1 and SMD2 and begin 
            % constructing the concatenated output SMD.
            AllFieldNames = unique([FieldNames1; FieldNames2]);
            SharedFields = AllFieldNames(...
                ismember(AllFieldNames, FieldNames1) ...
                & ismember(AllFieldNames, FieldNames2));
            IsVectorField = ismember(SharedFields, VectorFields);
            for ff = 1:numel(SharedFields)
                % Place the current field in the output SMD, concatenating
                % when appropriate.
                if IsVectorField(ff)
                    % Vector fields should be concatenated directly.
                    SMD.(SharedFields{ff}) = [SMD1.(SharedFields{ff}); ...
                        SMD2.(SharedFields{ff})];
                elseif isequal(...
                        SMD1.(SharedFields{ff}), SMD2.(SharedFields{ff}))
                    % We want to ensure non-vector fields are consistent
                    % with one another before setting the output field.
                    SMD.(SharedFields{ff}) = SMD1.(SharedFields{ff});
                else
                    % This appears to be a non-vector field, but the field 
                    % is not equal in SMD1 and SMD2; we will trust the user 
                    % and concatenate these two as though it were a vector 
                    % field.
                    warning(['smi_core.SingleMoleculeData.catSMD() - ', ...
                        'SMD field ''%s'' is different in each of ', ...
                        'SMD1 and SMD2, attempting to concatenate.'], ...
                        SharedFields{ff})
                    try
                        SMD.(SharedFields{ff}) = ...
                            [SMD1.(SharedFields{ff}); ...
                            SMD2.(SharedFields{ff})];
                    catch Exception
                        warning(...
                            ['smi_core.SingleMoleculeData.catSMD() - ', ...
                            'Concatenation of field ''%s'' failed, ', ...
                            'storing as cell array.\n', ...
                            '%s, %s'], SharedFields{ff}, ...
                            Exception.identifier, Exception.message)
                        SMD.(SharedFields{ff}) = ...
                            {SMD1.(SharedFields{ff}); ...
                            SMD2.(SharedFields{ff})};
                    end
                end
            end
            
            % Revise fields that may or may not be vector fields (e.g.,
            % PSFSigma is either scalar or vector depending on the fit that
            % was used).
            if (numel(unique(SMD.PSFSigma)) == 1)
                SMD.PSFSigma = unique(SMD.PSFSigma);
                SMD.PSFSigma_SE = unique(SMD.PSFSigma_SE);
            end
            if (numel(unique(SMD.DatasetNum)) == 1)
                SMD.DatasetNum = unique(SMD.DatasetNum);
            end
        end
        
    end
end