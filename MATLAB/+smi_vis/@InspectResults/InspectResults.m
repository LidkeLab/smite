classdef InspectResults < handle
    %InspectResults contains methods useful for inspecting SR data.
    % This class provides an interface for inspecting super-resolution
    % results.  The primary intention is that this class can be used to
    % associate localizations from a reconstruction image (e.g., a Gaussian
    % SR image) and the entries of a Single Molecule Data structure
    % (see smi_core.SingleMoleculeData).
    %
    % The suggested usage of this class is through the GUI, which will be
    % opened by default when creating an instance of this class:
    %   Inspect = smi_vis.InspectResults();
    % From the GUI, you should load a super-resolution image reconstruction
    % (saved as, e.g., a .png) as well as a *_Results.mat file containing
    % the SMD structure (see usage of smi.SMLM, which generates a 
    % *_Results.mat file).  Once the image is loaded, it'll be displayed in
    % the GUI figure.  Hovering over the GUI figure will reveal a toolbar
    % with a '+' icon which can be used to highlight a ROI and display
    % info. about localizations from the ROI.
    %
    % REQUIRES:
    %   Image Processing Toolbox
    %   Statistics and Machine Learning Toolbox

    % Created by:
    %   David J. Schodt (Lidke lab, 2021)
    
    properties
        % Single Molecule Data structure. (see smi_core.SingleMoleculeData)
        SMD

        % Super-resolution reconstruction image. (YxXx3 float)
        SRImage

        % Figure containing the interactive GUI.
        GUIFigure

        % Axes containing the SR image.
        ImageAxes

        % ROI in SMD coordinates. ([YStart, XStart, YEnd, XEnd])
        ROI
    end

    properties (Dependent)
        % Isolated SMD within the ROI defined by ROIHandle.
        SMDIsolated
    end


    properties (Hidden, SetAccess = 'protected')
        % ROI handle selected in the GUI.
        ROIHandle
    end

    properties (Hidden)
        % List of fields that can be plotted from the GUI.
        DispPlotsOptions = {'DatasetNum', 'FrameNum', 'X', 'Y', 'Z', ...
            'X_SE', 'Y_SE', 'Z_SE', ...
            'Photons', 'Photons_SE', 'Bg', 'Bg_SE', ...
            'PSFSigma', 'PSFSigmaX', 'PSFSigmaY', ...
            'PSFSigma_SE', 'PSFSigmaX_SE', 'PSFSigmaY_SE', ...
            'PValue', 'LogLikelihood', 'ThreshFlag', ...
            'ConnectID', 'NCombined'}

        % Single Molecule Fitting structure loaded along with SMD.
        % NOTE: This is only kept as an aid to exporting, so that both SMD
        %       and SMF can be exported to a file.
        SMF
    end
    
    methods
        function obj = InspectResults(StartGUI)
            if (~exist('StartGUI', 'var') || isempty(StartGUI))
                StartGUI = true;
            end
            if StartGUI
                obj.gui()
            end
        end

        function SMDIsolated = get.SMDIsolated(obj)
            % get method for obj.SMDIsolated.
            if isempty(obj.ROI)
                SMDIsolated = obj.SMD;
                if ~isempty(SMDIsolated)
                    obj.ROI = [1, 1, SMDIsolated.YSize, SMDIsolated.XSize];
                end
            else
                SMDIsolated = smi_core.SingleMoleculeData.isolateSubROI(...
                    obj.SMD, obj.ROI);
            end
        end
        
        [] = gui(obj)
    end


end