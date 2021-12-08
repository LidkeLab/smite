classdef InspectResults < handle
    %InspectResults contains methods useful for inspecting SR data.
    % This class provides an interface for inspecting super-resolution
    % results.  The primary intention is that this class can be used to
    % associate localizations from a reconstruction image (e.g., a Gaussian
    % SR image) and the entries of a Single Molecule Data structure
    % (see smi_core.SingleMoleculeData).
    
    properties
        % Single Molecule Data structure. (see smi_core.SingleMoleculeData)
        SMD

        % Super-resolution reconstruction image. (YxXx3 float)
        SRImage

        % Figure containing the interactive GUI.
        GUIFigure

        % Axes containing the SR image.
        ImageAxes
    end

    properties (SetAccess = 'protected')
        % ROI handle selected in the GUI.
        ROIHandle

        % ROI in SMD coordinates. ([YStart, XStart, YEnd, XEnd])
        ROI

        % Isolated SMD within the ROI defined by ROIHandle.
        SMDIsolated
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
    end
    
    methods
        function obj = InspectResults()
        end
        
        [] = gui(obj)
    end


end