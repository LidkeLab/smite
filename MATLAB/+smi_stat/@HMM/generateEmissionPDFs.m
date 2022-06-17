function [PDFHandles] = generateEmissionPDFs(ModelSpecifier)
%generateEmissionPDFs creates array of emission density fxn. handles.
% This method will generate an array of function handles corresponding to 
% the model specified by ModelSpecifier.  These function handles can be
% used to generate the emisssion probability densities for
% each state given observations (as is done in 
% smi_stat.HMM.generateEmissionMatrix()).
%
% NOTE: These function handles are designed to take a single input X, which
%       is a cell array containing all parameters/data needed to
%       compute values for the probability densities.
%       X{1}: Separation in current frame (pixels)
%       X{2}: Standard error in current frame, [particle 1, particle 2]
%             (pixels)(NObservations x 2)
%       X{3}: Frames elapsed since previous observation (DeltaT)
%       X{4}: Sigma overlay, i.e., channel registration error (pixels)
%       X{5}: Dimer separation (pixels)
%       X{6}: Diffusion constant for each trajectory (organized as
%             [D trajectory 1, D trajectory 2]). (pixel^2 / frame)
%       X{7}: Maximum separation between two trajectories (used in
%             computing initial probability density of free state)(pixels)
%       X{8:end}: Additional parameters/data, case dependent, e.g.,
%                 X{8} is the domain separation if ModelSpecifier = 'DDF'
%
% INPUTS:
%   ModelSpecifier: An array specifying the HMM model being used.  The
%                   first character will always be the number of states in
%                   the model.
%                   'DF': two-state dimer/free model
%                   'DDF': three-state dimer/domain/free model
%                   (char array)
%
% OUTPUTS:
%   PDFHandles: A cell array of function handles, with each index
%               corresponding to a state in the model 'ModelSpecifier'.
%               ModelSpecifier 'DF':
%                   FunctionHandles{1} gives the probability density of
%                       observing a separation between two particles in
%                       a dimer.
%                   FunctionHandles{2} gives the probability density of
%                       observing a separation between two particles
%                       that are freely diffusing.
%               ModelSpecifier 'DDF':
%                   FunctionHandles{1} gives the probability density of
%                       observing a separation between two particles in
%                       a dimer.
%                   FunctionHandles{2} gives the probability density of
%                       observing a separation between two particles
%                       that are co-confined (in a "domain").
%                   FunctionHandles{3} gives the probability density of
%                       observing a separation between two particles
%                       that are freely diffusing.

% Created by:
%   David J. Schodt (Lidke Lab, 2020)


% Populate the output array based on the ModelSpecifier.
switch ModelSpecifier
    case 'DF'
        % Generate the function handle for the dimer state.
        PDFHandles{1, 1} = @(X) generateDimerPDF(X);
        
        % Generate the function handle for the free state.
        PDFHandles{2, 1} = @(X) generateFreePDF(X);
    case 'DDF'
        % Generate the function handle for the dimer state.
        PDFHandles{1, 1} = @(X) generateDimerPDF(X);
        
        % Generate the function handle for the domain state.
        PDFHandles{2, 1} = @(X) generateDomainPDF(X);
        
        % Generate the function handle for the free state.
        PDFHandles{3, 1} = @(X) generateFreePDF(X);
end

    function [DimerPDF] = generateDimerPDF(X)
        % This function generates the dimer state probability density.
        
        % Extract some arrays from X to make the code more readable.
        Separation = X{1};
        PositionSE = X{2};
        SigmaOverlay = X{4};
        SeparationDimer = X{5};
        
        % Define the dimer state pdf.
        VarianceDimer = sum(PositionSE.^2, 2) + SigmaOverlay^2;
        DimerPDF = (Separation./VarianceDimer) ...
            .* exp(-0.5*(Separation.^2+SeparationDimer^2)./VarianceDimer) ...
            .* besseli(0, Separation*SeparationDimer./VarianceDimer);
        IsInvalid = (isinf(DimerPDF) | isnan(DimerPDF));
        if any(IsInvalid)
            % The Bessel function is unstable for large arguments, so we
            % sometimes need to compute the integral numerically at certain
            % points.
            InvalidIndices = find(IsInvalid);
            Theta = linspace(0, 2*pi, 1e3);
            DimerPDFInt = zeros(sum(IsInvalid), 1);
            for ii = 1:sum(IsInvalid)
                Index = InvalidIndices(ii);
                ProbIntegrand = ...
                    exp(-0.5 * (Separation(Index)^2+SeparationDimer^2 ...
                    - 2*Separation(Index)*SeparationDimer...
                    *cos(Theta))/VarianceDimer(Index));
                DimerPDFInt(ii) = ...
                    (Separation(Index)/(2*pi*VarianceDimer(Index))) ...
                    * trapz(Theta, ProbIntegrand);
            end
            DimerPDF(IsInvalid) = DimerPDFInt;
        end
    end

    function [DomainPDF] = generateDomainPDF(X)
        % This function generates the domain state probability density.
        
        % Extract some arrays from X to make the code more readable.
        Separation = X{1};
        PositionSE = X{2};
        SigmaOverlay = X{4};
        SeparationDomain = X{8};
        
        % Define the domain state pdf.
        VarianceDomain = sum(PositionSE.^2, 2) + SigmaOverlay^2;
        DomainPDF = (Separation./VarianceDomain) ...
            .* exp(-0.5*(Separation.^2+SeparationDomain^2)./VarianceDomain) ...
            .* besseli(0, Separation*SeparationDomain./VarianceDomain);
        IsInvalid = (isinf(DomainPDF) | isnan(DomainPDF));
        if any(IsInvalid)
            % The Bessel function is unstable for large arguments, so we
            % sometimes need to compute the integral numerically for some
            % observations.
            InvalidIndices = find(IsInvalid);
            Theta = linspace(0, 2*pi, 1e3);
            DomainPDFInt = zeros(sum(IsInvalid), 1);
            for ii = 1:sum(IsInvalid)
                Index = InvalidIndices(ii);
                ProbIntegrand = ...
                    exp(-0.5 * (Separation(Index)^2+SeparationDomain^2 ...
                    - 2*Separation(Index)*SeparationDomain...
                    *cos(Theta))/VarianceDomain(Index));
                DomainPDFInt(ii) = ...
                    (Separation(Index)/(2*pi*VarianceDomain(Index))) ...
                    * trapz(Theta, ProbIntegrand);
            end
            DomainPDF(IsInvalid) = DomainPDFInt;
        end
    end

    function [FreePDF] = generateFreePDF(X)
        % This function generates the free state probability density.
        
        % Define useful indexing arrays.
        NObservations = numel(X{1});
        if (NObservations == 1)
            IndexArrayCurrent = 1;
            IndexArrayPrevious = 1;
        else
            IndexArrayCurrent = 2:NObservations;
            IndexArrayPrevious = 1:(NObservations-1);
        end
        
        % Extract some arrays from X to make the code more readable.
        InitialSeparation = X{1}(1);
        Separation = X{1}(IndexArrayCurrent);
        SeparationPrevious = X{1}(IndexArrayPrevious);
        PositionSE = X{2}(IndexArrayCurrent, :);
        PositionSEPrevious = X{2}(IndexArrayPrevious, :);
        DeltaT = X{3};
        DiffusionCoefficients = X{6}(IndexArrayPrevious, :);
        
        % Define the free state pdf.
        VarianceFree = sum(PositionSE.^2, 2) ...
            + sum(PositionSEPrevious.^2, 2) ...
            + 2*sum(DiffusionCoefficients, 2).*DeltaT;
        VarianceFreeInitial = 2*sum(PositionSEPrevious(1, :).^2, 2) ...
            + 2*sum(DiffusionCoefficients(1, :), 2);
        FreePDF = (Separation./VarianceFree) ...
            .* exp(-0.5*(Separation.^2+SeparationPrevious.^2)./VarianceFree) ...
            .* besseli(0, Separation.*SeparationPrevious./VarianceFree);
        FreePDFInitial = (InitialSeparation./VarianceFreeInitial) ...
            .* exp(-InitialSeparation^2./VarianceFreeInitial) ...
            * besseli(0, InitialSeparation^2./VarianceFreeInitial);
        IsInvalid = (isnan(FreePDF) | isinf(FreePDF));
        IsInvalidInitial = (isnan(FreePDFInitial) ...
            || isinf(FreePDFInitial));
        if any(IsInvalid)
            % The Bessel function is unstable for large arguments, so we
            % sometimes need to compute the integral numerically at certain
            % points.
            InvalidIndices = find(IsInvalid);
            Theta = linspace(0, 2*pi, 1e3);
            FreePDFInt = zeros(sum(IsInvalid), 1);
            for ii = 1:sum(IsInvalid)
                Index = InvalidIndices(ii);
                ProbIntegrand = exp(-0.5...
                    *(Separation(Index)^2+SeparationPrevious(Index)^2 ...
                    - 2*Separation(Index)*SeparationPrevious(Index)...
                    *cos(Theta))/VarianceFree(Index));
                FreePDFInt(ii) = ...
                    (Separation(Index)/(2*pi*VarianceFree(Index))) ...
                    * trapz(Theta, ProbIntegrand);
            end
            FreePDF(IsInvalid) = FreePDFInt;
        end
        if any(IsInvalidInitial)
            % The Bessel function is unstable for large arguments, so we
            % sometimes need to compute the integral numerically at certain
            % points.
            Theta = linspace(0, 2*pi, 1e3);
            ProbIntegrand = ...
                exp(-(InitialSeparation^2) * (1-cos(Theta)) ...
                ./ VarianceFreeInitial);
            FreePDFInitial = ...
                (InitialSeparation ./ (2*pi*VarianceFreeInitial)) ...
                * trapz(Theta, ProbIntegrand);
        end
        FreePDF = [FreePDFInitial; FreePDF];
    end


end