classdef ChangePoint < handle
    properties (SetAccess=protected)
        %Input variables
        data; % Input data given to constructor or set with setData
        logBayesThreshold; % Input threshold for deciding on signficance level of bayes factors
        Nobservations; % Length of data

        %Output variables
        NchangePoints; % Output - The number of change points detected
        changePoints; % Output - The discrete times at which a change point was detected (size=NchangePoints)
        intensity; % Output - The mean intensities for each subsequence  (size=NchangePoints+1)
        logBayesFactors; % Output - The bayes factors for each accepted change points (size=NchangePoints)
        NrejectedChangePoints; % Output - The number of rejected change points
        rejectedChangePoints;  % Output - The rejected change points (size=NrejectedChangePoints)
        rejectedLogBayesFactors;  % Output - The bayes factors for each rejected change points (size=NrejectedChangePoints)
    end

    methods
        function obj=IntensityCPA(data, logBayesThreshold)
            % Inputs:
            %   data: a vector of integer values representing poission
            %       distributed values from a sequence where the mean of the poisson
            %       process changes at discrete times
            %   logBayesThreshold: the log of the Bayes factor to accept a
            %       change point with.  Should be positive. Larger values will estimate 
            %       fewer change points.  Values will need to be larger when the counts 
            %       in the data are bigger.   Start with 1e1 and increase as needed for lower
            %       sensitivity.
            obj.setData(data, logBayesThreshold);
        end

        function setData(obj, data, logBayesThreshold)
            % Inputs:
            %   data: a vector of integer values representing poission
            %       distributed values from a sequence where the mean of the poisson
            %       process changes at discrete times
            %   logBayesThreshold: the log of the Bayes factor to accept a
            %       change point with.  Should be negative. Smaller values will estimate fewer change
            %       points.  Values will need to be smaller when the counts in the
            %       data are bigger.   Try something in the range of -10 to -1E6.
            if ~isvector(data)
                error('IntensityCPA:setData','data must be a vector')
            end
            if ~all(round(data)==data)
                error('IntensityCPA:setData','data must be integer valued')
            end
            obj.data=data;
            obj.Nobservations=length(obj.data);
            if nargin==3
                if logBayesThreshold<=0
                    error('IntensityCPA:setData','logBayesThreshold should be positive')
                end
                obj.logBayesThreshold=logBayesThreshold;
            end
            obj.estimateSequence();
        end
        
        function f=plotIntensityEstimate(obj)
            f=figure();
            xs=1:obj.Nobservations;
            subplot(2,1,1);
            stairs(xs, obj.data, '-ko');
            hold on;
            Is=IntensityCPA.modelIntensity(obj.Nobservations, obj.changePoints, obj.intensity);
            stairs(xs, Is, '-r', 'LineWidth', 2.0);
            yl=ylim;
            ylim([0,yl(2)]);
            hold off;
            xlabel('time');
            ylabel('intensity');
            legend('Data', 'Estimated', 'Location', 'North');
            subplot(2,1,2);
            xs=2:obj.Nobservations;
            pCP=arrayfun(@(cp) IntensityCPA.logP_H2(obj.data,cp), xs);
            plot(xs, pCP, '-ko');
            xlabel('time');
            ylabel('logP(cp|H_2)');
        end
    end %public methods

    methods(Static=true)
        function data=simulate(nObservations, changePoints, intensity)
            % Simulate an intensity sequence
            % Inputs:
            %   nObservations - Scalar integer: length of data sequence
            %   changePoints - vector: indexs of change point locations
            %   intensity - vector: length=length(changePoints)+1: mean intensities
            %       for each subinterval
            % Outputs:
            %   data - A 1xnObservations vector of integer intensity values for
            %      specified change point sequence
            assert(length(intensity)==length(changePoints)+1);%Correct number of intenisites
            assert(all(changePoints>1)); % Not too small
            assert(all(changePoints)<=nObservations); %Not too big
            assert(length(unique(changePoints))==length(changePoints)); %No identical change points
            assert(all(sort(changePoints)==changePoints)); % alrady sorted
            assert(all(intensity(2:end)~=intensity(1:end-1))); %Change points actually change intensity
            data=zeros(1,nObservations);
            nCP=length(changePoints);
            cp=1;
            for n=1:nObservations
                if cp<=nCP && changePoints(cp)==n
                    cp=cp+1;
                end
                data(n)=poissrnd(intensity(cp),1,1);
            end
        end

        function [icp,f]=plotSimulatedEstimate(nObservations, changePoints, intensity, logBayesThreshold)
            % Given a fixed set of changepoints and intensities, 
            % simulate, analyze and plot a data sets.  Arguments are similar to
            % the simulate method.
            % Inputs:
            %   nObservations - Scalar integer: length of data sequence
            %   changePoints - vector: indexs of change point locations
            %   intensity - vector: length=length(changePoints)+1: mean intensities
            %       for each subinterval
            %   logBayesThreshold - positive scalar.  Threshold for accepting
            %       change points
            % Outputs:
            %   icp - The IntensityCPA for this data
            %   f - figure handle
            data=IntensityCPA.simulate(nObservations, changePoints, intensity);
            icp=IntensityCPA(data,logBayesThreshold);
            f=figure();
            xs=1:nObservations;
            subplot(2,1,1);
            stairs(xs, data, '-ko');
            hold on;
            modelIs=IntensityCPA.modelIntensity(nObservations, changePoints, intensity);
            estimatedIs=IntensityCPA.modelIntensity(nObservations, icp.changePoints, icp.intensity);
            stairs(xs, modelIs , '-b','LineWidth', 2.0);
            stairs(xs, estimatedIs, '-r', 'LineWidth', 2.0);
            yl=ylim;
            ylim([0,yl(2)]);
            hold off;
            xlabel('time');
            ylabel('intensity');
            legend('Data', 'Model', 'Estimated', 'Location', 'North');
            subplot(2,1,2);
            xs=2:nObservations;
            pCP=arrayfun(@(cp) IntensityCPA.logP_H2(data,cp), xs);
            plot(xs, pCP, '-ko');
            xlabel('time');
            ylabel('logP(cp|H_2)');
        end

        function [icp,f]=plotRandSimulatedEstimate(nObservations, nChangePoints, meanIntensity, logBayesThreshold)
            % Generate a random sequence of change points and intensities, then 
            % simulate, analyze and plot the data.
            % Inputs:
            %   nObservations - Scalar integer: length of data sequence
            %   nChangePoints - Scalar integer: number of change points to
            %                   simulate
            %   meanIntensity - scalar for mean intensity.  Intensities are
            %                   uniformly distributed on [1, 2*meanIntensity]
            %   logBayesThreshold - positive scalar.  Threshold for accepting
            %                       change points
            % Outputs:
            %   icp - The IntensityCPA for this data
            %   f - figure handle
            changePoints=sort(round(rand(1,nChangePoints)*(nObservations-1)+2));
            while length(unique(changePoints))~=length(changePoints)
                changePoints=sort(round(rand(1,nChangePoints)*(nObservations-1)+2));
            end
            intensity=round(1+rand(1,nChangePoints+1)*2*meanIntensity);
            while any(intensity(2:end)==intensity(1:end-1))
                intensity=round(1+rand(1,nChangePoints+1)*2*meanIntensity);
            end
            [icp,f]=IntensityCPA.plotSimulatedEstimate(nObservations, changePoints, intensity, logBayesThreshold);
        end
    end %public static methods


    methods(Access=protected)
        function estimateSequence(obj)
            %
            % This is the primary estimation routine entry point
            %
            cP=obj.estimateChangePoints(obj.data, obj.logBayesThreshold);
            [obj.changePoints, obj.logBayesFactors, obj.rejectedChangePoints, ...
                obj.rejectedLogBayesFactors]=obj.filterChangePoints(obj.data,cP, obj.logBayesThreshold);
            obj.NchangePoints=numel(obj.changePoints);
            obj.NrejectedChangePoints=numel(obj.rejectedChangePoints);
            if obj.NchangePoints==0
                obj.intensity=obj.estimateIntensity(obj.data);
                return
            end
            obj.intensity=zeros(1, obj.NchangePoints+1);
            for n=1:obj.NchangePoints+1
                if n==1
                    first=1;
                else
                    first=obj.changePoints(n-1);
                end
                if n>obj.NchangePoints
                    last=length(obj.data);
                else
                    last=obj.changePoints(n)-1;
                end
                obj.intensity(n)=IntensityCPA.estimateIntensity(obj.data(first:last));
            end
        end
    end % Protected methods

    methods(Static=true)
        
        function [changePoints, accLogBayesFactors, rejLogBayesFactors]=estimateChangePoints(data,logBayesFactorThreshold)
            % This is the main recursive method for identifiying change points
            % since it and other helper methods work on parts of the data set
            % recursivly we make it a static method that takes in only the
            % portion of the data that it is to operate on
            logBF=IntensityCPA.logBayesFactor(data);
            if ~isempty(logBF) && logBF>=logBayesFactorThreshold && numel(data)>=2
                cp=IntensityCPA.estimateChangePointLocation(data);
                data1=data(1:cp-1);
                data2=data(cp:end);
                [cP1, accLogBF1, rejLogBF1]=IntensityCPA.estimateChangePoints(data1,logBayesFactorThreshold);
                [cP2, accLogBF2, rejLogBF2]=IntensityCPA.estimateChangePoints(data2,logBayesFactorThreshold);
                accLogBayesFactors=[accLogBF1 logBF accLogBF2];
                rejLogBayesFactors=[rejLogBF1 rejLogBF2];
                changePoints=[cP1 cp cP2+cp-1];
            else
                changePoints=[];
                rejLogBayesFactors=logBF;
                accLogBayesFactors=[];
            end
        end

        function [accepted, accLogBayesFactors, rejected, rejLogBayesFactors]=filterChangePoints(data,changePoints, logBayesFactorThreshold)
            nChangePoints=length(changePoints);
            
            if nChangePoints==0
                accepted=[];
                accLogBayesFactors=[];
                rejected=[];
                rejLogBayesFactors=[];
                return
            end
            lBF=zeros(1,nChangePoints);
            for n=1:nChangePoints
                if n==1
                    first=1;
                else
                    first=changePoints(n-1);
                end
                if n==nChangePoints
                    last=length(data);
                else
                    last=changePoints(n+1)-1;
                end
                lBF(n)=IntensityCPA.logBayesFactorPoint(data(first:last),changePoints(n)-(first-1));
            end
            rejIdx=find(lBF<logBayesFactorThreshold);
            accIdx=setdiff(1:nChangePoints,rejIdx);
            accepted=changePoints(accIdx);
            accLogBayesFactors=lBF(accIdx);
            rejected=changePoints(rejIdx);
            rejLogBayesFactors=lBF(rejIdx);
            if ~isempty(rejected)
                [accepted, accLogBayesFactors, new_rejected, new_rejLogBayesFactors]=IntensityCPA.filterChangePoints(data,accepted, logBayesFactorThreshold);
                rejected=[rejected new_rejected];
                rejLogBayesFactors=[rejLogBayesFactors new_rejLogBayesFactors];
                [rejected, shuffle]=sort(rejected);
                rejLogBayesFactors=rejLogBayesFactors(shuffle);
            end
        end

        function logP=logP_H1(data)
            % Compute the unnormalized log probability of hypothesis H1 (there
            % are no change points and all data comes from the same poisson
            % distribution)
            % This is unnormalized as we don't have the Jeffery's prior constant
            N=length(data);
            C=sum(data);
            logR=-logFactorialSum(data);
            logP=logR+(C-.5)*log(C)-C*(1+log(N))+.5*log(2*pi);
        end

        function logP=logP_H2(data, change_point)
            % Compute the unnormalized log probability of hypothesis H2 (there
            % is a change point at change_points from a poisson process with one mean
            % to a poisson process witha second mean after marginalizing accross all possible
            % intensities).
            % This is unnormalized as we don't have the Jeffery's prior constant
            N=length(data);
            assert(change_point>1 && change_point<=N);
            N1=change_point-1;
            N2=N-N1;
            C1=sum(data(1:N1));
            C2=sum(data(N1+1:end));
            
            term1=(2/pi)/((C1/N1)^2+(C2/N2)^2);
            
            logR=-logFactorialSum(data);
            LFC1=logFactorial(C1);
            LFC2=logFactorial(C2);

            logP=log(term1)+logR+LFC1+LFC2-(C1+1)*log(N1)-(C2+1)*log(N2);
        end

        function logP=mean_LogP_H2(data)
            logP=log(mean(arrayfun(@(cp) exp(IntensityCPA.logP_H2(data, cp)), 2:length(data))));
        end

        function logBF=logBayesFactor(data)
            N=length(data);
            C=sum(data);
            H2LogTerms=zeros(1,N-1);
            C1=data(1);
            C2=C-C1;
            N1=1;
            N2=N-N1;
            for cp=2:N
                H2LogTerms(cp-1)= logFactorial(C1)+logFactorial(C2)-(C1+1)*log(N1)-(C2+1)*log(N2) ...
                                    - log((C1/N1)^2+(C2/N2)^2);
                C1=C1+data(cp);
                C2=C2-data(cp);
                N1=N1+1;
                N2=N2-1;
            end
            logBF=log(2/pi) + C*log(N) - logFactorial(C-1)-log(N-1)+logSum(H2LogTerms);
        end

        function logBF=logBayesFactorPoint(data, changePoint)
            logBF=IntensityCPA.logP_H2(data,changePoint) - IntensityCPA.logP_H1(data);
        end

        function logPs=logP_ChangePoint(data)
            logPs=zeros(1,length(data));
            N=length(data);
            C=sum(data);
            C1=data(1);
            C2=C-C1;
            N1=1;
            N2=N-N1;
            for cp=2:length(data)
                logPs(cp)=logFactorial(C1)+logFactorial(C2)-(C1+1)*log(N1)-(C2+1)*log(N2)...
                                - log((C1/N1)^2+(C2/N2)^2);
                C1=C1+data(cp);
                C2=C2-data(cp);
                N1=N1+1;
                N2=N2-1;
            end
        end

        function cp=estimateChangePointLocation(data)
            logPs=IntensityCPA.logP_ChangePoint(data);
            [~,cp]=max(logPs(2:end)); %the first time cannot be a change point
            cp=cp+1; %correct for dropping first point
        end

        function logP=logPLambda(data, lambda)
            N=length(data);
            C=sum(data);
            logP=C*log(N) + (C-1)*log(lambda) - N*lambda - logFactorial(C-1);
        end

        function lambda=estimateIntensity(data)
            N=length(data);
            C=sum(data);
            lambda=C/N;
        end

        function model_intensity=modelIntensity(nObservations, changePoints, intensity)
            % Helper for plotting an input intensity sequence
            %
            model_intensity=zeros(1,nObservations);
            nCP=length(changePoints);
            cp=1;
            for n=1:nObservations
                if cp<=nCP && changePoints(cp)==n
                    cp=cp+1;
                end
                model_intensity(n)=intensity(cp);
            end
        end

    end %static protected methods
end
