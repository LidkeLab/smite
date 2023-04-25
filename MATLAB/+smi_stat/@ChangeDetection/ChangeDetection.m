classdef ChangeDetection < handle
    %ChangeDetection contains methods for change detection analysis.
    % This class contains several methods used to detect change points 
    % from a sequence of intensity data, e.g. intensity trace from 
    % single-particle tracking.
    %
    % See reference 
    % Daniel L. Ensign and Vijay S. Pande, Bayesian (2010), 
    % Detection of Intensity Changes in Single Molecule and 
    % Molecular Dynamics Trajectories, J. Phys. Chem. B 2010
    % https://doi.org/10.1021/jp906786b
    % 
    % Created by:
    %   Mark J. Olah (Lidke Lab, 2014)

    properties (SetAccess=protected)
        %Input variables
        Data; % Input data given to constructor or set with setData
        LogBayesThreshold; % Input threshold for deciding on signficance level of bayes factors
        Nobservations; % Length of data

        %Output variables
        NchangePoints; % Output - The number of change points detected
        ChangePoints; % Output - The discrete times at which a change point was detected (size=NchangePoints)
        Intensity; % Output - The mean intensities for each subsequence  (size=NchangePoints+1)
        IntensityModel; % Output - Intensity of model sequence
        LogBayesFactors; % Output - The bayes factors for each accepted change points (size=NchangePoints)
        NrejectedChangePoints; % Output - The number of rejected change points
        RejectedChangePoints;  % Output - The rejected change points (size=NrejectedChangePoints)
        RejectedLogBayesFactors;  % Output - The bayes factors for each rejected change points (size=NrejectedChangePoints)
    end

    methods
        function obj=ChangeDetection(Data, LogBayesThreshold)
            % Inputs:
            %   data: a vector of integer values representing poission
            %       distributed values from a sequence where the mean of the poisson
            %       process changes at discrete times
            %   LogBayesThreshold: the log of the Bayes factor to accept a
            %       change point with.  Should be positive. Larger values will estimate 
            %       fewer change points.  Values will need to be larger when the counts 
            %       in the data are bigger.   Start with 1e1 and increase as needed for lower
            %       sensitivity.
            obj.setData(Data, LogBayesThreshold);
        end

        function setData(obj, Data, LogBayesThreshold)
            % Inputs:
            %   Data: a vector of integer values representing poission
            %       distributed values from a sequence where the mean of the poisson
            %       process changes at discrete times
            %   LogBayesThreshold: the log of the Bayes factor to accept a
            %       change point with.  Should be negative. Smaller values will estimate fewer change
            %       points.  Values will need to be smaller when the counts in the
            %       data are bigger.   Try something in the range of -10 to -1E6.
            if ~isvector(Data)
                error('ChangePoint:setData','data must be a vector')
            end
            if ~all(round(Data)==Data)
                error('ChangePoint:setData','data must be integer valued')
            end
            obj.Data=Data;
            obj.Nobservations=length(obj.Data);
            if nargin==3
                if LogBayesThreshold<=0
                    error('ChangePoint:setData','logBayesThreshold should be positive')
                end
                obj.LogBayesThreshold=LogBayesThreshold;
            end
            obj.estimateSequence();
        end
        
        function F=plotIntensityEstimate(obj)
            % plot estimated step trace and probablity of change point at
            % each observation
            F=figure();
            xs=1:obj.Nobservations;
            subplot(2,1,1);
            stairs(xs, obj.Data, '-ko');
            hold on;
            Is=obj.IntensityModel;
            stairs(xs, Is, '-r', 'LineWidth', 2.0);
            yl=ylim;
            ylim([0,yl(2)]);
            hold off;
            xlabel('time');
            ylabel('intensity');
            legend('Data', 'Estimated', 'Location', 'North');
            subplot(2,1,2);
            xs=2:obj.Nobservations;
            pCP=arrayfun(@(CP) smi_stat.ChangeDetection.logP_H2(obj.Data,CP), xs);
            plot(xs, pCP, '-ko');
            xlabel('time');
            ylabel('logP(cp|H_2)');
        end
    end %public methods

    methods(Static=true)
        success = unitTest()

        function Data=simulate(NObservations, ChangePoints, Intensity)
            % Simulate an intensity sequence
            % Inputs:
            %   NObservations - Scalar integer: length of data sequence
            %   ChangePoints - vector: indexs of change point locations
            %   intensity - vector: length=length(ChangePoints)+1: mean intensities
            %       for each subinterval
            % Outputs:
            %   Data - A 1xNObservations vector of integer intensity values for
            %      specified change point sequence
            assert(length(Intensity)==length(ChangePoints)+1);%Correct number of intenisites
            assert(all(ChangePoints>1)); % Not too small
            assert(all(ChangePoints)<=NObservations); %Not too big
            assert(length(unique(ChangePoints))==length(ChangePoints)); %No identical change points
            assert(all(sort(ChangePoints)==ChangePoints)); % alrady sorted
            assert(all(Intensity(2:end)~=Intensity(1:end-1))); %Change points actually change intensity
            Data=zeros(1,NObservations);
            nCP=length(ChangePoints);
            CP=1;
            for n=1:NObservations
                if CP<=nCP && ChangePoints(CP)==n
                    CP=CP+1;
                end
                Data(n)=poissrnd(Intensity(CP),1,1);
            end
        end

        function [Icp,F]=plotSimulatedEstimate(NObservations, ChangePoints, Intensity, LogBayesThreshold)
            % Given a fixed set of changepoints and intensities, 
            % simulate, analyze and plot a data sets.  Arguments are similar to
            % the simulate method.
            % Inputs:
            %   NObservations - Scalar integer: length of data sequence
            %   ChangePoints - vector: indexs of change point locations
            %   Intensity - vector: length=length(ChangePoints)+1: mean intensities
            %       for each subinterval
            %   LogBayesThreshold - positive scalar.  Threshold for accepting
            %       change points
            % Outputs:
            %   Icp - The object instance of ChangeDetection
            %   f - figure handle
            Data=smi_stat.ChangeDetection.simulate(NObservations, ChangePoints, Intensity);
            Icp=smi_stat.ChangeDetection(Data,LogBayesThreshold);
            F=figure();
            xs=1:NObservations;
            subplot(2,1,1);
            stairs(xs, Data, '-ko');
            hold on;
            modelIs=smi_stat.ChangeDetection.modelIntensity(NObservations, ChangePoints, Intensity);
            estimatedIs=Icp.IntensityModel;
            stairs(xs, modelIs , '-b','LineWidth', 2.0);
            stairs(xs, estimatedIs, '-r', 'LineWidth', 2.0);
            yl=ylim;
            ylim([0,yl(2)]);
            hold off;
            xlabel('time');
            ylabel('Intensity');
            legend('Data', 'Model', 'Estimated', 'Location', 'North');
            subplot(2,1,2);
            xs=2:NObservations;
            pCP=arrayfun(@(CP) smi_stat.ChangeDetection.logP_H2(Data,CP), xs);
            plot(xs, pCP, '-ko');
            xlabel('time');
            ylabel('logP(cp|H_2)');
        end

        function [Icp,F]=plotRandSimulatedEstimate(NObservations, NChangePoints, meanIntensity, LogBayesThreshold)
            % Generate a random sequence of change points and intensities, then 
            % simulate, analyze and plot the data.
            % Inputs:
            %   NObservations - Scalar integer: length of data sequence
            %   NChangePoints - Scalar integer: number of change points to
            %                   simulate
            %   meanIntensity - scalar for mean intensity.  Intensities are
            %                   uniformly distributed on [1, 2*meanIntensity]
            %   LogBayesThreshold - positive scalar.  Threshold for accepting
            %                       change points
            % Outputs:
            %   Icp - The ChangePoint for this data
            %   F - figure handle
            ChangePoints=sort(round(rand(1,NChangePoints)*(NObservations-1)+2));
            while length(unique(ChangePoints))~=length(ChangePoints)
                ChangePoints=sort(round(rand(1,NChangePoints)*(NObservations-1)+2));
            end
            Intensity=round(1+rand(1,NChangePoints+1)*2*meanIntensity);
            while any(Intensity(2:end)==Intensity(1:end-1))
                Intensity=round(1+rand(1,NChangePoints+1)*2*meanIntensity);
            end
            [Icp,F]=smi_stat.ChangeDetection.plotSimulatedEstimate(NObservations, ChangePoints, Intensity, LogBayesThreshold);
        end
    end %public static methods


    methods(Access=protected)
        function estimateSequence(obj)
            %
            % This is the primary estimation routine entry point
            %
            cP=obj.estimateChangePoints(obj.Data, obj.LogBayesThreshold);
            [obj.ChangePoints, obj.LogBayesFactors, obj.RejectedChangePoints, ...
                obj.RejectedLogBayesFactors]=obj.filterChangePoints(obj.Data,cP, obj.LogBayesThreshold);
            obj.NchangePoints=numel(obj.ChangePoints);
            obj.NrejectedChangePoints=numel(obj.RejectedChangePoints);
            if obj.NchangePoints==0
                obj.Intensity=obj.estimateIntensity(obj.Data);
                return
            end
            obj.Intensity=zeros(1, obj.NchangePoints+1);
            for n=1:obj.NchangePoints+1
                if n==1
                    first=1;
                else
                    first=obj.ChangePoints(n-1);
                end
                if n>obj.NchangePoints
                    last=length(obj.Data);
                else
                    last=obj.ChangePoints(n)-1;
                end
                obj.Intensity(n)=smi_stat.ChangeDetection.estimateIntensity(obj.Data(first:last));
            end
            obj.IntensityModel=smi_stat.ChangeDetection.modelIntensity(obj.Nobservations, obj.ChangePoints, obj.Intensity);
        end
    end % Protected methods

    methods(Static=true)
        
        function [ChangePoints, AccLogBayesFactors, RejLogBayesFactors]=estimateChangePoints(Data,LogBayesFactorThreshold)
            % This is the main recursive method for identifiying change points
            % since it and other helper methods work on parts of the data set
            % recursivly we make it a static method that takes in only the
            % portion of the data that it is to operate on
            % Inputs:
            %   data: a vector of integer values representing poission
            %       distributed values from a sequence where the mean of the poisson
            %       process changes at discrete times
            %   LogBayesThreshold - positive scalar.  Threshold for accepting
            %                       change points
            % Outputs:
            %   ChangePoints - The ChangePoint for this data
            %   AccLogBayesFactors - accepted logBayesThreshold
            %   RejLogBayesFactors - rejected logBayesThreshold

            LogBF=smi_stat.ChangeDetection.logBayesFactor(Data);
            if ~isempty(LogBF) && LogBF>=LogBayesFactorThreshold && numel(Data)>=2
                %split data into two parts at the change point
                CP=smi_stat.ChangeDetection.estimateChangePointLocation(Data);
                Data1=Data(1:CP-1);
                Data2=Data(CP:end);
                [cP1, accLogBF1, rejLogBF1]=smi_stat.ChangeDetection.estimateChangePoints(Data1,LogBayesFactorThreshold);
                [cP2, accLogBF2, rejLogBF2]=smi_stat.ChangeDetection.estimateChangePoints(Data2,LogBayesFactorThreshold);
                AccLogBayesFactors=[accLogBF1 LogBF accLogBF2];
                RejLogBayesFactors=[rejLogBF1 rejLogBF2];
                ChangePoints=[cP1 CP cP2+CP-1];
            else
                ChangePoints=[];
                RejLogBayesFactors=LogBF;
                AccLogBayesFactors=[];
            end
        end

        function [Accepted, AccLogBayesFactors, Rejected, RejLogBayesFactors]=filterChangePoints(Data,ChangePoints, LogBayesFactorThreshold)
            %remove change points that are below the logBayesFactorThreshold
            % Inputs:
            %   Data: a vector of integer values representing poission
            %       distributed values from a sequence where the mean of the poisson
            %       process changes at discrete times
            %   ChangePoints - The ChangePoint for this data
            %
            %   LogBayesFactorThreshold - positive scalar.  Threshold for accepting
            %                             change points
            % Outputs:
            %   Accepted - accepted Change Points
            %   AccLogBayesFactors - accepted logBayesThreshold
            %   Rejected - rejected Change Points
            %   RejLogBayesFactors - rejected logBayesThreshold

            NChangePoints=length(ChangePoints);
            
            if NChangePoints==0
                Accepted=[];
                AccLogBayesFactors=[];
                Rejected=[];
                RejLogBayesFactors=[];
                return
            end
            lBF=zeros(1,NChangePoints);
            for n=1:NChangePoints
                if n==1
                    first=1;
                else
                    first=ChangePoints(n-1);
                end
                if n==NChangePoints
                    last=length(Data);
                else
                    last=ChangePoints(n+1)-1;
                end
                lBF(n)=smi_stat.ChangeDetection.logBayesFactorPoint(Data(first:last),ChangePoints(n)-(first-1));
            end
            rejIdx=find(lBF<LogBayesFactorThreshold);
            accIdx=setdiff(1:NChangePoints,rejIdx);
            Accepted=ChangePoints(accIdx);
            AccLogBayesFactors=lBF(accIdx);
            Rejected=ChangePoints(rejIdx);
            RejLogBayesFactors=lBF(rejIdx);
            if ~isempty(Rejected)
                [Accepted, AccLogBayesFactors, New_rejected, New_rejLogBayesFactors]=smi_stat.ChangeDetection.filterChangePoints(Data,Accepted, LogBayesFactorThreshold);
                Rejected=[Rejected New_rejected];
                RejLogBayesFactors=[RejLogBayesFactors New_rejLogBayesFactors];
                [Rejected, shuffle]=sort(Rejected);
                RejLogBayesFactors=RejLogBayesFactors(shuffle);
            end
        end

        function LogP=logP_H1(Data)
            % Compute the unnormalized log probability of hypothesis H1 (there
            % are no change points and all data comes from the same poisson
            % distribution)
            % This is unnormalized as we don't have the Jeffery's prior constant
            N=length(Data);
            C=sum(Data);
            logR=-logFactorialSum(Data);
            LogP=logR+(C-.5)*log(C)-C*(1+log(N))+.5*log(2*pi);
        end

        function LogP=logP_H2(Data, change_point)
            % Compute the unnormalized log probability of hypothesis H2 (there
            % is a change point at change_points from a poisson process with one mean
            % to a poisson process witha second mean after marginalizing accross all possible
            % intensities).
            % This is unnormalized as we don't have the Jeffery's prior constant
            N=length(Data);
            assert(change_point>1 && change_point<=N);
            N1=change_point-1;
            N2=N-N1;
            C1=sum(Data(1:N1));
            C2=sum(Data(N1+1:end));
            
            term1=(2/pi)/((C1/N1)^2+(C2/N2)^2);
            
            logR=-logFactorialSum(Data);
            LFC1=logFactorial(C1);
            LFC2=logFactorial(C2);

            LogP=log(term1)+logR+LFC1+LFC2-(C1+1)*log(N1)-(C2+1)*log(N2);
        end

        function LogP=mean_LogP_H2(Data)
            % compute the total log probability of hypothesis H2 (there
            % is a change point from the data
            LogP=log(mean(arrayfun(@(CP) exp(smi_stat.ChangeDetection.logP_H2(Data, CP)), 2:length(Data))));
        end

        function LogBF=logBayesFactor(Data)
            %compute the total Bayes Factor, which is total log probability
            %of hypothesis H2 (one change point) minus the log probability of
            %hypothesis H1 (no change point)
            N=length(Data);
            C=sum(Data);
            H2LogTerms=zeros(1,N-1);
            C1=Data(1);
            C2=C-C1;
            N1=1;
            N2=N-N1;
            for CP=2:N
                H2LogTerms(CP-1)= logFactorial(C1)+logFactorial(C2)-(C1+1)*log(N1)-(C2+1)*log(N2) ...
                                    - log((C1/N1)^2+(C2/N2)^2);
                C1=C1+Data(CP);
                C2=C2-Data(CP);
                N1=N1+1;
                N2=N2-1;
            end
            LogBF=log(2/pi) + C*log(N) - logFactorial(C-1)-log(N-1)+logSum(H2LogTerms);
        end

        function LogBF=logBayesFactorPoint(Data, ChangePoint)
            %compute the Bayes Factor at the change point, which is the
            %log probability of hypothesis H2 (one change point) minus the log probability of
            %hypothesis H1 (no change point)
            LogBF=smi_stat.ChangeDetection.logP_H2(Data,ChangePoint) - smi_stat.ChangeDetection.logP_H1(Data);
        end

        function LogPs=logP_ChangePoint(Data)
            %compute the log probability of the change point
            LogPs=zeros(1,length(Data));
            N=length(Data);
            C=sum(Data);
            C1=Data(1);
            C2=C-C1;
            N1=1;
            N2=N-N1;
            for CP=2:length(Data)
                LogPs(CP)=logFactorial(C1)+logFactorial(C2)-(C1+1)*log(N1)-(C2+1)*log(N2)...
                                - log((C1/N1)^2+(C2/N2)^2);
                C1=C1+Data(CP);
                C2=C2-Data(CP);
                N1=N1+1;
                N2=N2-1;
            end
        end

        function CP=estimateChangePointLocation(Data)
            % get the location of the change point
            LogPs=smi_stat.ChangeDetection.logP_ChangePoint(Data);
            [~,CP]=max(LogPs(2:end)); %the first time cannot be a change point
            CP=CP+1; %correct for dropping first point
        end

        function LogP=logPLambda(Data, Lambda)
            %compute log probability of the mean intensity (Lambda)
            N=length(Data);
            C=sum(Data);
            LogP=C*log(N) + (C-1)*log(Lambda) - N*Lambda - logFactorial(C-1);
        end

        function Lambda=estimateIntensity(Data)
            % compute the mean intensity of the data
            N=length(Data);
            C=sum(Data);
            Lambda=C/N;
        end

        function Model_intensity=modelIntensity(NObservations, ChangePoints, Intensity)
            % Helper for plotting an input intensity sequence
            %   NObservations - Scalar integer: length of data sequence
            %   ChangePoints - vector: indexs of change point locations
            %   intensity - vector: length=length(ChangePoints)+1: mean intensities
            %       for each subinterval
            % Outputs:
            %   Model_intensity -intensity of model sequence (no noise)
            Model_intensity=zeros(1,NObservations);
            nCP=length(ChangePoints);
            CP=1;
            for n=1:NObservations
                if CP<=nCP && ChangePoints(CP)==n
                    CP=CP+1;
                end
                Model_intensity(n)=Intensity(CP);
            end
        end


    end %static protected methods
end
