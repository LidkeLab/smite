function genMAPN(obj,Chain,ROIs,ii)
%genMAPN Generates the MAPN coordinates from the most repeated model in the chain.
% obj = obj.genMAPN(Chain,ROIs,ii)
%
% The states of the chain with most repeated number of emitters (model)
% are identified. The coordinates from the identified states are
% extracted and analyzed with the kmeans clustering algorithm to find the
% clusters within the retrieved coordinates. The center of mass of the
% found clusters and their standard deviations are taken as the emitter
% locations and precisions, respectively. The kmeans algorithm is run with
% multiple sets of seeds, which depends on the number of emitters, and the
% configuration with the least cost is used. For some of the seeds, kmeans
% may fail to converge and throws a warning which is not an issue. The
% locations inside the overlapping regions (see genROIs()) are then
% eliminated from the list of emitter coordinates and the reminder is
% returned.
%
% INPUTS:
%    Chain: RJMCMC chain for a cluster (See BaGoL_RJMCMC())
%    ROIs:  2D Array of SMD structures (see genROIs())
%    ii:    Index of cluster associated to the input chain
%
% OUTPUTS:
%    NONE
%

% Created by:
%    Mohamadreza Fazel (LidkeLab 2019)

TChain = Chain;
N = [];
XChain = [];
YChain = [];
AlphaX_Chain = [];
AlphaY_Chain = [];
NChain = [];
for nn = 1:length(TChain)
    if TChain(nn).N == 0
        continue;
    end
    N = cat(1,N,length(TChain(nn).X));
end
Frequency = hist(N,(0:max(N)));
Frequency(1) = [];
MostFrequent = find(Frequency==max(Frequency));
MostFrequent = MostFrequent(1);

% Get last positions from chain with N = MostFrequent
for nn = 1:length(TChain)
    if TChain(nn).N == MostFrequent
        XChain = TChain(nn).X;
        YChain = TChain(nn).Y;
        AlphaX_Chain = TChain(nn).AlphaX;
        AlphaY_Chain = TChain(nn).AlphaY;
    end
end

% Run MCMC chain (RJCMCMC without Add/Remove).

if isempty(obj.SMD.FrameNum)
    MaxAlpha = 0;
elseif max(obj.SMD.FrameNum~=0)
    MaxAlpha = obj.Drift;
else
    MaxAlpha = 0;
end

P_Jumps_MCMC = [0.5 0.5 0 0];
ChainLength = length(obj.XiChain);

Xi = mean(obj.XiChain(round(ChainLength/2)+1:end,:));

N_Burnin = 500;
N_Trials = 1000;
if MostFrequent > 20
    AnimFlag = 0;
else
    AnimFlag = 0;
end

[TChain]=smi.BaGoL.BaGoL_RJMCMC(obj.ClusterSMD(ii),Xi,MaxAlpha,P_Jumps_MCMC,N_Trials,N_Burnin, AnimFlag, ...
    XChain,YChain,AlphaX_Chain, AlphaY_Chain);

% Cut burnin
TChain = TChain(end-N_Trials+1:end)
TailSZ = 10;
MAPN =  [mode([TChain(end-TailSZ:end).X],2), mean([TChain(end-TailSZ:end).Y],2)];

if 0 % debug
    figure; plot([TChain(:).X]', [TChain(:).Y]','+')
    hold on
    plot(MAPN(:,1)', MAPN(:,2)', 'kx','MarkerSize',16,'LineWidth',3)
    title('MCMC Chain Raw Mode Tail')
end

% Unmix emitter IDs
N = MostFrequent;
CostMatrix = zeros(N,N);
ChainUnmixed = TChain;
MAPN =  [mode([ChainUnmixed(end-TailSZ:end).X],2), mean([ChainUnmixed(end-TailSZ:end).Y],2)];
% Do a first pass with aligning using MSE to everything before.
for kk = 1:length(ChainUnmixed)
    if 1 % align to mode of tail of chain
        X = MAPN(:,1);
        Y = MAPN(:,2);
    else
        X = [ChainUnmixed(1:kk-1).X];
        Y = [ChainUnmixed(1:kk-1).Y];
    end
    for jj = 1:N
        CostMatrix(:,jj) = mean( (X - ChainUnmixed(kk).X(jj)).^2 ,2) + ...
            mean( (Y - ChainUnmixed(kk).Y(jj)).^2 ,2);
    end
    [AssignIDs, Cost] = smi.SPT.solveLAP(CostMatrix);
    ChainUnmixed(kk).X(AssignIDs) = ChainUnmixed(kk).X;
    ChainUnmixed(kk).Y(AssignIDs) = ChainUnmixed(kk).Y;
    ChainUnmixed(kk).AlphaX(AssignIDs) = ChainUnmixed(kk).AlphaX;
    ChainUnmixed(kk).AlphaY(AssignIDs) = ChainUnmixed(kk).AlphaY;
end

if 0 % Debug
    figure; plot([ChainUnmixed(:).X]', [ChainUnmixed(:).Y]','+')
    hold on
    MAPN =  [mean([ChainUnmixed(:).X],2), mean([ChainUnmixed(:).Y],2)];
    plot(MAPN(:,1), MAPN(:,2), 'rx','MarkerSize',16,'LineWidth',3)
    title('Unmixed from Mixed')
end


Map.X = mean([ChainUnmixed(:).X],2);
Map.Y = mean([ChainUnmixed(:).Y],2);
Map.X_SE = std([ChainUnmixed(:).X],0,2);
Map.Y_SE = std([ChainUnmixed(:).Y],0,2);
Map.AlphaX = mean([ChainUnmixed(:).AlphaX],2);
Map.AlphaY = mean([ChainUnmixed(:).AlphaY],2);
Map.AlphaX_SE = std([ChainUnmixed(:).AlphaX],0,2);
Map.AlphaY_SE = std([ChainUnmixed(:).AlphaY],0,2);

% Remove overlap with other boxes by truncation of points outside box
Ind = obj.removeOverlap(ROIs,Map.X,Map.Y,ii);

% Add these new MAPN values to entire object MAPN
obj.MAPN.X = cat(1,obj.MAPN.X,Map.X(Ind));
obj.MAPN.Y = cat(1,obj.MAPN.Y,Map.Y(Ind));
obj.MAPN.X_SE = cat(1,obj.MAPN.X_SE,Map.X_SE(Ind));
obj.MAPN.Y_SE = cat(1,obj.MAPN.Y_SE,Map.Y_SE(Ind));
obj.MAPN.AlphaX = cat(1,obj.MAPN.AlphaX,Map.AlphaX(Ind));
obj.MAPN.AlphaY = cat(1,obj.MAPN.AlphaY,Map.AlphaY(Ind));
obj.MAPN.AlphaX_SE = cat(1,obj.MAPN.AlphaX_SE,Map.AlphaX_SE(Ind));
obj.MAPN.AlphaY_SE = cat(1,obj.MAPN.AlphaY_SE,Map.AlphaY_SE(Ind));
if ~isempty(obj.SMD.Z)
    obj.MAPN.Z = cat(1,obj.MAPN.Z,Map.Z(Ind));
    obj.MAPN.Z_SE = cat(1,obj.MAPN.Z_SE,Map.Z_SE(Ind));
end
obj.MAPN.Nmean = cat(1,obj.MAPN.Nmean,obj.N_Trials*ones(MostFrequent,1));
end
