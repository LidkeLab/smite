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
   
    for nn = 1:length(TChain)
        if TChain(nn).N == MostFrequent
            NSum = zeros(MostFrequent,1);
            XChain = cat(1,XChain,TChain(nn).X);
            YChain = cat(1,YChain,TChain(nn).Y);
            AlphaX_Chain = cat(1,AlphaX_Chain,TChain(nn).AlphaX);
            AlphaY_Chain = cat(1,AlphaY_Chain,TChain(nn).AlphaY);
            for pp = 1: MostFrequent(1)
                NSum(pp) = sum(single(TChain(nn).ID)==pp); 
            end
            NChain = cat(1,NChain,NSum);
        end
    end
    XChain = XChain(:);
    YChain = YChain(:);
    Points(:,1)=XChain;
    Points(:,2)=YChain;
    if ~isempty(obj.SMD.Z)
        Points(:,3)=ZChain;                     
    end
    Iters = ceil(MostFrequent/3)+1;
    NSeed = randi([0,length(XChain)/MostFrequent-1],[Iters,1]);
    KSeed = zeros(MostFrequent,2,Iters);
    for jj = 1:Iters
        Nt = NSeed(jj);
        KSeed(:,:,jj) = [XChain(Nt*MostFrequent+1:(Nt+1)*MostFrequent),YChain(Nt*MostFrequent+1:(Nt+1)*MostFrequent)];
    end
    ID = kmeans(Points,MostFrequent,'Replicates',Iters,'MaxIter',300);
    
    Map.X = [];
    Map.Y = [];
    Map.X_SE = [];
    Map.Y_SE = [];
    Map.AlphaX = [];
    Map.AlphaY = [];
    Map.AlphaX_SE = [];
    Map.AlphaY_SE = [];
    Map.Nmean = [];
    for nn = 1:MostFrequent
        TX = mean(XChain(ID==nn)); 
        TY = mean(YChain(ID==nn));
        TX_SE = std(XChain(ID==nn));
        TY_SE = std(YChain(ID==nn));
        TAlphaX = mean(AlphaX_Chain(ID==nn));
        TAlphaY = mean(AlphaY_Chain(ID==nn));
        TAlphaX_SE = std(AlphaX_Chain(ID==nn));
        TAlphaY_SE = std(AlphaY_Chain(ID==nn));
        TNmean = mean(NChain(ID==nn));
        Map.X = cat(1,Map.X,TX);
        Map.Y = cat(1,Map.Y,TY);
        Map.X_SE = cat(1,Map.X_SE,TX_SE);
        Map.Y_SE = cat(1,Map.Y_SE,TY_SE);
        Map.AlphaX = cat(1,Map.AlphaX,TAlphaX);
        Map.AlphaY = cat(1,Map.AlphaY,TAlphaY);
        Map.AlphaX_SE = cat(1,Map.AlphaX_SE,TAlphaX_SE);
        Map.AlphaY_SE = cat(1,Map.AlphaY_SE,TAlphaY_SE);
        Map.Nmean = cat(1,Map.Nmean,TNmean);
    end
    Ind = obj.removeOverlap(ROIs,Map.X,Map.Y,ii);
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
    obj.MAPN.Nmean = cat(1,obj.MAPN.Nmean,Map.Nmean(Ind));
end
