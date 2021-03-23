function [LPost,D]=testData(obj,M)
%testData Tests the quality of the raw input data
% [LPost,D]=obj.testData(M)
%
% analyze_all() is the main function and all the other functions
% are called inside this function. (genROIs, precluster, run 
% all clusters, generate all statistics) 
%
%    obj:   BaGoL object with an SMD structure as a field,  The SMD structure
%           contains the following fields:
%       X:        Vector of X localization positions (nm)(Nx1)
%       Y:        Vector of Y localization positions (nm)(Nx1)
%       X_SE:     Vector of X localization standard errors (nm)(Nx1)
%       Y_SE:     Vector of Y localization standard errors (nm)(Nx1)
%       FrameNum: Vector of localization frame numbers (Nx1)
%
% INPUTS:
%    M:     number of ROIs to test
%
% OUTPUTS:
%    LPost: cluster posteriors
%    D:     scaled distances 
%

% Created by:
%    Mohamadreza Fazel (LidkeLab 2020)

Ind = obj.SMD.X_SE == 0 | obj.SMD.Y_SE == 0;
obj.SMD.X(Ind) = [];
obj.SMD.Y(Ind) = [];
obj.SMD.X_SE(Ind) = [];
obj.SMD.Y_SE(Ind) = [];
obj.SMD.FrameNum(Ind) = [];
if isempty(obj.SMD.FrameNum)
    obj.SMD.FrameNum = zeros(size(obj.SMD.X));
else
    obj.SMD.FrameNum = single(obj.SMD.FrameNum);
end
ROIs = obj.genROIs();
ID = randi([1,length(ROIs(:))],[1,M]);
tROIs = ROIs(ID);
obj = obj.precluster(tROIs);

ClustNumHeirar = length(obj.ClusterSMD);
if isempty(obj.SMD.FrameNum)
    MaxAlpha = 0; 
elseif max(obj.SMD.FrameNum~=0)
    MaxAlpha = obj.Drift / max(obj.SMD.FrameNum); 
else
    MaxAlpha = 0; 
end
D=[];
LPost=[];
for nn = 1:ClustNumHeirar
    [TChain]=smi.BaGoL.BaGoL_RJMCMC(obj.ClusterSMD(nn),obj.Lambda,MaxAlpha,obj.P_Jumps,obj.N_Trials,obj.N_Burnin,0);
    N=[];
    for mm = 1:length(TChain)
        if TChain(nn).N == 0
            continue;
        end
        N = cat(1,N,length(TChain(mm).X));
    end
    Frequency = hist(N,(0:max(N)));
    Frequency(1) = [];
    MostFrequent = find(Frequency==max(Frequency));
    MostFrequent = MostFrequent(1);
    for ii = 1:length(TChain)
        Chain = TChain(ii);
        X = Chain.X;
        Y = Chain.Y;
        if length(X) ~= MostFrequent
           continue; 
        end
        NE = length(X);
        for pp = 1:NE
            Ind = Chain.ID==pp;
            NL=sum(Ind);
            %look at CDF of data 
            XDist = (obj.ClusterSMD(nn).X(Ind)-X(pp))./obj.ClusterSMD(nn).X_SE(Ind);
            YDist = (obj.ClusterSMD(nn).Y(Ind)-Y(pp))./obj.ClusterSMD(nn).Y_SE(Ind);
            D=cat(1,D,XDist,YDist);  
            LPost = cat(1,LPost,(sum(log(normpdf(obj.ClusterSMD(nn).X(Ind),X(pp)*ones(NL,1),obj.ClusterSMD(nn).X_SE(Ind))))...
                +sum(log(normpdf(obj.ClusterSMD(nn).Y(Ind),Y(pp)*ones(NL,1),obj.ClusterSMD(nn).Y_SE(Ind))))...
                +log(binopdf(NL,length(Chain.ID),1/NE)))/NL);
        end
        
    end
end

CDF=cumsum(ones(size(D)));
CDF=CDF/max(CDF);
Dsort = sort(D);

Ds=linspace(-4,4,100);
NormCDF=1/2*(1+erf(Ds/(sqrt(2))));

figure;
plot(Ds,NormCDF,'linewidth',3)
xlabel('Scaled Dist')
ylabel('CDF')
hold on
plot(Dsort,CDF,'--','linewidth',2)
legend('Theory','Data')

figure;hist(exp(LPost),50)
xlabel('Cluster Posteriors')
ylabel('Frequency')

obj.ClusterSMD=[];

end
