function [Pvalue, X_Point, Sigma_Point] = singleLabelTest(X, Sigma, Sigma_Reg)
%singleLabelTest Tests if a cluster of points came from point source 
%   This function calculates a p-value for a cluster of points. 
%   The meaning of of the p-value is the probabilty that a more extreme set
%   of N points came from a point source, where N is the number of observed
%   points. 
%
% INPUTS:
%   X,Y  NxM Vectors of positions. N is number of particles, M is dimension
%   X_Sigma, Y_Sigma    NxM Vectors of position uncertainty (1 STD)
%   Sigma_Reg           1xM array of registration error (1 sigma)
%
% OUTPUTS:
%   Pvalue:             Probability of more exterme cluster
%   X_Point:            Weighted mean location value of point source
%   Sigma_Point:        Weighted uncertainty of point source

% Created by
%    Keith Lidke (2021)

N=size(X,1);
M=size(X,2);
DOF = M*N-M;

if DOF==0
    DOF=1;
end

%Make uncertainty larger due to registration error
Sigma = sqrt(Sigma.^2+repmat(Sigma_Reg.^2,[N,1]));

% MLE of center position:
X_Point = sum(X./Sigma.^2,1)./sum(1./Sigma.^2,1);
Sigma_Point = sqrt(1./sum(1./Sigma.^2,1));

%Likelihood at MLE
L=normpdf(X,repmat(X_Point,[N,1]),Sigma);

%Likelihood at Null
L0 = normpdf(X,X,Sigma);

%Calculate likelihood ratio:
R=-2*sum(sum(log(L./L0)));

X2_CDF=inline('gammainc(x/2,k/2)','k','x');
Pvalue=1-X2_CDF(DOF,R);

end
