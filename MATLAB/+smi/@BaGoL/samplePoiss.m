function Xi = samplePoiss(NPoints,K,Xi,Alpha,Beta)
%samplePoiss() infers parameter of Poisson prior on number of locas per emitters
%
%Here, the prior parameters are learned using a hierarchical Bayes skim.
%
% INPUT:
%   NPoints: Number of localizations within ROIs
%   K:       Number of found emitters within ROIs
%   Xi:      Current parameter values of the Prior
%   Alpha:   Shape parameter of Lambda hyper-prior (Default = 1)
%   Bets:    Scale parameter of Lambda hyper-prior (Default = 50)
%
% OUTPUT: 
%   Xi:  Updated mean for Poisson prior on number of emitters
%

% Created by:
%   Mohamadreza Fazel (Lidke lab, 2022)
% 

if nargin < 4
    Alpha = 1;
end
if nargin < 5
    Beta = 100;
end
Alpha_Prop = 3000;

%10 samples are taken in a row and the last one is returned
for ii = 1:10
    
    %Sample Lambda
    Xi_Prop = gamrnd(Alpha_Prop,Xi/Alpha_Prop);
    LogLikeR = sum(log(poisspdf(NPoints,K*Xi_Prop)) - log(poisspdf(NPoints,K*Xi)));
    LogPriorR = log(gampdf(Xi_Prop,Alpha,Beta)) - log(gampdf(Xi,Alpha,Beta));
    LogPropR = log(gampdf(Xi,Alpha_Prop,Xi_Prop/Alpha_Prop)) ...
        - log(gampdf(Xi_Prop,Alpha_Prop,Xi/Alpha_Prop));

    if LogLikeR + LogPriorR + LogPropR > log(rand())
       Xi = Xi_Prop; 
    end

end

end