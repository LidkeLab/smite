function Xi = sampleGam(NPoints,K,Xi,Alpha,Beta)
%sampleGam() infers parameters of gamma prior on number of locs per emitters
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
%   Xi:  Updated shape and scale parameters for gamma prior on number
%            of emitters

% Created by:
%   Mohamadreza Fazel (Lidke lab, 2022)
% 

Eta = Xi(1);
Gamma = Xi(2);

if nargin < 4
    Alpha = 1;
end
if nargin < 5
    Beta = 50;
end
Alpha_Prop = 3000;

Ind = NPoints > 0;
NPoints = NPoints(Ind);
K = K(Ind);

%10 samples are taken in a row and only the last one is returned
for ii = 1:10
    
    %Sample Eta
    Eta_Prop = gamrnd(Alpha_Prop,Eta/Alpha_Prop);
    
    LogLikeR1 = log(gampdf(NPoints,K*Eta_Prop,Gamma)); 
    LogLikeR2 = log(gampdf(NPoints,K*Eta,Gamma));
    if sum(isinf(LogLikeR1) & isinf(LogLikeR2)) 
        Ind = isinf(LogLikeR1) & isinf(LogLikeR2);
        LogLikeR1(Ind) = 0;
        LogLikeR2(Ind) = 0;
    end
    LogLikeR2(isinf(LogLikeR1)) = 0;
    LogLikeR1(isinf(LogLikeR1)) = 0;
    LogLikeR = sum(LogLikeR1 - LogLikeR2);
    
    LogPriorR = log(gampdf(Eta_Prop,Alpha,Beta)) - log(gampdf(Eta,Alpha,Beta));
    LogPropR = log(gampdf(Eta,Alpha_Prop,Eta_Prop/Alpha_Prop)) ...
        - log(gampdf(Eta_Prop,Alpha_Prop,Eta/Alpha_Prop));

    if LogLikeR + LogPriorR + LogPropR > log(rand())
       Eta = Eta_Prop; 
    end

    %Sample Gamma
    Gamma_Prop = gamrnd(Alpha_Prop,Gamma/Alpha_Prop);
    LogLikeR1 = log(gampdf(NPoints,K*Eta,Gamma_Prop)); 
    LogLikeR2 = log(gampdf(NPoints,K*Eta,Gamma));
    if sum(isinf(LogLikeR1) & isinf(LogLikeR2)) 
        Ind = isinf(LogLikeR1) & isinf(LogLikeR2);
        LogLikeR1(Ind) = 0;
        LogLikeR2(Ind) = 0;
    end
    LogLikeR2(isinf(LogLikeR1)) = 0;
    LogLikeR1(isinf(LogLikeR1)) = 0;
    LogLikeR = sum(LogLikeR1 - LogLikeR2);
    LogPriorR = log(gampdf(Gamma_Prop,Alpha,Beta)) - log(gampdf(Gamma,Alpha,Beta));
    LogPropR = log(gampdf(Gamma,Alpha_Prop,Gamma_Prop/Alpha_Prop)) ...
        - log(gampdf(Gamma_Prop,Alpha_Prop,Gamma/Alpha_Prop));

    if LogLikeR + LogPriorR + LogPropR > log(rand())
       Gamma = Gamma_Prop; 
    end

end
Xi(1) = Eta;
Xi(2) = Gamma;

end