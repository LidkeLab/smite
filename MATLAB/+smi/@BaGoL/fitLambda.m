function [ParamPoiss,ParamGam]=fitLambda(NMean,Percent)
%fitLambda Fits Gamma and Poisson distributions to the input data
%
% To find the correct distribution of localizations per emitter, one way is
% running BaGoL with a very broad distribution for lambda and then use the
% returned values to compute the distribution of localizations per emitter.
% fitLambda() takes the list localizations per emitter found from an 
% initial run and fits that with both gamma and Poisson distributions. User
% can pick one of the returned distribution for a second run of BaGoL with
% a more informative prior on localizations per emitter. A normalized
% histogram of the inputs and plots of both distributions is also returned
% to help with picking an appropriate distribution.
%
% INPUT:
%    NMean: list of localizations per emitter (Kx1)
%    Percent: Percentage of small lambda values to be removed (Default = 0)
%
% OUTPUT:
%    ParamPoiss: Mean of fitted Poisson dist.
%    ParamGam: Structure containing shape and rate parameters of fitted gamma
%
% Created by:
%    Mohamadreza Fazel (Lidke lab 2020)
%

if nargin < 2 
   Percent = 0; 
end
L = quantile(NMean,Percent);
NMean = NMean(NMean>L);
ParamGam = fitdist(NMean,'gamma');
ParamPoiss = fitdist(NMean,'Poisson');

X = 0:max(NMean)+20;
figure;histogram(NMean,'normalization','pdf')
hold;plot(X,gampdf(X,ParamGam.a,ParamGam.b),'linewidth',1.5)
plot(X,poisspdf(X,ParamPoiss.lambda),'linewidth',1.5)
legend({'Lambda','Gamma dist.','Poisson dist.'},'FontSize',12)
xlabel('Lambda','FontSize',15)
ylabel('PDF','FontSize',15)

end
