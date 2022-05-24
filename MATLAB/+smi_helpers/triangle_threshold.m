function [Cutoff] = triangle_threshold(Data,Percentile,ShowFig)
%triangle_threshold Use the triangle method to find a threshold
%   Implements method by:
%   Zack GW, Rogers WE, Latt SA (1977), 
%   "Automatic measurement of sister chromatid exchange frequency", 
%   J. Histochem. Cytochem. 25 (7): 741–53
%
% INPUTS
%   Data: Array of any size
%   Percentile: Remove Percentile of extremum (1e-4)
%   ShowFig:    Show the triangle figure (false)
% OUTPUT
%   Cutoff:     Threhold value

if nargin <2
    Percentile=1e-4;
end

if nargin <3
    ShowFig=false;
end

Data=Data(:);

% Remove extremum
MinV = prctile(Data,Percentile/2);
MaxV = prctile(Data,100-Percentile/2);
Mask = (Data>MinV) & (Data<MaxV);
Data=Data(Mask);

% calculate histogram
[N,EDGES] = histcounts(Data,'BinMethod','fd');
DX=EDGES(2)-EDGES(1);
BinCenters = EDGES(1:end-1)+DX/2;


%find max
ID = find(N==max(N),true,'first');

%find if tail is left or right
if ID>length(N)/2 %peak on right, flip!
    N=fliplr(N);
    BinCenters=fliplr(BinCenters);
end

X=1:length(N);
N=N./max(N)*max(X); % for visualization 

%the line
XMax = find(N==max(N),true,'first');
XMin = find(N==min(N),true,'last');
YMax = N(XMax);
YMin = N(XMin);

Slope = (YMax-YMin)/(XMax-XMin);
Offset = YMax - XMax*Slope;
myline=@(X) Slope*X+Offset;

% distance to line is always proportional to height
h = 0;
Cutoff=0;
for x=X
    testh = myline(x)-N(x);
    if testh>h && x>XMax
        h=testh;
        Cutoff = BinCenters(x);
        IDCutoff=x;
    end
end

if ShowFig
    figure
    bar(X,N,1.0)
    hold on
    plot(X,myline(X),'r')
    Slope2 = tan(pi/2-atan(-Slope)); 
    plot([IDCutoff,max(X)],[N(IDCutoff),(max(X)-IDCutoff)*Slope2+N(IDCutoff)],'r')
    xlabel('Bin ID')
    ylabel('Normalized Counts')
    axis equal
    set(gca,'XTick',X)
    set(gca,'XTickLabel',BinCenters)
end

end

