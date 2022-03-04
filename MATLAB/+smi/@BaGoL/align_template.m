function [Alpha,XShift,YShift,Aligned,Chain] = align_template(Temp,Input,Start,Cutoff,NChain,PlotFlag)
%align_template Aligns a set of points to a template.
% [Alpha,XShift,YShift,Aligned,Chain] = BaGoL.align_template(Temp,Input,Start,Cutoff,NChain,PlotFlag)
%
% This method finds the shift and rotation of a set of points that minimizes
% the sum of distances of pairs from the sample to the template. Pair
% distances above a cutoff and unmatched points use the cost at the 
% cutoff distance. 
%
% INPUTS:
%   Temp:       SMD Structure of the template with fields
%       X:      X positions of the template (nm)(Nx1)
%       Y:      Y positions of the template (nm)(Nx1)
%   Input:      SMD Structure of the sample with fields
%       X:      X positions of the sample (nm)(Mx1)
%       Y:      Y positions of the sample (nm)(Mx1)
%   Start:      Structure of starting values for shift and rotation
%       Theta:  (radians)
%       ShiftX: (nm)
%       ShiftY: (nm)    
%    Cutoff:   Saturation distance for cost function (nm) 
%    NChain:   Number of jumps in MCMC algorithm. (Default = 5000)
%    PlotFlag: 0 or 1, Show an animation of the accepted jumps. (Default = 0)
%
% OUTPUTS:
%    Alpha:    The rotation angle. (radians)
%    XShift:   The translation along the X-axis. (nm)
%    YShift:   The translation along the Y-axis. (nm)
%    Aligned:  Structure containing the aligned coordinates
%       X:     (Mx1) (nm)
%       Y:     (Mx1) (nm)
%    Chain:    Structure array of accepted jumps in parameter space 
%       Theta:      Angle (radians) (NChain x 1)
%       XShift:     Sample shift in X (NChain x 1)
%       YShift:     Sample shift in Y (NChain x 1)

% Created by:
%    Mohamadreza Fazel (LidkeLab 2019)

if nargin < 4
   error('There must be at leat 4 inputs.') 
end

if nargin < 5
   NChain = 5000;
   PlotFlag = 0;
end

if nargin < 6
   PlotFlag = 0; 
end

Aligned.X = [];
Aligned.Y = [];

TempX = Temp.X - mean(Temp.X);
TempY = Temp.Y - mean(Temp.Y);

XLim = [min(TempX)-10 max(TempX)+10];
YLim = [min(TempY)-10 max(TempY)+10];

X = Input.X';
Y = Input.Y';
X = X - mean(X);
Y = Y - mean(Y);

Points = [X;Y];
[~,Dis] = knnsearch(Points',[TempX',TempY']);
Dis = abs(Dis);
Dis(Dis>Cutoff) = Cutoff;
Cost_Current = sum(Dis);

Start_Theta = Start.Theta;
Start_ShiftX = Start.ShiftX;
Start_ShiftY = Start.ShiftY;

Theta = zeros(NChain,1);
DelX = zeros(NChain,1);
DelY = zeros(NChain,1);

Theta(1) = Start_Theta;
DelX(1) = Start_ShiftX;
DelY(1) = Start_ShiftY;
if PlotFlag
    figure;hold;
    V1 = VideoWriter('Alignment','MPEG-4');
    V1.FrameRate = 10; 
    open(V1) 
end
for nn = 2:NChain
    
    if nn < NChain/2
        Theta_P = Theta(nn-1)+randn();
        DelX_P = DelX(nn-1) + 0.5*randn();
        DelY_P = DelY(nn-1) + 0.5*randn();
    else
        Theta_P = Theta(nn-1)+0.1*randn();
        DelX_P = DelX(nn-1) + 0.02*randn();
        DelY_P = DelY(nn-1) + 0.02*randn();
    end
    R = [cos(Theta_P),sin(Theta_P);-sin(Theta_P),cos(Theta_P)];
    
    Points = [X;Y];
    Points = R*Points;
    Points(1,:) = Points(1,:)+DelX_P;
    Points(2,:) = Points(2,:)+DelY_P;
    
    [~,Dis]=knnsearch(Points',[TempX',TempY']);
    Dis = abs(Dis);
    Dis(Dis>Cutoff) = 20;
    Cost_Proposed = sum(Dis);
    
    if Cost_Current - Cost_Proposed > -rand()
        if PlotFlag
            plot(TempX,TempY,'*')
            hold;
            plot(Points(1,:),Points(2,:),'*')
            title(sprintf('Step:%g',nn))
            xlim(XLim);ylim(YLim)
            SP = sprintf('Data: %g',nn);
            legend(SP)
            Frame = getframe(gcf);
            writeVideo(V1,Frame);
            hold off
        end
        Theta(nn) = Theta_P;
        DelX(nn) = DelX_P;
        DelY(nn) = DelY_P;
        Cost_Current = Cost_Proposed;
    else
        Theta(nn) = Theta(nn-1);
        DelX(nn) = DelX(nn-1);
        DelY(nn) = DelY(nn-1);
    end
    
end
if PlotFlag
    close(V1)
end
Chain.Theta = Theta;
Chain.XShift = DelX;
Chain.YShift = DelY;

Alpha = mean(Theta(round(NChain/2),:));
XShift = mean(DelX(round(NChain/2),:));
YShift = mean(DelY(round(NChain/2),:));

R = [cos(Alpha),sin(Alpha);-sin(Alpha),cos(Alpha)];
Points = [X;Y];
Points = R*Points;
Points(1,:) = Points(1,:)+XShift;
Points(2,:) = Points(2,:)+YShift;

Aligned.X = cat(1,Aligned.X,Points(1,:)'+mean(Temp.X));
Aligned.Y = cat(1,Aligned.Y,Points(2,:)'+mean(Temp.Y));

end
