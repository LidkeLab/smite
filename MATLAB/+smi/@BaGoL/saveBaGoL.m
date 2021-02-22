function saveBaGoL(obj,L,SaveDir,OverlayFlag)
%Saves plots of NND, precisions, lambda, SR, Filter-SR, MAPN, Posterior, Overlay images
%
%The first plot is the histogram of the nearest neighbor distances of MAPN
%coordinates. The second and third plots are, respectively, the histograms
%of MAPN X-coordinates and MAPN Y-coordinates. The fourth plot is the
%histogram of the number of combined localizations for each emitter.
%The first image is the image reconstructed using the MAPN coordinates. The
%second and third images are, respectively, reconstructed using the raw
%input coordinates before and after NND-filter. The fourth image is the
%landscape image of the posterior. The last image is the ovelay of filtered
%SR-image and the posterior. Finally, the MAPN-coordinates are saved in a
%mat-file.
%
%INPUT:
%   SaveDir: Save directory of plots and images (optional)
%   L: Length of scale bars on generated images (nm) (Default=100) (optional)
%   SaveDir: Saving directory (optional)
%   OverlayFlag: binary parameter that makes overlay color images (Default=0) 
%
%OUTPUT:
%   NONE
%
%Created by:
%   Mohamadreza Fazel (Lidke Lab, 2020)
%

if nargin < 4
   OverlayFlag = 0; 
end
if nargin < 3
    if  ~isdir('Result_BaGoL')
       mkdir('Result_BaGoL'); 
    end
    SaveDir = 'Result_BaGoL';
end

if nargin>1
   Length = L; %nm 
else
   Length = 100; %nm 
end

%Saving NND-plot
NBins=30;
[~,Dis]=knnsearch([obj.MAPN.X,obj.MAPN.Y],[obj.MAPN.X,obj.MAPN.Y],'k',2);
Dis = Dis(:,2);
P = prctile(Dis,99);
figure('Visible','off')
hist(Dis(Dis<P),NBins)
xlabel('NND(nm)','FontSize',18)
ylabel('Frequency','FontSize',18)
print(gcf,fullfile(SaveDir,'NND'),'-dpng')

%Saving precisions-plots
X_SE = obj.MAPN.X_SE;
figure('Visible','off')
P = prctile(X_SE,99);
hist(X_SE(X_SE<P),NBins)
xlabel('X-SE(nm)','FontSize',18)
ylabel('Frequency','FontSize',18)
print(gcf,fullfile(SaveDir,'BaGoL_X-SE'),'-dpng')

Y_SE = obj.MAPN.Y_SE;
figure('Visible','off')
P = prctile(Y_SE,99);
hist(Y_SE(Y_SE<P),NBins)
xlabel('Y-SE(nm)','FontSize',18)
ylabel('Frequency','FontSize',18)
print(gcf,fullfile(SaveDir,'BaGoL_Y-SE'),'-dpng')

%Saving lambda
Nmean = obj.MAPN.Nmean;
P = prctile(Nmean,99);
figure('Visible','off')
histogram(Nmean(Nmean<P),0:max(Nmean)+15,'Normalization','pdf')
hold;
if length(obj.Lambda)>1
    Xp = 0:0.2:max(Nmean)+15;
    plot(Xp,gampdf(Xp,obj.Lambda(1),obj.Lambda(2)),'r','linewidth',1.5)
else
    Xp = 0:ceil(max(Nmean)+15);
    plot(Xp,poisspdf(Xp,obj.Lambda),'r','linewidth',1.5)
end
legend({'Found \lambda','Used dist.'},'FontSize',15)
xlabel('\lambda','FontSize',18);ylabel('PDF','FontSize',18)
print(gcf,fullfile(SaveDir,'Lambda'),'-dpng')

%Saving MAPN-image
ImFlag = 1;
PixelSize = obj.PixelSize;
[MapIm]=obj.genMAPNIm(ImFlag);
MapIm = smi.BaGoL.scaleIm(MapIm,98);
tMapIm = MapIm;
MapIm = smi.BaGoL.scalebar(MapIm,PixelSize,Length);
imwrite(MapIm,hot(256),fullfile(SaveDir,'MAPN-Im.png'));

%Saving SR-image
ImFlag = 2;
[SRIm]=obj.genMAPNIm(ImFlag);
SRIm = smi.BaGoL.scaleIm(SRIm,98);
SRIm = smi.BaGoL.scalebar(SRIm,PixelSize,Length);
imwrite(SRIm,hot(256),fullfile(SaveDir,'SR-Im.png'));

%Filter SR-image
[SM,Ind] = obj.NNDfilter(obj.SMD);
SMD.X = SM(:,1);
SMD.Y = SM(:,2);
SMD.X_SE = obj.SMD.X_SE(Ind);
SMD.Y_SE = obj.SMD.Y_SE(Ind);
SRImFilt = smi.BaGoL.makeIm(SMD,obj.PImageSize,obj.PixelSize,obj.XStart,obj.YStart);    
SRImFilt = smi.BaGoL.scaleIm(SRImFilt,98);
tSRImFilt = SRImFilt;
SRImFilt = smi.BaGoL.scalebar(SRImFilt,PixelSize,Length);
imwrite(SRImFilt,hot(256),fullfile(SaveDir,'SR-Im-Filter.png'));

%Saving posterior image
if obj.PImageFlag == 1
    PIm = smi.BaGoL.scaleIm(obj.PImage,96);
    tPIm = PIm;
    PIm = smi.BaGoL.scalebar(PIm,PixelSize,Length);
    imwrite(PIm,hot(256),fullfile(SaveDir,'Post-Im.png'));
end

%Overlay of filtered SR image and posterior image
if OverlayFlag
    Scale = 255;  
    overlayIm = zeros([size(tPIm),3]);                  
    overlayIm(:,:,1) = (tPIm/Scale+1.5*tSRImFilt/Scale)/2.5;                  
    overlayIm(:,:,2) = 1.5*tSRImFilt/Scale/2.5;
    overlayIm(:,:,3) = 1.5*tSRImFilt/Scale/2.5;
    overlayIm(:,:,1) = smi.BaGoL.scalebar(overlayIm(:,:,1)*Scale,PixelSize,Length);
    overlayIm(:,:,2) = smi.BaGoL.scalebar(overlayIm(:,:,2)*Scale,PixelSize,Length);
    overlayIm(:,:,3) = smi.BaGoL.scalebar(overlayIm(:,:,3)*Scale,PixelSize,Length);
    overlayIm = overlayIm/Scale;
    imwrite(overlayIm, fullfile(SaveDir,'Overlay_SR_Post.png'), 'PNG');                  

    %Overlay of filtered SR image and MAPN image image
    overlayIm = zeros([size(tMapIm),3]);                  
    overlayIm(:,:,1) = (tMapIm/Scale+1.5*tSRImFilt/Scale)/2.5;                  
    overlayIm(:,:,2) = 1.5*tSRImFilt/Scale/2.5;
    overlayIm(:,:,3) = 1.5*tSRImFilt/Scale/2.5;
    overlayIm(:,:,1) = smi.BaGoL.scalebar(overlayIm(:,:,1)*Scale,PixelSize,Length);
    overlayIm(:,:,2) = smi.BaGoL.scalebar(overlayIm(:,:,2)*Scale,PixelSize,Length);
    overlayIm(:,:,3) = smi.BaGoL.scalebar(overlayIm(:,:,3)*Scale,PixelSize,Length);
    overlayIm = overlayIm/Scale;
    imwrite(overlayIm, fullfile(SaveDir,'Overlay_SR_Map.png'), 'PNG'); 
end

MAPN = obj.MAPN;
save(fullfile(SaveDir,'MAPN'),'MAPN')

end
