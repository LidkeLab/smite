function saveBaGoL(obj,L,SaveDir,OverlayFlag)
%Saves plots of NND, precisions, Xi, SR, MAPN, Posterior, Overlay images
%
%The first plot is the histogram of the nearest neighbor distances of MAPN
%coordinates. The second and third plots are, respectively, the histograms
%of MAPN X-precisions and MAPN Y-precisions. The fourth plot is the
%histogram of the number of localizations allocated to each emitter.
%The first image is the image reconstructed using the MAPN coordinates. The
%second image is reconstructed using the raw input coordinates. 
%The third image is the landscape image of the posterior. The last image is 
%the ovelay of SR-image and the posterior. Finally, the 
% MAPN-coordinates are saved in a mat-file.
%
%INPUT:
%   SaveDir: Save directory of plots and images (optional)
%   L: Length of scale bars on generated images (nm) (Default=100) (optional)
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

try
MAPN = obj.MAPN;
save(fullfile(SaveDir,'MAPN'),'MAPN')
catch ME
fprintf('### PROBLEM with saving MAPN ###\n');
fprintf('%s\n', ME.identifier);
fprintf('%s\n', ME.message);
end

%Saving NND-plot
try
NBins=30;
[~,Dis]=knnsearch([obj.MAPN.X,obj.MAPN.Y],[obj.MAPN.X,obj.MAPN.Y],'k',2);
Dis = Dis(:,2);
P = prctile(Dis,99);
figure('Visible','off')
hist(Dis(Dis<P),NBins)
xlabel('NND(nm)','FontSize',18)
ylabel('Frequency','FontSize',18)
print(gcf,fullfile(SaveDir,'NND'),'-dpng')
catch ME
fprintf('### PROBLEM with saving NND ###\n');
fprintf('%s\n', ME.identifier);
fprintf('%s\n', ME.message);
end

%Saving precisions-plots
try
X_SE = obj.MAPN.X_SE;
figure('Visible','off')
P = prctile(X_SE,99);
hist(X_SE(X_SE<P),NBins)
xlabel('X-SE(nm)','FontSize',18)
ylabel('Frequency','FontSize',18)
print(gcf,fullfile(SaveDir,'BaGoL_X-SE'),'-dpng')
catch ME
fprintf('### PROBLEM with saving BaGoL_X-SE ###\n');
fprintf('%s\n', ME.identifier);
fprintf('%s\n', ME.message);
end

try
Y_SE = obj.MAPN.Y_SE;
figure('Visible','off')
P = prctile(Y_SE,99);
hist(Y_SE(Y_SE<P),NBins)
xlabel('Y-SE(nm)','FontSize',18)
ylabel('Frequency','FontSize',18)
print(gcf,fullfile(SaveDir,'BaGoL_Y-SE'),'-dpng')
catch ME
fprintf('### PROBLEM with saving BaGoL_Y-SE ###\n');
fprintf('%s\n', ME.identifier);
fprintf('%s\n', ME.message);
end

%Saving hierarchical parameters
try
if obj.HierarchFlag == 1
   LChain = floor(obj.N_Trials/obj.NSamples);
   if length(obj.Xi)==2
       Lambda = obj.XiChain(end-LChain:end,2).*obj.XiChain(end-LChain:end,1);
   else
        Lambda = obj.XiChain(end-LChain:end);
   end
   figure('Visible','off')
   histogram(Lambda,'normalization','pdf')
   xlabel('\xi');ylabel('pdf')
   xlim([0 max(Lambda)+20])
   print(gcf,fullfile(SaveDir,'Xi'),'-dpng')
else
    Nmean = obj.MAPN.Nmean;
    P = prctile(Nmean,99);
    figure('Visible','off')
    histogram(Nmean(Nmean<P),0:max(Nmean)+15,'Normalization','pdf')
    hold;
    if length(obj.Xi)>1
        Xp = 0:0.2:max(Nmean)+15;
        plot(Xp,gampdf(Xp,obj.Xi(1),obj.Xi(2)),'r','linewidth',1.5)
    else
        Xp = 0:ceil(max(Nmean)+15);
        plot(Xp,poisspdf(Xp,obj.Xi),'r','linewidth',1.5)
    end
    legend({'Found \xi','Used dist.'},'FontSize',15)
    xlabel('\xi','FontSize',18);ylabel('PDF','FontSize',18)
    print(gcf,fullfile(SaveDir,'Xi'),'-dpng')
end
catch ME
fprintf('### PROBLEM with saving Xi ###\n');
fprintf('%s\n', ME.identifier);
fprintf('%s\n', ME.message);
end

%Saving MAPN-image
try
ImFlag = 1;
PixelSize = obj.PixelSize;
[MapIm]=obj.genMAPNIm(ImFlag);
MapIm = smi.BaGoL.scaleIm(MapIm,98);
tMapIm = MapIm;
MapIm = smi.BaGoL.scalebar(MapIm,PixelSize,Length);
imwrite(MapIm,hot(256),fullfile(SaveDir,'MAPN-Im.png'));
catch ME
fprintf('### PROBLEM with saving MAPN-Im.png ###\n');
fprintf('%s\n', ME.identifier);
fprintf('%s\n', ME.message);
end

%Saving SR-image
try
ImFlag = 2;
[SRIm]=obj.genMAPNIm(ImFlag);
SRIm = smi.BaGoL.scaleIm(SRIm,98);
tSRIm = SRIm;
SRIm = smi.BaGoL.scalebar(SRIm,PixelSize,Length);
imwrite(SRIm,hot(256),fullfile(SaveDir,'SR-Im.png'));
catch ME
fprintf('### PROBLEM with saving SR-Im.png ###\n');
fprintf('%s\n', ME.identifier);
fprintf('%s\n', ME.message);
end

%Saving posterior image
try
if obj.PImageFlag == 1
    PIm = smi.BaGoL.scaleIm(obj.PImage,98);
    tPIm = PIm;
    PIm = smi.BaGoL.scalebar(PIm,PixelSize,Length);
    imwrite(PIm,hot(256),fullfile(SaveDir,'Post-Im.png'));
end
catch ME
fprintf('### PROBLEM with saving Post-Im.png ###\n');
fprintf('%s\n', ME.identifier);
fprintf('%s\n', ME.message);
end

% Color coding: SR (green), Posterior (blue), MAPN (red).
% When there are just 2 colors with one of them green, the blue and red are
% mixed together to make magenta for easier visualization.
if OverlayFlag
    Scale = 255;

    try
    %Overlay of filtered SR image (green) and posterior image (magenta)
    overlayIm = zeros([size(SRIm),3]);
    %overlayIm(:,:,1) = (2*tPIm/Scale+tSRIm/Scale)/3;
    overlayIm(:,:,1) = tPIm/Scale/3;
    overlayIm(:,:,2) = tSRIm/Scale/3;
    overlayIm(:,:,3) = tPIm/Scale/3;
    overlayIm(:,:,1) = smi.BaGoL.scalebar(overlayIm(:,:,1)*Scale,PixelSize,Length);
    overlayIm(:,:,2) = smi.BaGoL.scalebar(overlayIm(:,:,2)*Scale,PixelSize,Length);
    overlayIm(:,:,3) = smi.BaGoL.scalebar(overlayIm(:,:,3)*Scale,PixelSize,Length);
    imwrite(overlayIm, fullfile(SaveDir,'Overlay_gSR_mPost.png'), 'PNG');
    catch ME
    fprintf('### PROBLEM with saving Overlay_SR_Post.png ###\n');
    fprintf('%s\n', ME.identifier);
    fprintf('%s\n', ME.message);
    end

    try
    %Overlay of filtered SR image (green) and MAPN image (magenta)
    overlayIm = zeros([size(tMapIm),3]);
    %overlayIm(:,:,1) = (2*tMapIm/Scale+tSRIm/Scale)/3;
    overlayIm(:,:,1) = tMapIm/Scale/3;
    overlayIm(:,:,2) = tSRIm/Scale/3;
    overlayIm(:,:,3) = tMapIm/Scale/3;
    overlayIm(:,:,1) = smi.BaGoL.scalebar(overlayIm(:,:,1)*Scale,PixelSize,Length);
    overlayIm(:,:,2) = smi.BaGoL.scalebar(overlayIm(:,:,2)*Scale,PixelSize,Length);
    overlayIm(:,:,3) = smi.BaGoL.scalebar(overlayIm(:,:,3)*Scale,PixelSize,Length);
    %overlayIm = 10*overlayIm/Scale;
    imwrite(overlayIm, fullfile(SaveDir,'Overlay_gSR_mMap.png'), 'PNG');
    catch ME
    fprintf('### PROBLEM with saving Overlay_SR_Map.png ###\n');
    fprintf('%s\n', ME.identifier);
    fprintf('%s\n', ME.message);
    end

    try
    %Overlay of posterior image (cyan) and MAPN image (red)
    overlayIm = zeros([size(tMapIm),3]);
    %overlayIm(:,:,1) = (2*tMapIm/Scale+tPIm/Scale)/3;
    overlayIm(:,:,1) = tMapIm/Scale/3;
    overlayIm(:,:,2) = tPIm/Scale/3;
    overlayIm(:,:,3) = tPIm/Scale/3;
    overlayIm(:,:,1) = smi.BaGoL.scalebar(overlayIm(:,:,1)*Scale,PixelSize,Length);
    overlayIm(:,:,2) = smi.BaGoL.scalebar(overlayIm(:,:,2)*Scale,PixelSize,Length);
    overlayIm(:,:,3) = smi.BaGoL.scalebar(overlayIm(:,:,3)*Scale,PixelSize,Length);
    imwrite(overlayIm, fullfile(SaveDir,'Overlay_cPost_rMap.png'), 'PNG');
    catch ME
    fprintf('### PROBLEM with saving Overlay_Post_Map.png ###\n');
    fprintf('%s\n', ME.identifier);
    fprintf('%s\n', ME.message);
    end

    try
    %Overlay of SR image (green), posterior image (blue) and MAPN image (red)
    overlayIm = zeros([size(tSRIm),3]);
    overlayIm(:,:,1) = tMapIm/Scale/3;
    overlayIm(:,:,2) = tSRIm/Scale/3;
    overlayIm(:,:,3) = tPIm/Scale/3;
    overlayIm(:,:,1) = smi.BaGoL.scalebar(overlayIm(:,:,1)*Scale,PixelSize,Length);
    overlayIm(:,:,2) = smi.BaGoL.scalebar(overlayIm(:,:,2)*Scale,PixelSize,Length);
    overlayIm(:,:,3) = smi.BaGoL.scalebar(overlayIm(:,:,3)*Scale,PixelSize,Length);
    imwrite(overlayIm, fullfile(SaveDir,'Overlay_gSR_bPost_rMap.png'), 'PNG');
    catch ME
    fprintf('### PROBLEM with saving Overlay_SR_Post_Map.png ###\n');
    fprintf('%s\n', ME.identifier);
    fprintf('%s\n', ME.message);
    end
end

%MAPN = obj.MAPN;
%save(fullfile(SaveDir,'MAPN'),'MAPN')
close all;

end
