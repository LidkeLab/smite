function RGBImage=colorOverlay(obj)
% colorOverlay Displays an overlay of the (green) model with (red) data
% 
% Usage: 
% myObj.colorOverlay()    diplays a RBG Stack in a window
% RGBImage = myObj.colorOverlay()  returns an scaled RGB image 




if isempty(obj.SMF.Data.DataROI)
    SZ=size(obj.ScaledData,[1,2]);
    obj.SMF.Data.DataROI=[1 1 SZ SZ];
end
[Model] = smi_sim.GaussBlobs.gaussBlobImage(obj.SMD,obj.SMF,0,'SMF');

minData = prctile(obj.ScaledData(:),1);
maxData = prctile(obj.ScaledData(:),99);
minModel = prctile(Model(:),1);
maxModel = prctile(Model(:),99);

MinRange = min(minData,minModel);
MaxRange = max(maxData,maxModel);

% Deal with the rare case that can occur with very sparse data (such as the
% unitTest).
if MinRange == MaxRange
   minData = min(obj.ScaledData(:));
   maxData = max(obj.ScaledData(:));
   minModel = min(Model(:));
   maxModel = max(Model(:));

   MinRange = min(minData,minModel);
   MaxRange = max(maxData,maxModel);
end

R = (obj.ScaledData-MinRange)/(MaxRange-MinRange);
G = (Model-MinRange)/(MaxRange-MinRange);

RGBImage=cat(4,R,G,G*0);

if nargout <1
    F=figure;
    F.Name='Data(red) Model(green)';
    sliceViewer(RGBImage);
end

end
