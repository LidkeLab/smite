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

MinRange = min(min(obj.ScaledData(:),min(Model(:))));
MaxRange = max(max(obj.ScaledData(:),max(Model(:))));

R = (obj.ScaledData-MinRange)/(MaxRange-MinRange);
G = (Model-MinRange)/(MaxRange-MinRange);

RGBImage=cat(4,R,G,G*0);

if nargout <1
    F=figure
    F.Name='Data(red) Model(green)'
    sliceViewer(RGBImage)
end

end