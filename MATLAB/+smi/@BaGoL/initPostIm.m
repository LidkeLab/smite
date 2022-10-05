function initPostIm(obj)
%initPostIm Find the Posterior Image size and position

DX = obj.PixelSize;
obj.XStart = min(obj.SMD.X-3*obj.SMD.X_SE);
obj.YStart = min(obj.SMD.Y-3*obj.SMD.Y_SE);

X_max = max(obj.SMD.X+3*obj.SMD.X_SE);
Y_max = max(obj.SMD.Y+3*obj.SMD.Y_SE);

obj.PImageSize=ceil(max((X_max-obj.XStart)/DX,(Y_max-obj.YStart)/DX));

end