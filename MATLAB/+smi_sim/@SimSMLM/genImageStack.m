function [Model,Data]=genImageStack(obj)

%generate images

% To generate the blobs without noise, execute the following:
[Model] = smi_sim.GaussBlobs.gaussBlobImage(obj.SZ,obj.NFrames,obj.SMD_Model,0,0,0);
Model = Model+obj.Bg;

if nargout>1
    Data = poissrnd(single(Model));
end

end
