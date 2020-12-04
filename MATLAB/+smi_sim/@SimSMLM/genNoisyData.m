function [Data]=genNoisyData(obj,Model)
[Model] = Model+obj.Bg;
[Data] = poissrnd(single(Model));
NoisyImage = zeros(size(Data(:,:,1)));
[Data] = Data+randn(size(Data)).*repmat(NoisyImage,[1 1 obj.NFrames]);
end