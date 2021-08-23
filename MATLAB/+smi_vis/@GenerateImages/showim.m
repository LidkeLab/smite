function showim(Image)
%showim Displays a nicely formatted images of a grayscale or color image
%
%
%   INPUT
%      Image:   2D or 3D Array
%   OUTPUT
%
%   REQUIRES
%      Matlab 2014b or higher

%make 2^N bigger to make bigger that 256
SZ=max(size(Image,[1,2]));
Mag = 2^ceil(log(256/SZ)/log(2))

figure
imshow(Image,[],'border','tight','initialmagnification',Mag*100)
set(gca,'visible','off')
axis equal
axis tight

end

