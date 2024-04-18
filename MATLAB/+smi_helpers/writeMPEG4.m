function writeMPEG4(ResultsDir, FileBaseName, RGBImage)
%writeMPEG4 writes a .mp4 image constructed from the RGB image input on all OS.
%
% writeMPEG4 uses VideoWriter to generate video files.  Note that
% VideoWriter in Linux cannot generate .mp4 files, so it necessary to
% save in a different format and convert to .mp4 via external software
% (ffmpeg), which, of course, must be installed.  ffmpeg can convert
% between a variety of video formats.
%
% INPUT:
%    ResultsDir     output directory
%    FileBaseName   output filename with no extension
%    RGBImage       RGB image stack (4 dimensions: R, G, B, time sequence)
%
% OUTPUT:
%    ResultsDir/FileBaseName.mp4 video file
%
% REQUIREMENTS:
%    ffmpeg installed on Linux systems (https://ffmpeg.org)

% Created by
%    Michael J. Wester (Lidke Lab, 2022)

   islinux = isunix && ~ismac;
   RGBout = fullfile(ResultsDir, FileBaseName);
   RGBImageReordered = permute(RGBImage, [1, 2, 4, 3]);
   % Be sure RGBImageReordered is properly scaled if it is a float array.
   maxRGBImageReordered = max(max(max(max(RGBImageReordered(:)))));
   if isfloat(RGBImageReordered) & maxRGBImageReordered > 1
      RGBImageReordered = RGBImageReordered ./ maxRGBImageReordered;
   end
   % VideoWriter in Linux cannot generate .mp4 files, so it necessary to
   % save in a different format and convert to .mp4 via external software
   % (ffmpeg), which, of course, must be installed.  ffmpeg can convert
   % between a variety of video formats.
   if ~islinux
      v = VideoWriter(RGBout, 'MPEG-4');
   else
      v = VideoWriter(RGBout, 'Motion JPEG AVI');
   end
   open(v);
   writeVideo(v, RGBImageReordered);
   close(v);
   if islinux
      cmd = sprintf('ffmpeg -y -i "%s.avi" "%s.mp4"', RGBout, RGBout);
      [status, result] = system(cmd);
      if status ~= 0
         result
      end
   end

end
