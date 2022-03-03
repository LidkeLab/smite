function dispIm()
%dispIm GUI to display BaGoL output images with some useful tools.
% BaGoL.dispIm()
%
% When called, a GUI pops up and multiple images can be uploaded using the
% "Choose Image" button. The "Backward" and "Forward" buttons allow going
% back and forth between the loaded images. Tools from the menu bar, such as
% zoom in and zoom out can be used to adjust the displayed images. The
% adjustments are maintained when moving back and forth between the images.
%

% Created by:
%    Mohamadreza Fazel (LidkeLab 2019)

Im={};
map={};
Num = 1;
guiFig = figure('Visible','on','Position',[250,100,710,800],'Name',...
    'Image display','NumberTitle','off','Interruptible','off');

HAxis1 = axes('Units','normalize','YTick',[],'XTick',[],'Position',[.05,0.02,0.9,0.88]);

fileButton1 = uicontrol('Style', 'pushbutton','String','Choose Image',...
    'Unit','normalize','Enable','on','Position', [0.1 0.92 0.15 0.06],'Callback',@getFileList1);

 ForwardButton = uicontrol('Style', 'pushbutton','String','Forward',...
    'Unit','normalize','Enable','on','Position', [0.5 0.92 0.15 0.06],'Callback',@forward);
 BackwardButton = uicontrol('Style', 'pushbutton','String','Backward',...
    'Unit','normalize','Enable','on','Position', [0.75 0.92 0.15 0.06],'Callback',@backward);

function getFileList1(~,~) % Choose data file(s) and list them in the File List box
[filename, pathname]=uigetfile('Y:\*.png;*.tif','MultiSelect','on');
Im = cell(1,length(filename));
map = Im;
for ii = 1:length(filename)
[ImT,MapT] = imread(fullfile(pathname,filename{ii}));
Im{ii}=ImT;
map{ii}=MapT;
end
axes(HAxis1);imshow(Im{1},map{1});
end

function forward(~,~)
Num = Num+1;
if Num > length(Im)
    Num = 1; 
end
HAxis1.Children.CData = Im{Num};
end

function backward(~,~)
        Num = Num-1;
if Num == 0
    Num = length(Im);
end
HAxis1.Children.CData = Im{Num};
end
end
