% ==============================================
% function ReadImage
% ==============================================
function [file_name,Img,winsize,winstepsize,gridxy,LoadImgMethod] =ReadImage(varargin)

fprintf('Choose method to load images:  \n')
fprintf('     0: Select images folder;  \n')
fprintf('     1: Use prefix of image names;  \n')
fprintf('     2: Manually select images.  \n')
prompt = 'Input here: ';
LoadImgMethod = input(prompt);

switch LoadImgMethod 
    case 0
        % ==============================================
        imgfoldername = uigetdir(pwd,'Select images folder');
        img1 = dir(fullfile(imgfoldername,'*.jpg'));
        img2 = dir(fullfile(imgfoldername,'*.jpeg'));
        img3 = dir(fullfile(imgfoldername,'*.tif'));
        img4 = dir(fullfile(imgfoldername,'*.tiff'));
        img5 = dir(fullfile(imgfoldername,'*.bmp'));
        img6 = dir(fullfile(imgfoldername,'*.png'));
        img7 = dir(fullfile(imgfoldername,'*.jp2'));
        file_name = [img1;img2;img3;img4;img5;img6;img7];
        file_name = struct2cell(file_name);
    case 1
        % ==============================================
        fprintf('What is prefix of DIC images? E.g. img_0*.tif.   \n')
        prompt = 'Input here: ';
        file_name = input(prompt,'s');
        [~,imgname,imgext] = fileparts(file_name);
        file_name = dir([imgname,imgext]);
        file_name = struct2cell(file_name);
    otherwise
        % ==============================================
        disp('--- Please load First image ---')
        file_name{1,1} = uigetfile('*.tif','Select reference Image (Deformed)');
        disp('--- Please load Second image ---')
        file_name{1,2} = uigetfile('*.tif','Select deformed Image (Reference)');
        prompt = 'Do you want to load more deformed images? (0-yes; 1-no)';
        DoYouWantToLoadMoreImages = input(prompt); imageNo = 2;
        while ( DoYouWantToLoadMoreImages == 0 )   
            imageNo = imageNo + 1;
            file_name{1,imageNo} = uigetfile('*.tif','Select Deformed Image');
            prompt = 'Do you want to load more deformed images? (0-yes; 1-no)';
            DoYouWantToLoadMoreImages = input(prompt);
        end
end

% ==============================================
% The following codes only consider two images comparasion
numImages = size(file_name,2);
for i = 1:numImages
    Img{i} = imread(file_name{1,i});
    % Change color RGB images to grayscale images
    [~, ~, numberOfColorChannels] = size(Img{i});
    if (numberOfColorChannels==3)
        Img{i} = rgb2gray(Img{i});
    end
    Img{i} = double(Img{i})';
end
%f = imread(file_name{1});
%g = imread(file_name{2});


%if (numberOfColorChannels==3)
%    f = rgb2gray(f);
%    g = rgb2gray(g);
%end

% % Field of the undeformed and deformed greyscale;
%f = double(f'); g = double(g');

% images are:
% --------------------
% --  This is image --
% |                  |
% y                  |
% |                  |
% |  --> x direction |
% |                  |
% --------------------
% after transforming,  MatLab matrix direction:
% --  This is matrix in Matlab --
% |                             |
% x                             |
% |                             |
% |  --> y direction            |
% |                             |
% --------------------------------
 
% ==============================================
% Choose ZOI
disp('--- Choose ZOI two boundary points from the left-top to the right-bottom ---')
imshow( (imread(file_name{1})));

gridx = zeros(1,2); gridy = zeros(1,2);
[gridx(1), gridy(1)] = ginput(1);
fprintf('The left-top coordinates are (%4.3f,%4.3f)\n',gridx(1), gridy(1))

[gridx(2), gridy(2)] = ginput(1);
fprintf('The right-bottom coordinates are (%4.3f,%4.3f)\n',gridx(2), gridy(2))

gridxy.gridx = round(gridx); gridxy.gridy = round(gridy);

% Choose subset size
prompt = '--- What is the subset size? --- Input here: ';
winsize = input(prompt);

% Choose subset size
prompt = '--- What is the subset step? --- Input here: ';
winstepsize = input(prompt);
 


end