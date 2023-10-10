% Image processing pipeline for Cluster Sequencing on Epi microscope

%% User configuration to read all images
clear
close all
clc
% Read image folder
folderPath = 'C:\Users\yujingsong\OneDrive\Desktop\3_aligned'; %Project folder name

%file name detection
uniqueName_bk="_000_";  
uniqueName2 = {'_645_','_590_','_525_','_445_'};
filetype='*.tif'; 
f=fullfile(folderPath,filetype);
d=dir(f);
filename=fullfile(folderPath,{d(:).name});
filename_bk=filename(contains(filename, uniqueName_bk));
filename = filename(arrayfun(@(x) any(contains(x, uniqueName2)), filename));
totImg=length(filename);

%crop factors
sizeX=700;   %Crop size, X is col     Default 1600
sizeY=1500;   %Y is row                Default 1500
crop=[750 250 (sizeX-1) (sizeY-1)];

% Read a 365 or 445 image and check crop
I = imread(filename{3});
I = imcrop(I,crop);
figure, imshow(I,[0,5000]);

%read background image
I_bk = imread(filename_bk{1});
I_bk_crop = imcrop(I_bk,crop);
%figure, imshow(I_bk,[]);

%% batch processing for just median filter and crop
%Create a result saving folder
cd(folderPath);
outputFolder='AlignOnly';
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end
cd(outputFolder);
for ind_Img=1:totImg
    Img = imread(filename{ind_Img});
    Img = imcrop(Img,crop);
    Img = medfilt2(Img, [2,2]);  %either [3,3] will be slightly better than [2,2]
    [~, imageName, imageExt] = fileparts(filename{ind_Img});
    imwrite(Img, [imageName, imageExt]);
    %imwrite(I2, [imageName,'_test', imageExt]);
end
cd(folderPath);
% consider do registration here or use imageJ to do registration

%% bk normalization
for i=1:totImg
    Img = imread(filename{i});
    %Img = imcrop(Img,crop)-I_bk_crop;
    Img = imcrop(Img,crop);
    bk(i)=mean(mean(Img));
end
corr_factor=max(bk)./bk;

%% batch processing for bk norm and median filter
%Create a result saving folder
cd(folderPath);
outputFolder='Test';
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end
cd(outputFolder);
for ind_Img=1:totImg
    Img = imread(filename{ind_Img});
    %Img = imcrop(Img,crop)-I_bk_crop;
    %Img = medfilt2(Img, [3,3]);  %either [3,3] will be slightly better than [2,2]
    Img2 = double(Img).*corr_factor(ind_Img);
    Img3=uint16(Img2);
    %Img3(BW2)=mean(mean(Img3));
    [~, imageName, imageExt] = fileparts(filename{ind_Img});
    imwrite(Img3, [imageName, imageExt]);
    %imwrite(I2, [imageName,'_test', imageExt]);
end
cd(folderPath);
% consider do registration here or use imageJ to do registration

%% Scatter detect and correction
clear
close all
clc
% Read image folder
folderPath = 'C:\Users\yujingsong\OneDrive\Desktop\230929_S0486_R1387_onepot\6_bk_corr_aligned'; %Project folder name

%file name detection
filetype='*.tif'; 
f=fullfile(folderPath,filetype);
d=dir(f);
filename=fullfile(folderPath,{d(:).name});
totImg=length(filename);
I = imread(filename{4});
figure, imshow(I,[0,5000]);

%% Thresholding and generating scatter removal mask
%Threshold ratio
ImgSetNum=4;
BW = imbinarize(I);
BW(:,:)=0;
for i=1:4
    I = imread(filename{i});
    T1 = graythresh(I/3);
    if T1<0.1
        T1=0.1;   %incase graythresh goes wrong
    end
    BW1 = imbinarize(I,T1);
    BWdfill = imfill(BW1, 'holes');
    %figure, imshow(BW1);
    se2 = strel('square', 2);
    se5= strel('disk',5);
    BW2=imdilate(imerode(BWdfill,se2),se5);
    BW = BW | BW2;
end

figure, imshow(BW);
I2=I;
I2(BW)=mean(mean(I));
figure, imshow(I2,[0,5000]);

%% apply scatter removal mask and output images
%Create a result saving folder
cd(folderPath);
outputFolder2='scatterRemoval';
if ~exist(outputFolder2, 'dir')
    mkdir(outputFolder2);
end
cd(outputFolder2);
cd('scatterRemoval');
for ind_Img=1:totImg
    Img = imread(filename{ind_Img});
    Img(BW)=mean(mean(Img));
    [~, imageName, imageExt] = fileparts(filename{ind_Img});
    imwrite(Img, [imageName, imageExt]);
end
cd(folderPath);

%% Light source correction
clear
close all
clc
% Read image folder
folderPath = 'C:\Users\yujingsong\OneDrive\Desktop\230929_S0486_R1387_onepot\6_bk_corr_aligned\scatterRemoval'; %Project folder name

%file name detection
filetype='*.tif'; 
f=fullfile(folderPath,filetype);
d=dir(f);
filename=fullfile(folderPath,{d(:).name});
totImg=length(filename);

ImgSetNum=4;
outputFolder3 = 'LScorr_images_v3';  % Folder to save aligned images;

% Create a folder to store the aligned images
cd(folderPath);
if ~exist(outputFolder3, 'dir')
    mkdir(outputFolder3);
end
cd(outputFolder3);

for i=1:ImgSetNum
    Img_ls(:,:,i) = double(imread(filename{i})).*0.5; %default is 0.8
end

for i=1:(totImg-ImgSetNum)
    Img=imread(filename{i+ImgSetNum});
    if mod(i,ImgSetNum)==0
        %Img_corr=double(Img)-Img_ls(:,:,ImgSetNum);
        Img_corr=double(Img);
    elseif mod(i,ImgSetNum)==1
        Img_corr=double(Img);
    else
        Img_corr=(double(Img)-Img_ls(:,:,mod(i,ImgSetNum))).*2;
    end
    Img_corr(Img_corr < 0) = 0;
    Img_corr=uint16(Img_corr);
    [~, imageName, imageExt] = fileparts(filename{i});
    imwrite(Img_corr, [imageName, imageExt]);
end
cd(folderPath);

