%Created by Yujing Song 231007, 454 auto-ROI
%Last updated...

%% User configuration, Open ImageJ
clear
close all
clc

javaaddpath 'C:\Users\yujingsong\OneDrive\Documents\MATLAB\454 ROI\ij.jar';
javaaddpath 'C:\Users\yujingsong\OneDrive\Documents\MATLAB\454 ROI\mij.jar';
MIJ.start;

% Read image folder
folder='C:\Users\yujingsong\OneDrive\Desktop\230929_S0486_R1387_onepot\6_bk_corr_aligned\scatterRemoval';     %Input your folder's path here!
uniqueName="000";    %Qred unique name Should be either 'Qred', 'QRed', 'Q'
filetype='*.tif';     %Make sure you save in JPG format!
outputpath='C:\Users\yujingsong\OneDrive\Desktop\230929_S0486_R1387_onepot\6_bk_corr_aligned\scatterRemoval';    %path to save your results

CutOff=1000;  %modify the threshold here   default 200-500 
%sizeX=3400;   %Crop size, X is col  3400
%sizeY=2930;   %Y is row   2800
%crop=[350 30 (sizeX-1) (sizeY-1)];  %modify this crop parameter based on array design

f=fullfile(folder,filetype);
d=dir(f);
filename=fullfile(folder,{d(:).name});
filename=filename(contains(filename, uniqueName));
totImg=length(filename);
seq = strel("disk",2);
disp('Init completed!');

%% Read a image and check crop
Img_path=filename{5}
I = imread(Img_path);
%I = imcrop(I,crop);
figure, imshow(I,[0,5000]);


DotRec_result=struct;
ij.IJ.run("Clear Results");
ij.IJ.open(Img_path);
%ij.IJ.setTool("rectangle");
%ij.IJ.makeRectangle(crop(1), crop(2), crop(3), crop(4));
%ij.IJ.run("Median...", "radius=0.1");   %remove Salt+pepper noise
ij.IJ.run("Find Maxima...","prominence="+num2str(CutOff)+" output=[Point Selection]");
ij.IJ.run("Set Measurements...","area mean min redirect=None decimal=3");
ij.IJ.run("Measure");
Measurements=MIJ.getResultsTable();
if length(Measurements)==4
    DotRec_result.Int=0;
    DotRec_result.DotsNumber=0;
    DotRec_result.Int_Mean=0;
    DotRec_result.Int_Std=0;
    DotRec_result.Int_CV=0;
    DotRec_result.DotsPosition=[];
else
    DotRec_result.Int=Measurements(:,2);
    DotRec_result.DotsNumber=length(DotRec_result.Int);
    DotRec_result.Int_Mean=mean(Measurements(:,2));
    DotRec_result.Int_Std=std(Measurements(:,2));
    DotRec_result.Int_CV=DotRec_result.Int_Std/DotRec_result.Int_Mean;
    %DotRec_result.DotsPosition=Measurements(:,5:6)-[crop(1)-1, crop(2)-1];
    DotRec_result.DotsPosition=Measurements(:,5:6);
end
ij.IJ.run("Clear Results");
MIJ.run("Close All"," ");
MIJ.closeAllWindows;
MIJ.exit;

[~, imageName, imageExt] = fileparts(Img_path);
DotRec_result.imageName=imageName;
BW=imbinarize(I);
BW(:,:)=0;

if ~isempty(DotRec_result.DotsPosition)
    Cen_W=DotRec_result.DotsPosition;
    for i=1:DotRec_result.DotsNumber
        BW(Cen_W(i,2),Cen_W(i,1))=1;
    end
    BW2 = imdilate(BW, seq);
    figure, imshow(BW2);

    disp([imageName ' Total DotsNumber=' num2str(DotRec_result.DotsNumber)]);
else
    DotRec_result.DotsNumber_corr=0;
    disp('Warning: No well has been detected in this image!');
end
planeRGB = uint8(0);
planeRGB = cat(3, (uint8(I./256).*10), (BW2.*255), BW);
figure, imshow(planeRGB,[]);
%% Shut down imageJ
%ij.IJ.run("Clear Results");
%MIJ.run("Close All"," ");
%MIJ.closeAllWindows;
%MIJ.exit;