function [folderNameUFO2, prefixGrayscaleImage, imgWidth, imgHeight, numFrames, startFrame, endFrame, time_per_frame, ...
    gapClosingDist, gapClosingWindow, minimum_frames_presented] = setup_RTMR_fast()

folderNameUFO2 = 'D:\UFO2\_20230802atxxxxxx_HN031'; 

files = dir([folderNameUFO2 '\637nm\' '*.tif']); name = strsplit(files(1).name, '_frame'); 
prefixGrayscaleImage = [name{1} '_frame'];

UFO = 2; 

if UFO == 2
    imgHeight = 5120;
    imgWidth = 5120;  
else
    imgHeight = 3000; %3000
    imgWidth = 4096; %4096
end

numFrames = 30;

time_per_frame = 4; 

gapClosingDist = 15;

gapClosingWindow = 2;

startFrame = 0; 
endFrame = numFrames - 1;

minimum_frames_presented = 3;


disp('Pls DOUBLE-CHECK if you are running: ');
disp(['Experiment from FOLDER: ' folderNameUFO2]); disp(['Image-PREFIX: ' prefixGrayscaleImage]);
disp(['Using UFO' num2str(UFO)]);

end