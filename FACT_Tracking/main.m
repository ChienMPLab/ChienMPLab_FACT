%% clear and enable parloops ...
clear; close all; enableParloop(); clc; ROI = 0; cd (['C:\2D_TGMM\UFO2\Realtime_MultiROIs_fast2\ROI' num2str(ROI)]); region = sprintf('_ROI%01d', ROI); 

%%
autoTagging = 1; target = 'FMCs'; % 'SMCs'; disp(['In this experiment of ROI-' num2str(ROI) ', your target cells are ' target]);
[folderNameUFO2, prefixGrayscaleImage, imgWidth, imgHeight, numFrames, startFrame, endFrame, time_per_frame, ...
    gapClosingDist, gapClosingWindow, minimum_frames_presented] = setup_paths();

myRT = realtime(); myParser = xmlParser();

[path, folderFakeGMs, folderRawImages, folderXML, folderMamut, folderTB, grayscaleFormat, folderGrayscaleImages, ...
                folderForPhototag, imgPath, folderForAutotag] = myRT.initializer(folderNameUFO2, prefixGrayscaleImage, numFrames, ROI);    

myRT.refreshFolders(folderTB, folderMamut, folderFakeGMs, 1);
     
%%
t_start = tic; numCorrected = 0; numDivision = 0; frame = startFrame; 
GMMs = cell(1, numFrames);

while(frame <= endFrame)    
    tic; cd (path);    
    gsFile = [folderGrayscaleImages '\' prefixGrayscaleImage sprintf('%05d_ROI%01d.tif', frame, ROI)];    
    
    FLAG = 0;
    while (FLAG == 0)
        t = toc;
        if isfile(gsFile)
            disp(['Segmentation of frame ' num2str(frame) ' is detected ... start tracking ...']);
            ImgGS = imread(gsFile);
            FLAG = 1;
        else
            disp(['Segmentation of frame ' num2str(frame) ' is not generated ... waiting ...time elapsed ... ' num2str(floor(t)) ' sec']);
        end
        pause(1);
    end
    
    if FLAG == 1
        pause(2);
        fakeGMFile = [folderFakeGMs '\' prefixGrayscaleImage sprintf('%05d_ROI%01d.tif', frame, ROI)];
        if ~isfile(fakeGMFile)
            myRT. getFakeGM3(ImgGS, folderFakeGMs, folderRawImages, prefixGrayscaleImage, frame, ROI); pause(1);     
        end
        Img = imread(fakeGMFile);    
        [Area,Centroid, MeanIntensity, MaxIntensity] = myRT.getRegionprops(ImgGS, Img);
        
        for i = 1:size(Area, 1)
            m = [Centroid(i, :) 0]; id = i-1; lineage = id; parent = -1;
            area = Area(i); meanIns = MeanIntensity(i); maxIns = MaxIntensity(i);
            gmm = GMM(); gmm = gmm.initializer(m, id, lineage, parent,area, meanIns, maxIns);
            GMMs{1, frame+1}(id+1) = gmm;
        end
            
        if frame == startFrame
            myRT.write2xml(GMMs{1, frame+1}, folderXML, sprintf('GMEMfinalResult_frame%04d.xml', frame));
        else
            % to track
            gmm_1 = GMMs{1, frame};
            gmm_2 = GMMs{1, frame+1};            
            coordinates1 = reshape([gmm_1.m], [3, size(gmm_1, 2)]); coordinates1 = coordinates1(1:2, :)';
            coordinates2 = reshape([gmm_2.m], [3, size(gmm_2, 2)]); coordinates2 = coordinates2(1:2, :)';
            dist = pdist2(coordinates2, coordinates1);      
            
            orphanCell_idx = [];            
            for d = 1:size(dist, 1)
                tmp = dist(d, :);
                tmp_idx = find(tmp <= gapClosingDist);
                if isempty(tmp_idx) % assign new cell
                    orphanCell_idx = [orphanCell_idx; d];
                else
                    tmp_dist = tmp(tmp_idx);
                    if length(tmp_idx) == 1
                        gmm_2(d).parent = gmm_1(tmp_idx).id;
                    else
                        minIns2choose = [];
                        maxIns2choose = [];
                        areas2choose = [];
                        for t = 1:length(tmp_idx)
                            minIns2choose = [minIns2choose; gmm_1(tmp_idx(t)).meanIns];
                            maxIns2choose = [maxIns2choose; gmm_1(tmp_idx(t)).maxIns];
                            areas2choose = [areas2choose; gmm_1(tmp_idx(t)).area];
                        end                        
                        myMeanIns = gmm_2(d).meanIns; myMaxIns = gmm_2(d).maxIns; myArea = gmm_2(d).area;
                        
                        delta_dist = tmp_dist'/gapClosingDist; 
                        delta_meanIns = double(abs(minIns2choose)/myMeanIns); 
                        delta_maxIns = double(abs(maxIns2choose)/myMaxIns); 
                        delta_area = abs(areas2choose)/myArea;
                        
                        cost_weights = [0.7 0.10 0.10 0.10]'; % dist, meanIns, maxIns, area
                        cost_matrix = [delta_dist delta_meanIns delta_maxIns delta_area];                        
                        cost = cost_matrix*cost_weights;
                        
                        bestId = find(cost == min(cost)); bestId = bestId(1);
                        gmm_2(d).parent = gmm_1(tmp_idx(bestId)).id;
                    end
                end
            end
            
            GMMs{1, frame+1} = gmm_2;
            myRT.write2xml(GMMs{1, frame+1}, folderXML, sprintf('GMEMfinalResult_frame%04d.xml', frame));        
            
        end
        
    end
    delta_t = toc; 
    
    frame = frame + 1; disp(['Took ' num2str(delta_t) ' seconds.']);
end

%% Starting here are post analysis 
pause(2); cd (folderXML);

GMMs_mamut2 = myParser.readXMLs('GMEMfinalResult_frame', startFrame+1, endFrame+1); % If this crash, switch to .readXMLs2
myCellTracks = cellTracks;
myCellTracks.fullTracks = myCellTracks.getFullCellTracks(startFrame+1, endFrame+1, GMMs_mamut2); fTracks = myCellTracks.fullTracks;
myCellTracks.partialTracks = myCellTracks.getPartialCellTracks(startFrame+1, endFrame+1); pTracks = myCellTracks.partialTracks;
myCellTracks.deadTracks = myCellTracks.identifyDeadCellTracks(startFrame+1, endFrame+1, GMMs_mamut2); dTracks = myCellTracks.deadTracks;
tic
rescueDist = gapClosingDist*2; gapClosingWindow = 2; [myCellTracks, numRescued] = myCellTracks.rescue2(gapClosingWindow, rescueDist, 1); % slidingWindow, distance, maxNumLoop
toc

[idxMatrix, coordMatrix2X, coordMatrix2Y, SVidxMatrix] = myCellTracks.getMatrix(startFrame+1, endFrame+1);
save([path '\tables\coordMatrix2X.mat'], 'coordMatrix2X'); save([path '\tables\coordMatrix2Y.mat'], 'coordMatrix2Y'); save([path '\tables\idxMatrix.mat'], 'idxMatrix');
csvwrite([path '\tables\coordMatrix2X.csv'], coordMatrix2X); csvwrite([path '\tables\coordMatrix2Y.csv'], coordMatrix2Y); 
csvwrite([path '\tables\idxMatrix.csv'], idxMatrix);
disp('Finished retrieving coordinate tables, MOVE TO NEXT SECTION!')

% retrieving superpixels folder
cd (path);

close all; clear matrixDistanceSpeed;
matrixDistanceSpeed = myRT.getDistanceSpeedTable(coordMatrix2X, coordMatrix2Y, minimum_frames_presented, time_per_frame);
save([folderTB '\matrixDistanceSpeed.mat'],'matrixDistanceSpeed');
table_matrixDistanceSpeed = array2table(matrixDistanceSpeed,'VariableNames',{'CellIndex', 'X','Y','No_Frames', 'AccumulatedDist', 'SpeedByAccuDist', 'NetDist', 'SpeedByNetDist', 'DispRatio', 'maxiJump'});
writetable(table_matrixDistanceSpeed, [folderTB '\tb_matrixDistanceSpeed_ROI' num2str(ROI) '.xlsx']);
disp('Calculation finishes and the table is saved, MOVE TO NEXT SECTION!');

% t_end = toc(t_start); disp(['Realtime tracking on ROI-' num2str(ROI) ' used ' num2str(t_end) ' seconds.']);

validId = [];
for i = 1:size(coordMatrix2X, 1)
    tmp = find(coordMatrix2X(i, :) == 0);
    if isempty(tmp)
        validId = [validId; i];
    end
end

coordMatrix2X_valid = coordMatrix2X(validId, :);
coordMatrix2Y_valid = coordMatrix2Y(validId, :);
csvwrite([folderTB '\coordMatrix2X_valid_ROI' num2str(ROI) '.csv'], coordMatrix2X_valid); 
csvwrite([folderTB '\coordMatrix2Y_valid_ROI' num2str(ROI) '.csv'], coordMatrix2Y_valid);

disp([num2str(size(coordMatrix2X_valid, 1)) '/' num2str(size(coordMatrix2X, 1)) ' cells are valid.']);

%% Sorting the table by ... AND GET SPEED DISTRIBUTION PLOT
sort_by_accumulatedSpeed = sortrows(matrixDistanceSpeed,6,{'descend'}); sort_by_netSpeed = sortrows(matrixDistanceSpeed,8,{'descend'});
sort_by_accumulatedDist = sortrows(matrixDistanceSpeed,5,{'descend'}); sort_by_netDist = sortrows(matrixDistanceSpeed,7,{'descend'});

% Plot the histgram to determine 
close all; [figure1, cutoff, data1, N, ss, distri] = myRT.plotHist([matrixDistanceSpeed], 8); 
col = 8; 
switch col % 6: accumulatedSpeed; 8: netSpeed; 5: accumuDist; 7: netDist;
    case 5
        sort_by_what = sort_by_accumulatedDist;
    case 6
        sort_by_what = sort_by_accumulatedSpeed;
    case 7
        sort_by_what = sort_by_netDist;
    otherwise
        sort_by_what = sort_by_netSpeed;
end
sort_by_what = sort_by_what(1:ss, :); 

selectionUsing = 'ratio'; cutoff = 1; target = target;
if isequal(target, 'SMCs')
    sort_by_what =  sortrows(sort_by_what, col, {'ascend'});
    selectionUsing = 'ratio'; 
end

if isequal(selectionUsing, 'cutoff') % In this case we use ratio, e.g., 5% to determine
    speedCutoff = cutoff;    
else 
    r = 0.01;
    speedCutoff = sort_by_what(round(ss*r), col);
end

if isequal(target, 'SMCs')
    % write down here the SMCs selection
    fastOnes = intersect(find(sort_by_what(:, col) < speedCutoff), find(sort_by_what(:, col) >0));  % *(endFrame - 1)*time_per_frame
    str1 = ['Slow cell count: ' num2str(length(fastOnes))]; 
else    
    fastOnes = find(sort_by_what(:, col) > speedCutoff);  % *(endFrame - 1)*time_per_frame
    str1 = ['Fast cell count: ' num2str(length(fastOnes))]; 
end
fastRatio = length(fastOnes)/ss;
str2 = ['Total cell count: ' num2str(ss)];
str3 = ['Ratio: ' num2str(round(fastRatio*100)) '%']; str4 = ['Chosen cutoff: ' num2str(speedCutoff)];
xline(speedCutoff, 'k--'); hold on;
text(cutoff, 0.8*N, str1); text(cutoff, 0.7*N, str2); text(cutoff, 0.6*N, str3); text(cutoff, 0.5*N, str4); hold off; disp([str1 '   ' str2 '   ' str3 '   ' str4]);
saveas(figure1, [folderTB '\speed_ROI' num2str(ROI) '.png']); saveas(figure1, [folderTB '\speed_ROI' num2str(ROI) '.fig']);
% t_end = toc(t_start); disp(['Realtime tracking on ROI-' num2str(ROI) ' used ' num2str(t_end) ' seconds.']);

%% Get coordinates of interest and send to Galvo
% sort_by_what = sort_by_netDist; 
Idx = sort_by_what(1:length(fastOnes), 1); coordinate = sort_by_what(1:length(fastOnes), 2:3); disp(['You found ' num2str(length(Idx)) ' cells of interest.']);

write2csv4imageJ(sort_by_what(1:length(fastOnes), 2:3), [folderTB '\cellsOfInterest']) ;
savename2 = ['cells4tagging' region sprintf('_CS%05d', length(Idx))]; copyfile([folderTB '\cellsOfInterest.csv'], [folderForPhototag '\ImageJ_' savename2 '.csv']);

%% Send for tagging
if autoTagging == 1
    tagTime = 2000; coordinate(:,3) = tagTime; % time = 2000 ms
    csvwrite([folderForAutotag '\Cells4Tagging' region '_galvo_' sprintf('CS%05d', length(Idx)) '.csv'], coordinate);
end

%% till here you are done with UFO2
% Visualization of selection
img = imread(imgPath); 
figure; imshow(img); hold on;
Idx = sort_by_what(1:length(fastOnes), 1);
X = coordMatrix2X; Y = coordMatrix2Y;
for k = 1:length(Idx)
    i = Idx(k);
    disp([num2str(i) '/' num2str(size(X, 1))]);    
    for j = 1:size(X, 2) - 1
        x0 = X(i, j);
        y0 = Y(i, j);        
        x = X(i, j+1);
        y = Y(i, j+1);        
        if (x0*y0*x*y == 0)
            continue;
        else
            plot([x0 x], [y0 y], 'y','LineWidth', 1); hold on;
        end
    end
    x = X(i, end);
    y = Y(i, end);
    plot(x, y, '.m', 'MarkerSize', 5); hold on;     
end
saveas(gcf, [folderTB '\' 'trajectory_ROI' num2str(ROI) '.png']);

%% Saving log
tagTime = 2000;
fileID = myRT.getRunInfo(startFrame, endFrame, time_per_frame, ...
    minimum_frames_presented, speedCutoff, length(fastOnes), fastRatio, tagTime, gapClosingWindow, gapClosingDist, idxMatrix, ROI);
