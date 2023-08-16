classdef realtime
    
    properties
        maxiJumpThr = 50;        
    end
    
    methods(Static)
        
        function refreshFolders(folderTB, folderMamut, folderFakeGMs, toEmpty)
            if toEmpty == 1
                delete([folderTB '\*']);
                delete([folderMamut '\*']);
                delete([folderFakeGMs, '\*']);
            end
        end
        
        function [path, folderFakeGMs, folderRawImages, folderXML, folderMamut, folderTB, grayscaleFormat, folderGrayscaleImages, ...
                folderForPhototag, imgPath, folderForAutotag] = initializer(folderNameUFO2, prefixGrayscaleImage, numFrames, ROI)
            path = ['C:\2D_TGMM\UFO2\Realtime_MultiROIs_fast2\ROI' num2str(ROI)];  % path of the whole bunch thing

            folderXML = [path '\xmlOutput']; % folder saves all results from tracking
            if ~exist(folderXML, 'dir')
                mkdir(folderXML);
            end
            folderMamut = [folderXML '\xmlMamut'];
            if ~exist(folderMamut, 'dir')
                mkdir(folderMamut);
            end
            
            folderTB = [path '\tables'];
            if ~exist(folderTB, 'dir')
                mkdir(folderTB);
            end
            
            folderFakeGMs = [path '\fakeGMs'];
            if ~exist(folderFakeGMs, 'dir')
                mkdir(folderFakeGMs);
            end
            
            grayscaleFormat = '.tif'; % YOU MAY CHANGE HERE
            folderGrayscaleImages = [folderNameUFO2 '\ImgGS']; % folder saves all grayscale images (segmentation from Jason)
            folderForPhototag = [folderNameUFO2 '\results']; % folder saves results used for photo-tagging
            if ~exist(folderForPhototag, 'dir')
                mkdir(folderForPhototag);
            end
            folderRawImages = [folderNameUFO2 '\637nm'];
            folderForAutotag = 'D:\UFO2\Galvo';
            imgPath = [folderNameUFO2 '\637nm\' prefixGrayscaleImage sprintf('%05d_ROI%01d.tif', numFrames-1, ROI)]; % used for checking tracking trajectories
        end
        
               
        function [GMM_Rwl, GMM_Res, XML, GMMs_mamut] = setTmpMemories(endFrame, rescueWindowSize)
            % preset memories: GMM_Rwl, Res, XML are for tmp use
            GMM_Rwl = cell(2, 2); % save GMMs after removing wrong links
            GMM_Res = cell(1, 2*rescueWindowSize+1);
            XML = cell(1, 2*rescueWindowSize+1);
            % save final results per frame
            GMMs_mamut = cell(1, endFrame);
        end
        
        function pixels = getPixelsOfDisk(cx, cy, radius, width, height)
            % get range
            x0 = max(1, cx-radius);
            x1 = min(cx+radius, width);
            y0 = max(1, cy-radius);
            y1 = min(cy+radius, height);
            
            [x, y] = meshgrid(x0:x1, y0:y1);
            pixels = [];
            
            dist = sqrt((x-cx).^2 + (y-cy).^2);
            
            for i = 1:size(x, 1)
                for j = 1:size(x, 2)
                    if dist(i, j) <= radius
                        pixels = [pixels; [x(i,j)-1 y(i, j)-1]];
                    end
                end
            end            
        end

        function getFakeGM3(ImgGS, folderFakeGMs, folderRawImages, prefixGrayscaleImage, frame, ROI)
            ImgEdge = edge(ImgGS); %imgEdge = imdilate(imgEdge, strel('disk', 1));
            ImgEdge = ImgEdge > 0; 
            ImgGS2 = imsubtract(ImgGS > 0, ImgEdge > 0);
            ImgRaw = imread([folderRawImages '\' prefixGrayscaleImage sprintf('%05d_ROI%01d.tif', frame, ROI)]);
            ImgGM = uint16(ImgGS2).*ImgRaw;
            imwrite(ImgGM, [folderFakeGMs '\' prefixGrayscaleImage sprintf('%05d_ROI%01d.tif', frame, ROI)]);
        end
        
        function getFakeGM2(centroidsFile, width, height, frame, cropWidth, cropHeight, cropX0, cropY0, folderFakeGMs,folderRawImages, prefixGrayscaleImage, ROI)
            tb = readmatrix(centroidsFile);
            y = floor((tb-1)/width);
            x = tb-1-y*width;
            
            centroids_tmp = [x y];
            centroids = centroids_tmp;
            
            if cropWidth*cropHeight > 0
                tmp1 = centroids_tmp - [cropX0 cropY0];
                tmp2 = [cropX0+cropWidth  cropY0+cropHeight] - centroids_tmp;
                tmp = [tmp1 tmp2];
                for i = 1:size(tmp, 2)
                    tmp(:, i) = max(tmp(:, i), 0);
                end
                result = prod(tmp, 2);
                ids = find(result > 0);
                centroids = centroids_tmp(ids, :);
            end
            
            x = centroids(:, 1); y = centroids(:, 2);
            
            imgGM = zeros(height, width, 'uint16');
            
            for c = 1:size(centroids, 1)
                y(c) = max(1, y(c));
                x(c) = max(1, x(c));
                y(c) = min(height, y(c));
                x(c) = min(width, x(c));                
                imgGM(y(c), x(c)) = 1;
            end
            SE = strel("disk", 4);
            imgGM = imdilate(imgGM, SE); 
            imgGM = imgGM.*(imread([folderRawImages '\' prefixGrayscaleImage sprintf('%05d_ROI%01d.tif', frame, ROI)]));
            imwrite(imgGM, [folderFakeGMs '\' prefixGrayscaleImage sprintf('%05d_ROI%01d.tif', frame, ROI)]);
        end
        
        function getFakeGM(mycsv, width, height, frame, cropWidth, cropHeight, cropX0, cropY0, folderFakeGMs, folderRawImages, prefixGrayscaleImage, ROI)
            tb = readmatrix(mycsv);
            y = floor((tb-1)/width);
            x = tb-1-y*width;
            
            centroids_tmp = [x y];
            centroids = centroids_tmp;
            
            if cropWidth*cropHeight > 0
                tmp1 = centroids_tmp - [cropX0 cropY0];
                tmp2 = [cropX0+cropWidth  cropY0+cropHeight] - centroids_tmp;
                tmp = [tmp1 tmp2];
                for i = 1:size(tmp, 2)
                    tmp(:, i) = max(tmp(:, i), 0);
                end
                result = prod(tmp, 2);
                ids = find(result > 0);
                centroids = centroids_tmp(ids, :);
            end
            
            x = centroids(:, 1); y = centroids(:, 2);
            
            % taking 0.3 sec
            parfor i = 1:size(centroids, 1)
                cx = x(i);     cy = y(i);
                pixels = getPixelsOfDisk(cx, cy, 5, width, height);
                val = width*pixels(:, 2) + pixels(:, 1) + 1;
                pixelList{i, 1} = val;
            end
            
            imgGM = zeros(height, width, 'uint16');
            
            for i = 1:size(pixelList, 1)
                pixels = pixelList{i, :}
                aa = find(pixels > 0);
                pixels = pixels(1:max(aa));
                %     intensity = 200 + floor(rand*55);
                intensity = 1;
                for j = 1:size(pixels, 2)
                    y = floor((pixels(j)-1)/width);
                    x = pixels(j)-1-width*y;
                    y = y + 1; x = x + 1;
                    %         imgGS(y, x) = i + 100;
                    imgGM(y, x) = intensity;
                end
            end
            imgGM = imgaussfilt(imgGM, 1);
            imgGM = imgGM.*(imread([folderRawImages '\' prefixGrayscaleImage sprintf('%05d_ROI%01d.tif', frame, ROI)]));            
           
            imwrite(imgGM, [folderFakeGMs '\' sprintf('file%05d.tif', frame)]);
            
        end
        
        function [Area,Centroid, MeanIntensity, MaxIntensity] = getRegionprops(ImgGS, Img)
            tb_GM = regionprops('table', Img > 0, Img, 'Centroid', 'Area', 'MeanIntensity', 'MaxIntensity');           
            tb_GS = regionprops('table', Img > 0, ImgGS, 'Centroid', 'Area', 'MeanIntensity', 'MaxIntensity', 'PixelList');
            a = tb_GS.MeanIntensity - double(tb_GS.MaxIntensity);      
            b = find(tb_GM.Area <= 50);            
            idx2remove = sort([find(a ~= 0); b], 'descend');            
            pixels = [];            
            for b = 1:length(idx2remove)
                id = idx2remove(b);
                tb_GM(id, :) = [];
                x = round(tb_GS.Centroid(id, 1)); y = round(tb_GS.Centroid(id, 2));
                pl = tb_GS.PixelList{id};
                pixels = [pixels; pl];
            end             
            array_GM = table2array(tb_GM);
            if ~isempty(pixels)
                tmp = zeros(size(Img));
                for p = 1:size(pixels, 1)
                    tmp(pixels(p, 2), pixels(p, 1)) = 1;
                end
                tmp = uint16(tmp).*ImgGS;
                cellIdx = unique(tmp);                
                toAdd = [];
                for b = 2:length(cellIdx)
                    id = cellIdx(b);
                    tmp2 = (tmp == id);
                    tmp_tb = regionprops('table', tmp2, Img, 'Centroid', 'Area', 'MeanIntensity', 'MaxIntensity');
                    toAdd = [toAdd; tmp_tb.Area tmp_tb.Centroid tmp_tb.MeanIntensity tmp_tb.MaxIntensity];
                end
                array_GM = [array_GM; toAdd];
            end
            array_GM = double(array_GM);
            Area = array_GM(:, 1); Centroid = array_GM(:, 2:3); MeanIntensity = array_GM(:, 4); MaxIntensity = array_GM(:, 5);

        end
        
        function write2xml(gmm, path2save, xmlFilename)
            doc = com.mathworks.xml.XMLUtils.createDocument('document');
            document = doc.getDocumentElement;
            
            for ww = 1:size(gmm, 2)
                myFile = doc.createElement('GaussianMixtureModel');
                myFile.setAttribute('id', num2str(gmm(ww).id));
                myFile.setAttribute('lineage',num2str(gmm(ww).lineage));
                myFile.setAttribute('parent',num2str(gmm(ww).parent));
                myFile.setAttribute('dims',num2str(gmm(ww).dims));
%                 myFile.setAttribute('area',num2str(gmm(ww).area));
%                 myFile.setAttribute('minIns',num2str(gmm(ww).meanIns));
%                 myFile.setAttribute('maxIns',num2str(gmm(ww).maxIns));
                myFile.setAttribute('splitScore',num2str(gmm(ww).splitScore));
                myFile.setAttribute('svIdx',num2str(gmm(ww).svIdx));
                str = '';
                for i = 1:size(gmm(ww).m, 2)
                    str = [str num2str(gmm(ww).m(i)) ' '];
                end
                myFile.setAttribute('m', str);               
                
                str = '';
                for i = 1:size(gmm(ww).scale, 2)
                    str = [str num2str(gmm(ww).scale(i)) ' '];
                end
                myFile.setAttribute('scale',str);
                myFile.setAttribute('nu',num2str(gmm(ww).nu));
                myFile.setAttribute('beta',num2str(gmm(ww).beta));
                myFile.setAttribute('alpha',num2str(gmm(ww).alpha));                
                str = '';
                for i = 1:size(gmm(ww).W, 2)
                    str = [str num2str(gmm(ww).W(i)) ' '];
                end                
                myFile.setAttribute('W',str);
                str = '';
                for i = 1:size(gmm(ww).WPrior, 2)
                    str = [str num2str(gmm(ww).WPrior(i)) ' '];
                end
                myFile.setAttribute('WPrior',str);
                
                myFile.setAttribute('nuPrior',num2str(gmm(ww).nuPrior));
                myFile.setAttribute('betaPrior',num2str(gmm(ww).betaPrior));
                myFile.setAttribute('alphaPrior',num2str(gmm(ww).alphaPrior));
                myFile.setAttribute('distMRFPrior',num2str(gmm(ww).distMRFPrior));
                str = '';
                for i = 1:size(gmm(ww).mPrior, 2)
                    str = [str num2str(gmm(ww).mPrior(i)) ' '];
                end
                myFile.setAttribute('mPrior',str);
                
                
                document.appendChild(myFile);
            end
            
            xmlwrite([path2save '\' xmlFilename], doc);
            clear doc;
        end

        
        
       
    end
    
    methods            
        
        function writingMask2CSV(obj, filename, folderFakeGMs, folderPixelList, frame)
            tic
            disp('writing fakeGM and pixellist ...');
            img = imread(filename); level = graythresh(img);
            bw = im2bw(img, 0.1*level); bw = bw';
            I = imgaussfilt(im2double(bw), 1);
            fakeGM = I;
            imwrite(fakeGM',[folderFakeGMs '\' sprintf('file0%04d',frame) '.tif']);            
        end
        
        function matrixDistanceSpeed = getDistanceSpeedTable(obj, coordMatrix2X, coordMatrix2Y, ...
                minimum_frames_presented, time_per_frame)
            totalCellCount = size(coordMatrix2X, 1);
            matrixDistanceSpeed = zeros(totalCellCount, 10);
            %             T = size(coordMatrix2X, 2);
            
            for i = 1:totalCellCount
                xs = coordMatrix2X(i, :);
                ys = coordMatrix2Y(i, :);       
                xpos = find(xs > 0); 
                ypos = find(ys > 0);
                pos = intersect(xpos, ypos);
                
                if isempty(pos)
                    continue;
                end
                
                if length(xpos) > length(ypos)
                    pos = ypos;
                end
                xs = xs(pos);
                ys = ys(pos);

                matrixDistanceSpeed(i, 1) = i; % the index
                matrixDistanceSpeed(i, 2) = xs(end); % x
                matrixDistanceSpeed(i, 3) = ys(end); % y
                noFrames = size(xs, 2);
                matrixDistanceSpeed(i, 4) = noFrames; % number of frames the given cell persists
                
                xy = [xs; ys];
                accumulatedDistance = 0;
                accumulatedSpeed = 0;
                netDistance = 0;
                netSpeed = 0;
                displacementRatio = 0;
                
                maxiJump = 0;
                
                if noFrames >= minimum_frames_presented
                    
                    for j = 1:noFrames - 1
                        a = (xy(:, j+1) - xy(:, j)).^2;
                        dist = sqrt(sum(a));
                        if dist > maxiJump
                            maxiJump = dist;
                        end
                        accumulatedDistance = accumulatedDistance + dist;
                    end
                    
                    accumulatedSpeed = accumulatedDistance /(time_per_frame * (noFrames - 1));
                    netDistance = sqrt(sum((xy(:, 1) - xy(:, noFrames)).^2));
                    netSpeed = netDistance / (time_per_frame * (noFrames - 1));
                    displacementRatio = maxiJump / netDistance;
                end
                
                if maxiJump > obj.maxiJumpThr
                    accumulatedDistance = 0;
                    accumulatedSpeed = 0;
                    netDistance = 0;
                    netSpeed = 0;
                    displacementRatio = 0;
                end
                matrixDistanceSpeed(i, 5) = accumulatedDistance;
                matrixDistanceSpeed(i, 6) = accumulatedSpeed;
                matrixDistanceSpeed(i, 7) = netDistance;
                matrixDistanceSpeed(i, 8) = netSpeed;
                matrixDistanceSpeed(i, 9) = displacementRatio;
                matrixDistanceSpeed(i, 10) = maxiJump;                
            end
        end
        
     
        function [f, cutoff, data, N, s, distri] = plotHist(obj, matrixDistanceSpeed, col)
            data = matrixDistanceSpeed(:, col);
            data1 = [];
            for d = 1:length(data)
                if data(d) > 0
                    data1 = [data1; d data(d)];
                end
            end
            data1 = data1(:, 2); data1 = sort(data1, 'descend');
            s = size(data1, 1);
            
            ratio = 0.05; num = round(s*ratio);            
            data2 = data1(1:num);
            
            
            numClusters = 3; 
            
            %Clustering using GMM
            distri = fitgmdist(data2, numClusters);
            c = sort(distri.mu); id = find(distri.mu == c(1)); 
            
            cutoff = distri.mu(id)+1*distri.Sigma(id);
            
%             % plotHist(data2);
%             distri = fitgmdist(data2, 3);
%             c = sort(distri.mu); id = find(distri.mu == c(2));
%                        
%             cutoff = distri.mu(id)-3*distri.Sigma(id);
            cutoff = round(cutoff*100)/100;
            
            f = figure;
            H = histogram(data1, ceil(sqrt(s)+10), 'FaceColor', [0.4660 0.6740 0.1880], 'FaceAlpha', 0.5); hold on;
            N = max(histcounts(data1, ceil(sqrt(s)+10)));
            axis([-0.05, 2, 0 N]); box off;  xlabel('Pixels per minute'); ylabel('Counts');

            title(['Estimated speed cutoff is ' sprintf('%1.2f', cutoff)]);
%             selectedCellCount = size(find(data1 > cutoff), 1);
%             
%             text(cutoff + 0.05, 0.8*N, ['Total cell count: ' num2str(s)]);
%             text(cutoff + 0.05, 0.7*N, ['Fast cell count: ' num2str(selectedCellCount)]);
%             text(cutoff + 0.05, 0.6*N, ['Fast cell ratio: ' num2str(round(100*(selectedCellCount/s))) '%']); hold on;
%             plot([[cutoff cutoff], [0 100]], '-r', 'Linewidth', 6);
%             xline(cutoff, 'k--'); hold on;
%             plot(conv(H.BinEdges, [0.5 0.5], 'valid'), H.BinCounts, '-r', 'Color', [0.6350 0.0780 0.1840], 'Linewidth', 1.5);

           
%             f = figure; histogram(data1(:, 2), ceil(sqrt(s)+10));
%             N = max(histcounts(data1(:, 2), ceil(sqrt(s)+10)));
%             distri = fitgmdist(data1(:, 2), 2); hold on;
%             cutoff = distri.mu(2) - distri.Sigma(2);
%             xlabel('Pixels per minute'); ylabel('Counts'); 
% %             title(['Estimated speed cut-off is ' num2str(cutoff)]); 
%             axis([-0.05 2 0 N]); box off;
% %             hold off;
        end
        
        function pressEnter(Obj, event)
            import java.awt.*;
            import java.awt.event.*;
            rob = Robot;
            rob.keyPress(KeyEvent.VK_ENTER)
            rob.keyRelease(KeyEvent.VK_ENTER)
        end

        function fileID = getRunInfo(obj, startFrame, endFrame, time_per_frame, ...
                  minimum_frames_presented, speedCutoff, numFastCells, fastRatio, tagTime, gapClosingWindow, gapClosingDist, idxMatrix, ROI) 
              
                 currentTime = datestr(clock,'YYYYmmddatHHMMFFF');
             fileID = fopen(['infoTracking_ROI' num2str(ROI) '_' currentTime '.txt'],'wt');
             fprintf(fileID, 'startFrame = %6.2f\nendFrame = %6.2f\ntime_per_frame = %6.2f\n', startFrame, endFrame, time_per_frame);
%              fprintf(fileID, 'rescueWindowSize = %6.2f\nrescueDist = %6.2f\nthrWrongLink = %6.2f\n', rescueWindowSize, rescueDist, thrWrongLink);
             fprintf(fileID, 'gapClosingWindow= %6.2f\ngapClosingDist = %6.2f\n', gapClosingWindow, gapClosingDist);
             fprintf(fileID, 'minimum_frames_presented = %6.2f\n', minimum_frames_presented);
             fprintf(fileID, 'speedCutoff = %6.2f\nnum_fastCells = %6.2f\nratio_fastCells = %6.2f\n', speedCutoff, numFastCells, fastRatio);
             fprintf(fileID, 'phototag_time_per_cell = %6.2f\n', tagTime);
             n = size(idxMatrix, 1);
             k = 0;
             for i = 1:n
                 tmp = find(idxMatrix(i, :) == 0);
                 if isempty(tmp) %|| length(tmp) < endFrame*0.3
                     k = k + 1;
                 end
             end
             fprintf(fileID, 'estimated_tracking_accuracy = %6.2f\n', k/n);
             fprintf(fileID, 'cell_count = %6.2f\n', n);
             fclose(fileID);             
        end
         

        
        

    end
    
end