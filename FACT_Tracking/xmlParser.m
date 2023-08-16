 classdef xmlParser
    
    methods
        
        function GMMs = readXMLs(obj, basename, startFrame, endFrame) % to load with MEX
            numOfFrames = endFrame - startFrame + 1;
            GMMs = cell(1, numOfFrames);
            for frame = startFrame:endFrame
                disp(['Frame ... ' num2str(frame)]);
                [gmms] = readXMLmixtureGaussians([basename num2str(frame - 1,'%.4d') '.xml']);
                GMMs{1, frame - startFrame + 1} = gmms; % 
            end
        end
        
        function GMMs = readXMLs2(obj, basename, startFrame, endFrame) % to load if MEX crashes my computer ... when dataset is tooooo large
            numOfFrames = endFrame - startFrame + 1; % if starting image is not 0 (if 0-based) or 1 (if 1-based), then we shift the image id
            GMMs = cell(1, numOfFrames);
            parfor frame = 1:numOfFrames 
                disp(['Parsing XML of frame ... ' num2str(frame)]);
                [gmms] = xml2gmm([basename num2str(startFrame + frame-2,'%.4d') '.xml']); % was with  - 2
                L = length(gmms.document.GaussianMixtureModel);
                for l = 1:L
                    gmm = gmms.document.GaussianMixtureModel(l);
                    blob = gmm{1,1}.Attributes;
                    GMMs{1, frame}(l).m = str2num(blob.m);
                    GMMs{1, frame}(l).id = str2num(blob.id);
                    GMMs{1, frame}(l).lineage = str2num(blob.lineage);
                    GMMs{1, frame}(l).parent = str2num(blob.parent);
                    GMMs{1, frame}(l).svIdx = str2num(blob.svIdx);
                    GMMs{1, frame}(l).dims = str2num(blob.dims);
                    GMMs{1, frame}(l).splitScore = str2num(blob.splitScore);
                    GMMs{1, frame}(l).scale = str2num(blob.scale);
                    GMMs{1, frame}(l).nu = str2num(blob.nu);
                    GMMs{1, frame}(l).beta = str2num(blob.beta);
                    GMMs{1, frame}(l).alpha = str2num(blob.alpha);
                    GMMs{1, frame}(l).W = str2num(blob.W);                    
                    GMMs{1, frame}(l).nuPrior = str2num(blob.nuPrior);
                    GMMs{1, frame}(l).betaPrior = str2num(blob.betaPrior);
                    GMMs{1, frame}(l).alphaPrior = str2num(blob.alphaPrior);
                    GMMs{1, frame}(l).mPrior = str2num(blob.mPrior);
                    GMMs{1, frame}(l).WPrior = str2num(blob.WPrior);
                    GMMs{1, frame}(l).distMRFPrior = str2num(blob.distMRFPrior);                    
                end
                GMMs{1, frame} = GMMs{1, frame}';
            end
        end
        
        function write2xml(obj, file2, xmlFilename)
            doc = com.mathworks.xml.XMLUtils.createDocument('document');
            document = doc.getDocumentElement;            
            for ww = 1:size(file2, 1)
                if ~isempty(file2(ww).m)
                gmm = doc.createElement('GaussianMixtureModel');
                gmm.setAttribute('id',num2str(file2(ww).id));
                gmm.setAttribute('lineage',num2str(file2(ww).lineage));
                if isempty(file2(ww).parent)
                    file2(ww).parent = -1;
                end
                gmm.setAttribute('parent',num2str(file2(ww).parent));     
                gmm.setAttribute('dims',num2str(file2(ww).dims));
                gmm.setAttribute('splitScore',num2str(file2(ww).splitScore));
                str = '';
                for i = 1:size(file2(ww).scale, 2)
                    str = [str num2str(file2(ww).scale(i)) ' '];
                end
                gmm.setAttribute('scale',str);
                gmm.setAttribute('nu',num2str(file2(ww).nu));
                gmm.setAttribute('beta',num2str(file2(ww).beta));
                gmm.setAttribute('alpha',num2str(file2(ww).alpha));
                str = '';
                for i = 1:size(file2(ww).m, 2)
                    str = [str num2str(file2(ww).m(i)) ' '];
                end
                gmm.setAttribute('m',str); 
                if isempty(file2(ww).W)                   
                    file2(ww).W = [0.0230,0,0,0,0.0230,0,0,0,0.0452];
                end
                str = '';
                for i = 1:size(file2(ww).W, 2)
                    str = [str num2str(file2(ww).W(i)) ' '];
                end                
                gmm.setAttribute('W',str);
                
                gmm.setAttribute('nuPrior',num2str(file2(ww).nuPrior));
                gmm.setAttribute('betaPrior',num2str(file2(ww).betaPrior));
                gmm.setAttribute('alphaPrior',num2str(file2(ww).alphaPrior));
                gmm.setAttribute('distMRFPrior',num2str(file2(ww).distMRFPrior));
                str = '';
                for i = 1:size(file2(ww).mPrior, 2)
                    str = [str num2str(file2(ww).mPrior(i)) ' '];
                end
                gmm.setAttribute('mPrior',str);
                if isempty(file2(ww).WPrior)
                    file2(ww).WPrior = [0.05555, 0, 0, 0, 0.05555, 0, 0, 0, 0.10888];
                end
                str = '';
                for i = 1:size(file2(ww).WPrior, 2)
                    str = [str num2str(file2(ww).WPrior(i)) ' '];
                end
                gmm.setAttribute('WPrior',str);
                str = '';
                for i = 1:size(file2(ww).svIdx, 2)
                    str = [str num2str(file2(ww).svIdx(i)) ' '];
                end
                gmm.setAttribute('svIdx',str);
                document.appendChild(gmm);
                end
            end
            xmlwrite(xmlFilename,doc);
            clear doc;
            %             type('myfile.xml');
        end
        
        function file2 = linkVertical(obj, file1, file2) % for details pls check note on April 7, 2020
%             matchedId = []; % this stores the matched id, one row is one match ... col 1 for i-th and col 2 for (i+1)-th segments
            %/* this part needs to speed up, was so slow: e.g., for 8000
            %cells taking around 30 sec 
            file1_svIdx = zeros(size(file1, 1), 1);
            for i = 1:size(file1, 1)
                svid = file1(i).svIdx;
                file1_svIdx(i, 1) = svid(1);
            end
            file2_svIdx = zeros(size(file2, 1), 1);
            file2Id = [];
            for j = 1:size(file2, 1)
                svid = file2(j).svIdx;
                file2_svIdx(j, 1) = svid(1);
                file2Id = [file2Id; file2(j).id];
            end

            inter = intersect(file1_svIdx, file2_svIdx);

            matchedId = zeros(size(inter, 2), 2);
            for m = 1:size(inter, 1)
                p1 = find(inter(m) == file1_svIdx);
                p2 = find(inter(m) == file2_svIdx);
                matchedId(m, :) = [p1 p2];
            end
%             for i = 1:size(file1, 1)
%                 match = [i 0]; % 0 means no match, id counting from 1 ...
%                 for j = 1:size(file2, 1)
%                     if isequal(file1(i).svIdx, file2(j).svIdx) % isequal([file1(i).m  file1(i).svIdx], [file2(j).m file2(j).svIdx])
%                         match = [i j]; 
%                         break;
%                     end
%                 end
%                 matchedId = [matchedId; match];
%             end
            % */
            count = size(file2, 1);
            c1 = 0; c2 = 0;
%             disp([num2str(len) ' cells are detected in file first, BUT not in file second, we will start adding cells in file second, from ' num2str(count)]);
            for k = 1:size(matchedId, 1)
                id1 = matchedId(k, 1); id2 = matchedId(k, 2);
                if id1*id2 ~= 0 % files 1&2 are identical, if cells are detected both in files 1&2, we just link them
                    c1 = c1 + 1;
                    file2(id2).parent = file1(id1).id;
                else % files 1&2 are identical, but if in file 2 we miss cells, we need to create them first ...
                    count = count + 1;
                    c2 = c2 + 1;
                    file2(count).m = file1(id1).m;
                    file2(count).id = c2 + max(file2Id); %count - 1;
                    file2(count).lineage = file2(count).id;
                    file2(count).parent = file1(id1).id;
                    file2(count).svIdx = file1(id1).svIdx;
                    file2(count).dims = file1(id1).dims;
                    file2(count).splitScore = file1(id1).splitScore;
                    file2(count).nu = file1(id1).nu;
                    file2(count).beta = file1(id1).beta;
                    file2(count).alpha = file1(id1).alpha;
                    file2(count).W = file1(id1).W;
                    file2(count).nuPrior = file1(id1).nuPrior;
                    file2(count).betaPrior = file1(id1).betaPrior;
                    file2(count).alphaPrior = file1(id1).alphaPrior;
                    file2(count).WPrior = file1(id1).WPrior;
                    file2(count).mPrior = file1(id1).mPrior;
                    file2(count).scale = file1(id1).scale;
                end
            end             
%             disp([num2str(c1) ' cells are kept, while ' num2str(c2) ' cells are added.']);
        end
        
        function file2 = linkVertical2(obj, file1, file2) % for details pls check note on April 7, 2020
%             matchedId = []; % this stores the matched id, one row is one match ... col 1 for i-th and col 2 for (i+1)-th segments
            matchedId = zeros(size(file1, 1), 2);
            file2Id = [];
            parfor i = 1:size(file1, 1) %size(file1, 1)
                matchedId(i,:) = [i 0]; % 0 means no match, id counting from 1 ...
                for j = 1:size(file2, 1)
                    file2Id = [file2Id file2(j).id];
                    if isequal(file1(i).svIdx, file2(j).svIdx) % isequal([file1(i).m  file1(i).svIdx], [file2(j).m file2(j).svIdx])
                        matchedId(i,:) = [i j];
                    end
                end
            end
            len = length(find(matchedId(:, 2) == 0));            
            
            count = size(file2, 1);
            disp([num2str(len) ' cells are detected in file first, BUT not in file second, we will start adding cells in file second, from ' num2str(count)]);

            c1 = 0; c2 = 0;
            for k = 1:size(matchedId, 1)
                id1 = matchedId(k, 1); id2 = matchedId(k, 2);
                if id1*id2 ~= 0 % files 1&2 are identical, if cells are detected both in files 1&2, we just link them
                    %                     file2(id2).parent = file1(id1).id;
                    file2(id2).parent = file1(id1).id;
                    c1 = c1 + 1;                
                else % files 1&2 are identical, but if in file 2 we miss cells, we need to create them first ...
% problem here: count is conflicting with tracking idx ... 
                    count = count + 1;
                    c2 = c2 + 1;
                    file2(count).m = file1(id1).m;
                    file2(count).id = c2 + max(file2Id);
                    file2(count).lineage = file2(count).id;
                    file2(count).parent =  file1(id1).id; % -1
                    file2(count).svIdx = file1(id1).svIdx;
                    file2(count).dims = file1(id1).dims;
                    file2(count).splitScore = file1(id1).splitScore;
                    file2(count).nu = file1(id1).nu;
                    file2(count).beta = file1(id1).beta;
                    file2(count).alpha = file1(id1).alpha;
                    file2(count).W = file1(id1).W;
                    file2(count).nuPrior = file1(id1).nuPrior;
                    file2(count).betaPrior = file1(id1).betaPrior;
                    file2(count).alphaPrior = file1(id1).alphaPrior;
                    file2(count).WPrior = file1(id1).WPrior;
                    file2(count).mPrior = file1(id1).mPrior;
                    file2(count).scale = file1(id1).scale;
                end
              end
             disp([num2str(c1) ' cells are kept, while ' num2str(c2) ' cells are added.']);
          end
        
        function file3 = linkHorizontal(obj, file2, file3)
            file2Id = [];
            for j = 1:size(file2, 1)
               file2Id = [file2Id file2(j).id];
            end
            rowNumUsed = [];
            for i = 1:size(file3, 1)
                oldPidx = file3(i).parent;
                if oldPidx >= 0
                    rowInFile2 = find(oldPidx == file2Id);
                    if ~isempty(rowInFile2)
                        rowInFile2 = rowInFile2(1);
                        rowNumUsed = [rowNumUsed rowInFile2];
                        file3(i).parent = file2(rowInFile2).parent;
                    end
                end
            end      
            rows = 1:size(file2, 1);
            rowNumUnused = setdiff(rows, rowNumUsed);
            % /* adjusted on 06/18
            xInFile2 = zeros(length(rowNumUnused), 1); yInFile2 = zeros(length(rowNumUnused), 1);
            for j = 1:length(rowNumUnused)   
                xInFile2(j, 1) = file2(rowNumUnused(j)).m(1);
                yInFile2(j, 1) = file2(rowNumUnused(j)).m(2);
            end
            % */
            for i = 1:size(file3, 1)
                if file3(i).parent > -1
                    continue;
                end
                xInFile3 = ones(length(rowNumUnused), 1) * file3(i).m(1);
                yInFile3 = ones(length(rowNumUnused), 1) * file3(i).m(2);                
                dist = sqrt((xInFile2-xInFile3).^2+(yInFile2-yInFile3).^2);                
                pos = intersect(find(dist < 10), find(dist == min(dist)));
                if isempty(pos)
                    continue;
                end
                pos = pos(1);
                file3(i).parent = file2(rowNumUnused(pos)).parent;
            end
        end
        
        function file3 = linkHorizontal3(obj, file2, file3)
            file2Id = [];
            for j = 1:size(file2, 1)
               file2Id = [file2Id file2(j).id];
            end
            rowNumUsed = [];
            for i = 1:size(file3, 1)
                oldPidx = file3(i).parent;
                if oldPidx >= 0
                    rowInFile2 = find(oldPidx == file2Id);
                    if ~isempty(rowInFile2)
                        if length(rowInFile2) > 1
                            i
                            rowInFile2
                        end
                        rowInFile2 = rowInFile2(1);
                        rowNumUsed = [rowNumUsed rowInFile2];
                        file3(i).parent = file2(rowInFile2).parent;
                    end
                end
            end      
            rows = 1:size(file2, 1);
            rowNumUnused = setdiff(rows, rowNumUsed);
%             for i = 1:size(file3, 1)
%                 for j = 1:length(rowNumUnused)
%                     if file3(i).parent == -1 
%                         xInFile3 = file3(i).m(1);
%                         yInFile3 = file3(i).m(2);
%                         xInFile2 = file2(rowNumUnused(j)).m(1); yInFile2 = file2(rowNumUnused(j)).m(2); 
%                         if ((xInFile2-xInFile3)^2+(yInFile2-yInFile3)^2) < 10^2 
%                             file3(i).parent = file2(rowNumUnused(j)).parent;
%                             rowNumUsed = [rowNumUsed rowNumUnused(j)];
%                         end
%                     end
%                 end
%             end
        end
        
        function maxId = getMaxIdPerFrame(obj, gmm)
            ids = zeros(size(gmm, 1), 1);
            for i = 1:size(gmm, 1)
                ids(i, 1) = gmm(i).id;
            end
%             ids
            maxId = max(ids);
% maxId = vertcat(gmm.id);
        end
        
        function [file3, GMMs] = linkHorizontal2(obj, GMMs, thisFrame, window, resWindowSize, minDist, file2, file3)
            % this function shall consider rescue function ... 
            file2Id = [];
            for j = 1:size(file2, 1)
               file2Id = [file2Id file2(j).id];
            end
            rowNumUsed = [];
            for i = 1:size(file3, 1)
                oldPidx = file3(i).parent;
                if oldPidx >= 0
                    rowInFile2 = find(oldPidx == file2Id);                   
                    if ~isempty(rowInFile2)
                        rowInFile2 = rowInFile2(1);
                        rowNumUsed = [rowNumUsed rowInFile2];
                        file3(i).parent = file2(rowInFile2).parent;
                    end
                end
            end            
            rows = 1:size(file2, 1);
            rowNumUnused = setdiff(rows, rowNumUsed);
            %/* adjusted June 19
            xInFile2 = zeros(length(rowNumUnused), 1); yInFile2 = zeros(length(rowNumUnused), 1);
            for j = 1:length(rowNumUnused)
                xInFile2(j, 1) = file2(rowNumUnused(j)).m(1);
                yInFile2(j, 1) = file2(rowNumUnused(j)).m(2);
            end
            %*/
            for i = 1: size(file3, 1)
%                 for j = 1:length(rowNumUnused)
                    if file3(i).parent == -1 
                        xInFile3 = ones(length(rowNumUnused), 1) * file3(i).m(1);
                        yInFile3 = ones(length(rowNumUnused), 1) * file3(i).m(2);
                        dist = sqrt((xInFile2-xInFile3).^2+(yInFile2-yInFile3).^2);
                        
                        pos = find(dist < 10);
                        
                        if ~isempty(pos)
                            file3(i).parent = file2(rowNumUnused(pos(1))).parent;
                             rowNumUsed = [rowNumUsed rowNumUnused(pos(1))];
                        end
                        
%                         xInFile2 = file2(rowNumUnused(j)).m(1); yInFile2 = file2(rowNumUnused(j)).m(2); 
%                         if ((xInFile2-xInFile3)^2+(yInFile2-yInFile3)^2) < 10^2 
% %                             disp(['Found xInFile2 ' num2str(j) '-th row']);
%                             file3(i).parent = file2(rowNumUnused(j)).parent;
% %                             /* added June 13
%                             rowNumUsed = [rowNumUsed rowNumUnused(j)];
% %                             */
%                         else 
% %                              disp(['No Found xInFile2 ']);
%                         end
                    end
%                 end
            end
            % considering the rescue part ... find out the GMM backwards
            % based on window ...
            c = 0;
            GMMs_bw = [];
            for frame = max(1, thisFrame-1-window-1):thisFrame-2 % thisFrame-1-window-1
                c = c + 1;
                GMMs_bw{1,c} = GMMs{1, frame};
            end            
%             disp(['Frames being considered are ' num2str(thisFrame-1-window-1) ' to ' num2str(thisFrame-2)]);

            for i = 1:size(file3, 1) 
                if file3(i).parent == -1
                    candidate = []; % saves frame id angle dist coordinates ..
                   
                    if isempty(GMMs_bw)
                        continue;
                    end
                    
                    for k = size(GMMs_bw, 2):-1:2 % starting frome 2, as the first one is used as reference to generate an vector for calculating angles
                        file4 = GMMs_bw{1,k};
%                         ids = []; 
                        xInFile4 = zeros(size(file4, 2), 1); yInFile4 = zeros(size(file4, 2), 1);                       
                        for j = 1:size(file4, 2) % size(file4, 1)
                            xInFile4(j, 1) = file4(j).m(1);
                            yInFile4(j, 1) = file4(j).m(2);
                        end
                        xInFile3 = ones(size(file4, 2), 1)*file3(i).m(1);
                        yInFile3 = ones(size(file4, 2), 1)*file3(i).m(2);
                        dist = sqrt((xInFile4-xInFile3).^2+(yInFile4-yInFile3).^2);
                        
                        pos = find(dist < minDist);
                        
                        for po = 1:length(pos)
                            P0 = [xInFile4(pos(po), 1) yInFile4(pos(po), 1)]; P1 = [xInFile3(1,1) yInFile3(1,1)];
                            if file4(pos(po)).parent+1 >= 1 &&  file4(pos(po)).parent+1 <=  size(GMMs_bw{1,k-1}, 2)
                                P2 = [GMMs_bw{1,k-1}(file4(pos(po)).parent+1).m(1) GMMs_bw{1,k-1}(file4(pos(po)).parent+1).m(2)];
                            else
                                P2 = P0;
                            end
                            ang = atan2(abs(det([P2-P0;P1-P0])),dot(P2-P0,P1-P0))*180/pi;
                            if ang > 0
                                candidate = [candidate; thisFrame-window-3+k, file4(pos(po)).id  ang dist(pos(po)) P0 P1];
                            end
                        end
                            
%                             xInFile4 = file4(j).m(1); yInFile4 = file4(j).m(2);
%                             dist = sqrt((xInFile4-xInFile3)^2+(yInFile4-yInFile3)^2);
%                             if dist < minDist
%                                  P0 = [xInFile4 yInFile4]; P1 = [xInFile3 yInFile3];   
%                                  if file4(j).parent+1 >= 1 &&  file4(j).parent+1 <=  size(GMMs_bw{1,k-1}, 2) 
%                                      P2 = [GMMs_bw{1,k-1}(file4(j).parent+1).m(1) GMMs_bw{1,k-1}(file4(j).parent+1).m(2)]; 
%                                  else
%                                      P2 = P0;
%                                  end
%                                  ang = atan2(abs(det([P2-P0;P1-P0])),dot(P2-P0,P1-P0))*180/pi;
%                                  if ang > 0
%                                      candidate = [candidate; thisFrame-window-3+k, file4(j).id  ang dist P0 P1];           
%                                  end
%                             end
%                     end
                        
                    end
                    if isempty(candidate)
                        continue;
                    else                  
                        % /* added on June 12
                        candidate = sortrows(candidate, 1, 'ascend');
                        frameOfGaps = unique(candidate(:, 1));
                        sortedCandidate = zeros(length(frameOfGaps), 8);
                        
                        for fog = 1:length(frameOfGaps)
                            seq = find(candidate(:, 1) == frameOfGaps(fog));
                            tmp = candidate(min(seq):max(seq), :);
                            tmp = sortrows(tmp, 3, 'descend');
                            sortedCandidate(fog, :) = tmp(1, :);
                        end
                        % till here */
                        
                        %                         candidate = sortrows(candidate, 3, 'descend');
                        %                         theCandidate = candidate(1,:);
                        %                         fakeParent = theCandidate(2);
                        
                        % /* added on June 12
                        sortedCandidate = sortrows(sortedCandidate, 1, 'descend');
%                         sortedCandidate
                        theCandidate = sortedCandidate(1,:);
                        fakeParent = theCandidate(1, 2);
%                         disp(['fakeParentId = ' num2str(fakeParent)]);
                        
                        %                         for sc = 1:size(sortedCandidate, 1)
                        %                             theCandidate = sortedCandidate(sc, :);
                        % */
%                         disp(['length(frameOfGaps) = ' num2str(length(frameOfGaps)) ', using the ' num2str(sc) '-th candidate ...']);
                        for frameToGenerateFakeCells = theCandidate(1)+1:thisFrame-1
%                             frameToGenerateFakeCells = theCandidate(1)+1;
%                             disp(['frameToGenerateFakeCells = ' num2str(frameToGenerateFakeCells)]);
                            %                             disp(['frameToGenerateFakeCells = ' num2str(frameToGenerateFakeCells)]);
                            xInFile4 = theCandidate(5); yInFile4 = theCandidate(6);
                            xInFile3 = theCandidate(7); yInFile3 = theCandidate(8);
                            %                         /* added on June 12
                            s = size(GMMs{1, frameToGenerateFakeCells}', 1); % size indicates the id goes sequentially, but if id = 1 3 5 7 9, ...
                            % while size is 5 ... so use s = max(all_ids);
                            maxId = getMaxIdPerFrame(obj, GMMs{1, frameToGenerateFakeCells}');
                            %                            */
%                             disp(['size(GMMs{1, frameToGenerateFakeCells}, 1) = ' num2str(s)]);
%                             disp(['xInFile4=' num2str(xInFile4) ', yInFile4=' num2str(yInFile4)]);
%                             disp(['xInFile3=' num2str(xInFile3) ', yInFile3=' num2str(yInFile3)]);
%                             disp(['cell id = ' num2str(file3(i).id)]);
                            GMMs{1, frameToGenerateFakeCells}(s+1).m = [rand*(xInFile3-xInFile4)+xInFile4, rand*(yInFile3-yInFile4)+yInFile4, 0];
%                             disp(['m: ' num2str(GMMs{1, frameToGenerateFakeCells}(s+1).m(1)) ', ' num2str(GMMs{1, frameToGenerateFakeCells}(s+1).m(2))]);
                            GMMs{1, frameToGenerateFakeCells}(s+1).id = maxId + 1; % was s
                            GMMs{1, frameToGenerateFakeCells}(s+1).lineage = GMMs{1, frameToGenerateFakeCells}(s+1).id;
                            GMMs{1, frameToGenerateFakeCells}(s+1).parent = fakeParent;
                            fakeParent = GMMs{1, frameToGenerateFakeCells}(s+1).id;
                            GMMs{1, frameToGenerateFakeCells}(s+1).svIdx = [-1];
                            GMMs{1, frameToGenerateFakeCells}(s+1).dims = 3;
                            GMMs{1, frameToGenerateFakeCells}(s+1).splitScore = 1;
                            GMMs{1, frameToGenerateFakeCells}(s+1).scale = [1, 1, 1.4];
                            GMMs{1, frameToGenerateFakeCells}(s+1).nu = 1;
                            GMMs{1, frameToGenerateFakeCells}(s+1).beta = 1;
                            GMMs{1, frameToGenerateFakeCells}(s+1).alpha = 1;
                            GMMs{1, frameToGenerateFakeCells}(s+1).nuPrior = 1;
                            GMMs{1, frameToGenerateFakeCells}(s+1).betaPrior = 1;
                            GMMs{1, frameToGenerateFakeCells}(s+1).alphaPrior = 1;
                            GMMs{1, frameToGenerateFakeCells}(s+1).W = [0.0230,0,0,0,0.0230,0,0,0,0.0452];
                            GMMs{1, frameToGenerateFakeCells}(s+1).WPrior = [0.05555, 0, 0, 0, 0.05555, 0, 0, 0, 0.10888];
                            GMMs{1, frameToGenerateFakeCells}(s+1).mPrior =  GMMs{1, frameToGenerateFakeCells}(s+1).m;
                            GMMs{1, frameToGenerateFakeCells}(s+1).distMRFPrior = 1;
                            %                             logFakeCells = [logFakeCells '\nGenerating the fake cells of frame=' num2str(frameToGenerateFakeCells) ' at x=' ...
                            %                                 num2str(GMMs{1, frameToGenerateFakeCells}(s+1).m(1)) ', y=' num2str(GMMs{1, frameToGenerateFakeCells}(s+1).m(2)) ...
                            %                                 ' with id=' num2str(GMMs{1, frameToGenerateFakeCells}(s+1).id) ', its parent id=' num2str(GMMs{1, frameToGenerateFakeCells}(s+1).parent)];
                            %                         end
%                             str = ['Generating the fake cells of frame=' num2str(frameToGenerateFakeCells) ' at x=' ...
%                                 num2str(GMMs{1, frameToGenerateFakeCells}(s+1).m(1)) ', y=' num2str(GMMs{1, frameToGenerateFakeCells}(s+1).m(2)) ...
%                                 ' with id=' num2str(GMMs{1, frameToGenerateFakeCells}(s+1).id) ', its parent id=' num2str(GMMs{1, frameToGenerateFakeCells}(s+1).parent) ...
%                                 ' now cell count is ' num2str(size(GMMs{1, frameToGenerateFakeCells}', 1))];
%                             disp(str);
                        end % */ added on June 12
                        
                        %                         /* added on June 13
                        if frameToGenerateFakeCells + 1 == thisFrame
                            file3(i).parent = GMMs{1, frameToGenerateFakeCells}(s+1).id;
                        else
                            file3(i).parent =  -1;
                        end
                        %                         */ till here
                        
                    end
                end
            end
            
        end
        
        function [file3, GMMs] = linkHorizontal4(obj, GMMs, thisFrame, window, minDist, file2, file3)
            % this function shall consider rescue function ...             
            file2Id = [];
            for j = 1:size(file2, 1)
               file2Id = [file2Id file2(j).id];
            end
            rowNumUsed = [];
            for i = 1:size(file3, 1)
                oldPidx = file3(i).parent;
                if oldPidx >= 0
                    rowInFile2 = find(oldPidx == file2Id);                   
                    if ~isempty(rowInFile2)
                        rowInFile2 = rowInFile2(1);
                        rowNumUsed = [rowNumUsed rowInFile2];
                        file3(i).parent = file2(rowInFile2).parent;
                    end
                end    
                
    
            end            
            rows = 1:size(file2, 1);
            rowNumUnused = setdiff(rows, rowNumUsed);
            %/* adjusted June 19
            xInFile2 = zeros(length(rowNumUnused), 1); yInFile2 = zeros(length(rowNumUnused), 1);
            for j = 1:length(rowNumUnused)
                xInFile2(j, 1) = file2(rowNumUnused(j)).m(1);
                yInFile2(j, 1) = file2(rowNumUnused(j)).m(2);                
            end
            %*/
   
            for i = 1: size(file3, 1)
%                 for j = 1:length(rowNumUnused)
                    if file3(i).parent == -1 
                        xInFile3 = ones(length(rowNumUnused), 1) * file3(i).m(1);
                        yInFile3 = ones(length(rowNumUnused), 1) * file3(i).m(2);
                        dist = sqrt((xInFile2-xInFile3).^2+(yInFile2-yInFile3).^2);                  
                        
                        pos = find(dist < 10);
                        
                        if ~isempty(pos)
                            file3(i).parent = file2(rowNumUnused(pos(1))).parent;
                             rowNumUsed = [rowNumUsed rowNumUnused(pos(1))];
                        end

                    end

            end

            % considering the rescue part ... find out the GMM backwards
            % based on window ...
            c = 0; GMMs_bw = []; framemap = [];
            for frame = max(1, thisFrame-1-window-1):thisFrame-1 % thisFrame-1-window-1
                c = c + 1;
                GMMs_bw{1,c} = GMMs{1, frame};
                framemap = [framemap frame];
            end

            disp(['Frames being considered are ' num2str(thisFrame-1-window-1) ' to ' num2str(thisFrame-1)]);
            
            xs = cell(1, window); ys = cell(1, window); %parents = cell(1, window); ids = cell(1, window);
            
            for k = size(GMMs_bw, 2):-1:1 % starting frome 2, as the first one is used as reference to generate an vector for calculating angles
                file4 = GMMs_bw{1,k};
                xInFile4 = zeros(size(file4, 2), 1); yInFile4 = zeros(size(file4, 2), 1);
                %parent = zeros(size(file4, 2), 1); id = zeros(size(file4, 2), 1);
                parfor j = 1:size(file4, 2) % size(file4, 1)
                    xInFile4(j, 1) = file4(j).m(1);
                    yInFile4(j, 1) = file4(j).m(2);
                    %parent(j, 1) = file4(j).parent;
                   % id(j, 1) = file4(j).id;
                end
                xs{1, k} = xInFile4; ys{1, k} = yInFile4;
                %parents{1, k} = parent; ids{1, k} = id;
            end
            
            orphanCells = [];
            for i = 1:size(file3, 1)
                if file3(i).parent == -1
                    orphanCells = [orphanCells; i];
                end
            end
            tic
            disp([num2str(length(orphanCells)) ' cells need synthetic positions.']);
%             orphanCells'
            for orphan = 1:length(orphanCells) 
                i = orphanCells(orphan);
                    candidate = []; % saves frame id angle dist coordinates ..    
                    flag = 0;
                    if isempty(GMMs_bw)
                        continue;
                    end
                    disp(['Generating synthetic cell ' num2str(orphan) '/' num2str(length(orphanCells))]);
                   
                    for k = size(GMMs_bw, 2)-1:-1:2 % starting frome 2, as the first one is used as reference to generate an vector for calculating angles
%                         disp(['For tmp frame ' num2str(k) ' ...']);
                        
                        file4 = GMMs_bw{1,k};    
                        file5 = GMMs_bw{1, k+1};
                        usedParents = [file5.parent];

                        xInFile4 = xs{1, k}; 
                        yInFile4 = ys{1, k}; 
                        
                        xInFile3 = ones(size(file4, 2), 1)*file3(i).m(1);
                        yInFile3 = ones(size(file4, 2), 1)*file3(i).m(2);
                        dist = sqrt((xInFile4-xInFile3).^2+(yInFile4-yInFile3).^2);
                        
                        pos = find(dist < minDist);        
                        
                        for po = 1:length(pos)
                            if ismember(pos(po)-1, usedParents)
%                                 disp(['POS = ' num2str(pos(po))]);
                                continue;
                            end
                            P0 = [xInFile4(pos(po), 1) yInFile4(pos(po), 1)]; P1 = [xInFile3(1,1) yInFile3(1,1)];
                            if file4(pos(po)).parent+1 >= 1  &&  file4(pos(po)).parent+1 <=  size(GMMs_bw{1,k-1}, 2)                     
                                P2 = [GMMs_bw{1,k-1}(file4(pos(po)).parent+1).m(1) GMMs_bw{1,k-1}(file4(pos(po)).parent+1).m(2)];
%                                 P2 = [GMMs_bw{1,k-1}(tmp3).m(1) GMMs_bw{1,k-1}(tmp3).m(2)];
%                                 ang = atan2(abs(det([P2-P0;P1-P0])),dot(P2-P0,P1-P0))*180/pi;

%                                 disp(['xInFile4 = ' num2str(P0(1)) ', yInFile4 = ' num2str(P0(2))]);
%                                 disp(['xInFile3 = ' num2str(P1(1)) ', yInFile3 = ' num2str(P1(2))]);
%                                 disp(['P2 = ' num2str(P2(1)) ', yP2 = ' num2str(P2(2)) ', ang = ' num2str(ang) ', pos(po) = ' num2str(pos(po)) ' , k = ' num2str(k)]);
%                                 candidate = [candidate; thisFrame-window-3+k, file4(pos(po)).id  ang dist(pos(po)) P0 P1]; % candidate as parents
                                 candidate = [candidate; framemap(k), file4(pos(po)).id dist(pos(po)) P0 P1]; % candidate as parents
                            end                            
                        end    

                     if ~isempty(candidate)
                         flag = 1;
                     end
                     if flag == 1
%                          disp(['Break at k = ' num2str(k)]);
                         break;
                     end
                        
                    end
                    
                    if isempty(candidate)
%                         disp('empty candidate')
                        continue;
                    else                  
                        sortedCandidate = sortrows(candidate, 1, 'descend');
                        theCandidate = sortedCandidate(1,:);
                        fakeParent = theCandidate(1, 2);

                        for frameToGenerateFakeCells = theCandidate(1)+1:thisFrame-1                        
%                            disp(['frameToGenerateFakeCells = ' num2str(frameToGenerateFakeCells)]);
                            xInFile4 = theCandidate(4); yInFile4 = theCandidate(5);
                            xInFile3 = theCandidate(6); yInFile3 = theCandidate(7);
                            %                         /* added on June 12
                            s = size(GMMs{1, frameToGenerateFakeCells}', 1); % size indicates the id goes sequentially, but if id = 1 3 5 7 9, ...
                            % while size is 5 ... so use s = max(all_ids);
                            maxId = getMaxIdPerFrame(obj, GMMs{1, frameToGenerateFakeCells}');
                            %                            */
%                             disp(['size(GMMs{1, frameToGenerateFakeCells}, 1) = ' num2str(s)]);
%                             disp(['xInFile4=' num2str(xInFile4) ', yInFile4=' num2str(yInFile4)]);
%                             disp(['xInFile3=' num2str(xInFile3) ', yInFile3=' num2str(yInFile3)]);
%                             disp(['cell id = ' num2str(file3(i).id)]);
                            GMMs{1, frameToGenerateFakeCells}(s+1).m = [rand*(xInFile3-xInFile4)+xInFile4, rand*(yInFile3-yInFile4)+yInFile4, 0];
%                             disp(['m: ' num2str(GMMs{1, frameToGenerateFakeCells}(s+1).m(1)) ', ' num2str(GMMs{1, frameToGenerateFakeCells}(s+1).m(2))]);
                            GMMs{1, frameToGenerateFakeCells}(s+1).id = maxId + 1; % was s
%                             disp(['maxId = ' num2str(maxId)]);
                            GMMs{1, frameToGenerateFakeCells}(s+1).lineage = GMMs{1, frameToGenerateFakeCells}(s+1).id;
                            GMMs{1, frameToGenerateFakeCells}(s+1).parent = fakeParent;
%                             disp(['fakeParent = ' num2str(fakeParent) ' s = ' num2str(s)]);
                            fakeParent = GMMs{1, frameToGenerateFakeCells}(s+1).id;
                            GMMs{1, frameToGenerateFakeCells}(s+1).svIdx = [-1];
                            GMMs{1, frameToGenerateFakeCells}(s+1).dims = 3;
                            GMMs{1, frameToGenerateFakeCells}(s+1).splitScore = 1;
                            GMMs{1, frameToGenerateFakeCells}(s+1).scale = [1, 1, 1.4];
                            GMMs{1, frameToGenerateFakeCells}(s+1).nu = 1;
                            GMMs{1, frameToGenerateFakeCells}(s+1).beta = 1;
                            GMMs{1, frameToGenerateFakeCells}(s+1).alpha = 1;
                            GMMs{1, frameToGenerateFakeCells}(s+1).nuPrior = 1;
                            GMMs{1, frameToGenerateFakeCells}(s+1).betaPrior = 1;
                            GMMs{1, frameToGenerateFakeCells}(s+1).alphaPrior = 1;
                            GMMs{1, frameToGenerateFakeCells}(s+1).W = [0.0230,0,0,0,0.0230,0,0,0,0.0452];
                            GMMs{1, frameToGenerateFakeCells}(s+1).WPrior = [0.05555, 0, 0, 0, 0.05555, 0, 0, 0, 0.10888];
                            GMMs{1, frameToGenerateFakeCells}(s+1).mPrior =  GMMs{1, frameToGenerateFakeCells}(s+1).m;
                            GMMs{1, frameToGenerateFakeCells}(s+1).distMRFPrior = 1;
                         
                            str = ['Generating the fake cells of frame=' num2str(frameToGenerateFakeCells) ' at x=' ...
                                num2str(GMMs{1, frameToGenerateFakeCells}(s+1).m(1)) ', y=' num2str(GMMs{1, frameToGenerateFakeCells}(s+1).m(2)) ...
                                ' with id=' num2str(GMMs{1, frameToGenerateFakeCells}(s+1).id) ', its parent id=' num2str(GMMs{1, frameToGenerateFakeCells}(s+1).parent) ...
                                ' now cell count is ' num2str(size(GMMs{1, frameToGenerateFakeCells}', 1))];
                            
%                             disp(str);
                        end % */ added on June 12
                        
                        %                         /* added on June 13
                        if frameToGenerateFakeCells + 1 == thisFrame
                            file3(i).parent = GMMs{1, frameToGenerateFakeCells}(s+1).id;
                        else
                            file3(i).parent =  -1;
                        end    
%                         tfill = toc;
%                         disp(['tfill = ' num2str(tfill)]);
                    end
            end
            toc
        end
        
        function [file3, GMMs] = linkHorizontal5(obj, GMMs, thisFrame, window, resWindowSize, minDist, file2, file3)
            % this function shall consider rescue function ... 
            file2Id = [];
            for j = 1:size(file2, 1)
               file2Id = [file2Id file2(j).id];
            end
            rowNumUsed = [];
            for i = 1:size(file3, 1)
                oldPidx = file3(i).parent;
                if oldPidx >= 0
                    rowInFile2 = find(oldPidx == file2Id);                   
                    if ~isempty(rowInFile2)
                        rowInFile2 = rowInFile2(1);
                        rowNumUsed = [rowNumUsed rowInFile2];
                        file3(i).parent = file2(rowInFile2).parent;
                    end
                end
            end            
            rows = 1:size(file2, 1);
            rowNumUnused = setdiff(rows, rowNumUsed);
       
            %/* adjusted June 19
            xInFile2 = zeros(length(rowNumUnused), 1); yInFile2 = zeros(length(rowNumUnused), 1);
            for j = 1:length(rowNumUnused)
                xInFile2(j, 1) = file2(rowNumUnused(j)).m(1);
                yInFile2(j, 1) = file2(rowNumUnused(j)).m(2);
            end
            %*/
            for i = 1: size(file3, 1)
%                 for j = 1:length(rowNumUnused)
                    if file3(i).parent == -1 
                        xInFile3 = ones(length(rowNumUnused), 1) * file3(i).m(1);
                        yInFile3 = ones(length(rowNumUnused), 1) * file3(i).m(2);
                        dist = sqrt((xInFile2-xInFile3).^2+(yInFile2-yInFile3).^2);
                        
                        pos = find(dist < 10);
                        
                        if ~isempty(pos)
                            file3(i).parent = file2(rowNumUnused(pos(1))).parent;
                             rowNumUsed = [rowNumUsed rowNumUnused(pos(1))];
                        end
                        
                    end

            end
   
            % considering the rescue part ... find out the GMM backwards
            % based on window ...
            c = 0;
            GMMs_bw = [];
            for frame = max(1, thisFrame-1-window-1):thisFrame-2 % thisFrame-1-window-1
                c = c + 1;
                GMMs_bw{1,c} = GMMs{1, frame};
            end
            %             disp(['Frames being considered are ' num2str(thisFrame-1-window-1) ' to ' num2str(thisFrame-2)]);
            
            xs = cell(1, window); ys = cell(1, window); parents = cell(1, window); ids = cell(1, window);
            for k = size(GMMs_bw, 2):-1:1 % starting frome 2, as the first one is used as reference to generate an vector for calculating angles
                file4 = GMMs_bw{1,k};
                xInFile4 = zeros(size(file4, 2), 1); yInFile4 = zeros(size(file4, 2), 1); 
                parent = zeros(size(file4, 2), 1); id = zeros(size(file4, 2), 1);
                parfor j = 1:size(file4, 2) % size(file4, 1)
                    xInFile4(j, 1) = file4(j).m(1);
                    yInFile4(j, 1) = file4(j).m(2);
                    parent(j, 1) = file4(j).parent;
                    id(j, 1) = file4(j).id;
                end
                xs{1, k} = xInFile4; ys{1, k} = yInFile4;
                parents{1, k} = parent; ids{1, k} = id;
            end
            
            orphanCells = [];
            for i = 1:size(file3, 1)
                if file3(i).parent == -1
                    orphanCells = [orphanCells; i];
                end
            end
            tic
            disp([num2str(length(orphanCells)) ' cells need synthetic positions.']);
            for orphan = 1:min(200, length(orphanCells))
                i = orphanCells(orphan);           
                    if isempty(GMMs_bw)
                        continue;
                    end
%                     disp(['Generating synthetic cell ' num2str(orphan) '/' num2str(length(orphanCells))]);
                    
                    for k = size(GMMs_bw, 2):-1:2 % starting frome 2, as the first one is used as reference to generate an vector for calculating angles
                        xInFile4 = xs{1, k}; yInFile4 = ys{1, k};
                        xInFile3 = ones(size(xInFile4, 1), 1)*file3(i).m(1);
                        yInFile3 = ones(size(xInFile4, 1), 1)*file3(i).m(2);
                        dist = sqrt((xInFile4-xInFile3).^2+(yInFile4-yInFile3).^2);
                        
                        pos = find(dist < minDist);
                           candidate = []; % saves frame id angle dist coordinates ..                
                        if isempty(pos)
                            continue;
                        end
%                         disp([num2str(length(pos)) ' candidates are found for cell ' num2str(orphan) '/' num2str(length(orphanCells) )]);
                            
                        pos = pos(1);
                        for po = 1:length(pos)
                            P0 = [xInFile4(pos(po), 1) yInFile4(pos(po), 1)]; P1 = [xInFile3(1,1) yInFile3(1,1)];
                            P2 = P0;
                            if parents{1, k}(pos(po))+1 >= 1  &&  parents{1, k}(pos(po))+1 <=  size(parents{1, k-1}, 1)
                                P2 = [xs{1, k-1}(parents{1, k}(pos(po))+1) ys{1, k-1}(parents{1, k}(pos(po))+1)];
%                                 disp('getting P2');
                            end
                            ang = atan2(abs(det([P2-P0;P1-P0])),dot(P2-P0,P1-P0))*180/pi;                    
                            if ang > 0
                                candidate = [candidate; thisFrame-window-3+k, ids{1, k}(pos(po))  ang dist(pos(po)) P0 P1];
                            end
                        end                          
                        
                    end
                    if isempty(candidate)
                        continue;
                    else                  
%                         /* added on June 12
                        candidate = sortrows(candidate, 1, 'ascend');
                        frameOfGaps = unique(candidate(:, 1));
                        sortedCandidate = zeros(length(frameOfGaps), 8);
                        
                        for fog = 1:length(frameOfGaps)
                            seq = find(candidate(:, 1) == frameOfGaps(fog));
                            tmp = candidate(min(seq):max(seq), :);
                            tmp = sortrows(tmp, 3, 'descend');
                            sortedCandidate(fog, :) = tmp(1, :);
                        end
                       
                        sortedCandidate = sortrows(sortedCandidate, 1, 'descend');

                        theCandidate = sortedCandidate(1,:);
                        fakeParent = theCandidate(1, 2);
                       
%                         disp(['length(frameOfGaps) = ' num2str(length(frameOfGaps)) ', using the ' num2str(sc) '-th candidate ...']);
%                         tic
                        for frameToGenerateFakeCells = theCandidate(1)+1:thisFrame-1                        
                            %                             disp(['frameToGenerateFakeCells = ' num2str(frameToGenerateFakeCells)]);
                            xInFile4 = theCandidate(5); yInFile4 = theCandidate(6);
                            xInFile3 = theCandidate(7); yInFile3 = theCandidate(8);
                            %                         /* added on June 12
                            s = size(GMMs{1, frameToGenerateFakeCells}', 1); % size indicates the id goes sequentially, but if id = 1 3 5 7 9, ...
                            % while size is 5 ... so use s = max(all_ids);
                            
                            maxId = getMaxIdPerFrame(obj, GMMs{1, frameToGenerateFakeCells}');

                            %                            */
%                             disp(['size(GMMs{1, frameToGenerateFakeCells}, 1) = ' num2str(s)]);
%                             disp(['xInFile4=' num2str(xInFile4) ', yInFile4=' num2str(yInFile4)]);
%                             disp(['xInFile3=' num2str(xInFile3) ', yInFile3=' num2str(yInFile3)]);
%                             disp(['cell id = ' num2str(file3(i).id)]);
                            GMMs{1, frameToGenerateFakeCells}(s+1).m = [rand*(xInFile3-xInFile4)+xInFile4, rand*(yInFile3-yInFile4)+yInFile4, 0];
%                             disp(['m: ' num2str(GMMs{1, frameToGenerateFakeCells}(s+1).m(1)) ', ' num2str(GMMs{1, frameToGenerateFakeCells}(s+1).m(2))]);
                            GMMs{1, frameToGenerateFakeCells}(s+1).id = maxId + 1; % was s
                            GMMs{1, frameToGenerateFakeCells}(s+1).lineage = GMMs{1, frameToGenerateFakeCells}(s+1).id;
                            GMMs{1, frameToGenerateFakeCells}(s+1).parent = fakeParent;
                            fakeParent = GMMs{1, frameToGenerateFakeCells}(s+1).id;
                            GMMs{1, frameToGenerateFakeCells}(s+1).svIdx = [-1];
                            GMMs{1, frameToGenerateFakeCells}(s+1).dims = 3;
                            GMMs{1, frameToGenerateFakeCells}(s+1).splitScore = 1;
                            GMMs{1, frameToGenerateFakeCells}(s+1).scale = [1, 1, 1.4];
                            GMMs{1, frameToGenerateFakeCells}(s+1).nu = 1;
                            GMMs{1, frameToGenerateFakeCells}(s+1).beta = 1;
                            GMMs{1, frameToGenerateFakeCells}(s+1).alpha = 1;
                            GMMs{1, frameToGenerateFakeCells}(s+1).nuPrior = 1;
                            GMMs{1, frameToGenerateFakeCells}(s+1).betaPrior = 1;
                            GMMs{1, frameToGenerateFakeCells}(s+1).alphaPrior = 1;
                            GMMs{1, frameToGenerateFakeCells}(s+1).W = [0.0230,0,0,0,0.0230,0,0,0,0.0452];
                            GMMs{1, frameToGenerateFakeCells}(s+1).WPrior = [0.05555, 0, 0, 0, 0.05555, 0, 0, 0, 0.10888];
                            GMMs{1, frameToGenerateFakeCells}(s+1).mPrior =  GMMs{1, frameToGenerateFakeCells}(s+1).m;
                            GMMs{1, frameToGenerateFakeCells}(s+1).distMRFPrior = 1;
                         
%                             str = ['Generating the fake cells of frame=' num2str(frameToGenerateFakeCells) ' at x=' ...
%                                 num2str(GMMs{1, frameToGenerateFakeCells}(s+1).m(1)) ', y=' num2str(GMMs{1, frameToGenerateFakeCells}(s+1).m(2)) ...
%                                 ' with id=' num2str(GMMs{1, frameToGenerateFakeCells}(s+1).id) ', its parent id=' num2str(GMMs{1, frameToGenerateFakeCells}(s+1).parent) ...
%                                 ' now cell count is ' num2str(size(GMMs{1, frameToGenerateFakeCells}', 1))];
%                             disp(str);
                        end % */ added on June 12
                        
                        %                         /* added on June 13
                        if frameToGenerateFakeCells + 1 == thisFrame
                            file3(i).parent = GMMs{1, frameToGenerateFakeCells}(s+1).id;
                        else
                            file3(i).parent =  -1;
                        end    
%                         tfill = toc;
%                         disp(['tfill = ' num2str(tfill)]);
                    end
            end
            toc
           end
        
        

        function [GMMs] = updateGMMwithSyntheticCells(obj, GMMs, frameToGenerateFakeCells, m, id, fakeParent, s)
            GMMs{1, frameToGenerateFakeCells}(s+1).m = m;
%                             disp(['m: ' num2str(GMMs{1, frameToGenerateFakeCells}(s+1).m(1)) ', ' num2str(GMMs{1, frameToGenerateFakeCells}(s+1).m(2))]);
                            GMMs{1, frameToGenerateFakeCells}(s+1).id = id; % was s
                            GMMs{1, frameToGenerateFakeCells}(s+1).lineage = GMMs{1, frameToGenerateFakeCells}(s+1).id;
                            GMMs{1, frameToGenerateFakeCells}(s+1).parent = fakeParent;
%                             fakeParent = GMMs{1, frameToGenerateFakeCells}(s+1).id;
                            
                            GMMs{1, frameToGenerateFakeCells}(s+1).svIdx = [-1];
                            GMMs{1, frameToGenerateFakeCells}(s+1).dims = 3;
                            GMMs{1, frameToGenerateFakeCells}(s+1).splitScore = 1;
                            GMMs{1, frameToGenerateFakeCells}(s+1).scale = [1, 1, 1.4];
                            GMMs{1, frameToGenerateFakeCells}(s+1).nu = 1;
                            GMMs{1, frameToGenerateFakeCells}(s+1).beta = 1;
                            GMMs{1, frameToGenerateFakeCells}(s+1).alpha = 1;
                            GMMs{1, frameToGenerateFakeCells}(s+1).nuPrior = 1;
                            GMMs{1, frameToGenerateFakeCells}(s+1).betaPrior = 1;
                            GMMs{1, frameToGenerateFakeCells}(s+1).alphaPrior = 1;
                            GMMs{1, frameToGenerateFakeCells}(s+1).W = [0.0230,0,0,0,0.0230,0,0,0,0.0452];
                            GMMs{1, frameToGenerateFakeCells}(s+1).WPrior = [0.05555, 0, 0, 0, 0.05555, 0, 0, 0, 0.10888];
                            GMMs{1, frameToGenerateFakeCells}(s+1).mPrior =  GMMs{1, frameToGenerateFakeCells}(s+1).m;
                            GMMs{1, frameToGenerateFakeCells}(s+1).distMRFPrior = 1;
        end
            
        function ori_idxMatrix2 = matrix2xmlMamut_seudofill(obj, ori_idxMatrix, ori_coordMatrix2X, ori_coordMatrix2Y, ori_SVidxMatrix, xmlPathRearranged)
            ori_idxMatrix2 = ori_idxMatrix; 
            for i = size(ori_idxMatrix2, 2)-1:-1:2
                counter_discontinuity = 0;
                for j = 1:size(ori_idxMatrix2, 1)
                    if  ori_idxMatrix2(j, i) == 0 && ori_idxMatrix2(j, i - 1)*ori_idxMatrix2(j, i + 1) > 0
                        counter_discontinuity = counter_discontinuity + 1;
                        ori_idxMatrix2(j, i) = max(ori_idxMatrix(:, i)) + counter_discontinuity;
                    end
                end
            end
            for col = size(ori_idxMatrix, 2):-1:1
                doc = com.mathworks.xml.XMLUtils.createDocument('document');
                document = doc.getDocumentElement;
%                 counter_discontinuity = 0;
                for row = 1:size(ori_idxMatrix, 1)
                    gmm = doc.createElement('GaussianMixtureModel');   
                    if ori_idxMatrix(row, col) > 0
                        gmm.setAttribute('id',num2str(ori_idxMatrix(row, col)-1));
                        gmm.setAttribute('lineage',num2str(ori_idxMatrix(row, col)-1));
                        pIdx = ori_idxMatrix(row, max(1, col-1))-1;
                        if pIdx == -1 && ori_idxMatrix(row, max(1, col-2)) > 0 % dealing with discontinuity ...
                            pIdx = ori_idxMatrix2(row, max(1, col-1))-1;
                        end
                        gmm.setAttribute('parent',num2str(pIdx));
                        gmm.setAttribute('dims', num2str(3));
                        gmm.setAttribute('splitScore', num2str(1));
                        gmm.setAttribute('scale', [num2str(1) ' ' num2str(1) ' ' num2str(1.4) ' ']); % for full field of view, aspect ratio = 1.4
                        gmm.setAttribute('nu', num2str(1)); % random value
                        gmm.setAttribute('beta', num2str(1));
                        gmm.setAttribute('alpha', num2str(1));
                        gmm.setAttribute('m', [num2str(ori_coordMatrix2X(row, col)) ' ' num2str(ori_coordMatrix2Y(row, col)) ' ' num2str(0) ' ']);
                        W = [0.0230,0,0,0,0.0230,0,0,0,0.0452]; str = '';
                        for i = 1:size(W, 2)
                            str = [str num2str(W(i)) ' '];
                        end
                        gmm.setAttribute('W',str);
                        gmm.setAttribute('nuPrior',num2str(1));
                        gmm.setAttribute('betaPrior',num2str(1));
                        gmm.setAttribute('alphaPrior',num2str(1));
                        gmm.setAttribute('distMRFPrior',num2str(1));
                        gmm.setAttribute('mPrior', [num2str(ori_coordMatrix2X(row, col)) ' ' num2str(ori_coordMatrix2Y(row, col)) ' ' num2str(0) ' ']);
                        WPrior = [0.05555, 0, 0, 0, 0.05555, 0, 0, 0, 0.10888];
                        str = '';
                        for i = 1:size(WPrior, 2)
                            str = [str num2str(WPrior(i)) ' '];
                        end
                        gmm.setAttribute('WPrior',str);
                        str = '';
                        for i = 1:size(ori_SVidxMatrix{row, col}, 2)
                            str = [str num2str(ori_SVidxMatrix{row, col}(i)-1) ' '];
                        end
                        gmm.setAttribute('svIdx',str);
                        document.appendChild(gmm);
                    else
                        % for case when discontinuity occurs in the middle course after appending, 
                        % such as [1 1 2 3 5 4 0 1 5 9 2] of a same track ...
                        if col > 1 && col < size(ori_idxMatrix, 2) && ori_idxMatrix(row, col - 1)*ori_idxMatrix(row, col + 1) > 0
%                             counter_discontinuity = counter_discontinuity + 1;
                            %fakeId = max(ori_idxMatrix(:, col)) + counter_discontinuity; % guess an id which does not interfere with other cells' ids...
                            fakeId = ori_idxMatrix2(row, col);
                            gmm.setAttribute('id',num2str(fakeId - 1));
                            disp(['Creating fake cells of row ' num2str(row) ' col ' num2str(col) ' with fakeId = ' num2str(fakeId)]);
                            gmm.setAttribute('lineage',num2str(fakeId - 1));
                            gmm.setAttribute('parent',num2str(ori_idxMatrix(row, max(1, col-1))-1));
                            gmm.setAttribute('dims', num2str(3));
                            gmm.setAttribute('splitScore', num2str(1));
                            gmm.setAttribute('scale', [num2str(1) ' ' num2str(1) ' ' num2str(1.4) ' ']); % for full field of view, aspect ratio = 1.4
                            gmm.setAttribute('nu', num2str(1)); % random value
                            gmm.setAttribute('beta', num2str(1));
                            gmm.setAttribute('alpha', num2str(1));
                            x1 = ori_coordMatrix2X(row, col-1); x2 = ori_coordMatrix2X(row, col+1); 
                            y1 = ori_coordMatrix2Y(row, col-1); y2 = ori_coordMatrix2Y(row, col+1);
                            x = rand*(x2-x1)+x1; y = rand*(y2-y1)+y1;
                            gmm.setAttribute('m', [num2str(x) ' ' num2str(y) ' ' num2str(0) ' ']);
                            W = [0.0230,0,0,0,0.0230,0,0,0,0.0452]; str = '';
                            for i = 1:size(W, 2)
                                str = [str num2str(W(i)) ' '];
                            end
                            gmm.setAttribute('W',str);
                            gmm.setAttribute('nuPrior',num2str(1));
                            gmm.setAttribute('betaPrior',num2str(1));
                            gmm.setAttribute('alphaPrior',num2str(1));
                            gmm.setAttribute('distMRFPrior',num2str(1));
                            gmm.setAttribute('mPrior', [num2str(ori_coordMatrix2X(row, col)) ' ' num2str(ori_coordMatrix2Y(row, col)) ' ' num2str(0) ' ']);
                            WPrior = [0.05555, 0, 0, 0, 0.05555, 0, 0, 0, 0.10888];
                            str = '';
                            for i = 1:size(WPrior, 2)
                                str = [str num2str(WPrior(i)) ' '];
                            end
                            gmm.setAttribute('WPrior',str);
                            str = '10000 ';
                            gmm.setAttribute('svIdx',str);
                            document.appendChild(gmm);
                        end        
                    end     
                end
                
                xmlFilename = [xmlPathRearranged '\GMEMfinalResult_frame' num2str(col-1,'%.4d') '.xml'];
                xmlwrite(xmlFilename,doc);
                clear doc;
                disp(['Finished writing frame ' num2str(col-1)]);
            end
        end
        
        function  matrix2xmlMamut(obj, ori_idxMatrix, ori_coordMatrix2X, ori_coordMatrix2Y, ori_SVidxMatrix, xmlPathRearranged)           
            for col = size(ori_idxMatrix, 2):-1:1
                doc = com.mathworks.xml.XMLUtils.createDocument('document');
                document = doc.getDocumentElement;
                for row = 1:size(ori_idxMatrix, 1)
                    gmm = doc.createElement('GaussianMixtureModel');
                    if ori_idxMatrix(row, col) > 0
                        gmm.setAttribute('id',num2str(ori_idxMatrix(row, col)-1));
                        gmm.setAttribute('lineage',num2str(ori_idxMatrix(row, col)-1));
                        pIdx = ori_idxMatrix(row, max(1, col-1))-1;
                        gmm.setAttribute('parent',num2str(pIdx));
                        gmm.setAttribute('dims', num2str(3));
                        gmm.setAttribute('splitScore', num2str(1));
                        gmm.setAttribute('scale', [num2str(1) ' ' num2str(1) ' ' num2str(1.4) ' ']); % for full field of view, aspect ratio = 1.4
                        gmm.setAttribute('nu', num2str(1)); % random value
                        gmm.setAttribute('beta', num2str(1));
                        gmm.setAttribute('alpha', num2str(1));
                        gmm.setAttribute('m', [num2str(ori_coordMatrix2X(row, col)) ' ' num2str(ori_coordMatrix2Y(row, col)) ' ' num2str(0) ' ']);
                        W = [0.0230,0,0,0,0.0230,0,0,0,0.0452]; str = '';
                        for i = 1:size(W, 2)
                            str = [str num2str(W(i)) ' '];
                        end
                        gmm.setAttribute('W',str);
                        gmm.setAttribute('nuPrior',num2str(1));
                        gmm.setAttribute('betaPrior',num2str(1));
                        gmm.setAttribute('alphaPrior',num2str(1));
                        gmm.setAttribute('distMRFPrior',num2str(1));
                        gmm.setAttribute('mPrior', [num2str(ori_coordMatrix2X(row, col)) ' ' num2str(ori_coordMatrix2Y(row, col)) ' ' num2str(0) ' ']);
                        WPrior = [0.05555, 0, 0, 0, 0.05555, 0, 0, 0, 0.10888];
                        str = '';
                        for i = 1:size(WPrior, 2)
                            str = [str num2str(WPrior(i)) ' '];
                        end
                        gmm.setAttribute('WPrior',str);
                        str = '';
                        for i = 1:size(ori_SVidxMatrix{row, col}, 2)
                            str = [str num2str(ori_SVidxMatrix{row, col}(i)-1) ' '];
                        end
                        gmm.setAttribute('svIdx',str);
                        document.appendChild(gmm);                        
                    end
                end                
                xmlFilename = [xmlPathRearranged '\GMEMfinalResult_frame' num2str(col-1,'%.4d') '.xml'];
                xmlwrite(xmlFilename,doc);
                clear doc;
                disp(['Finished writing frame ' num2str(col-1)]);
            end
        end
        
        function GMM_Res = matrix2gmm(obj, ori_idxMatrix, ori_coordMatrix2X, ori_coordMatrix2Y, ori_SVidxMatrix, GMM_Res, se)
         
            for col = size(ori_idxMatrix, 2):-1:1
                counter = 0;
                for row = 1:size(ori_idxMatrix, 1)
                    if ori_idxMatrix(row, col) > 0
                        counter = counter + 1;
                        GMM_Res{se, col}(counter).m = [ori_coordMatrix2X(row, col) ori_coordMatrix2Y(row, col) 0];
                        GMM_Res{se, col}(counter).id = ori_idxMatrix(row, col)-1;
                        GMM_Res{se, col}(counter).lineage = ori_idxMatrix(row, col)-1;
                        GMM_Res{se, col}(counter).parent = ori_idxMatrix(row, max(1, col-1))-1;
                        GMM_Res{se, col}(counter).svIdx = ori_SVidxMatrix{row, col} - 1;
                        GMM_Res{se, col}(counter).dims = 3;
                        GMM_Res{se, col}(counter).splitScore = 1;
                        GMM_Res{se, col}(counter).scale = [1 1 1.4];
                        GMM_Res{se, col}(counter).nu = 1;
                        GMM_Res{se, col}(counter).beta = 1;
                        GMM_Res{se, col}(counter).alpha = 1;
                        GMM_Res{se, col}(counter).W = [0.0230,0,0,0,0.0230,0,0,0,0.0452];
                        GMM_Res{se, col}(counter).nuPrior = 1;
                        GMM_Res{se, col}(counter).betaPrior = 1;
                        GMM_Res{se, col}(counter).alphaPrior = 1;
                        GMM_Res{se, col}(counter).mPrior = [ori_coordMatrix2X(row, col) ori_coordMatrix2Y(row, col) 0];
                        GMM_Res{se, col}(counter).WPrior = [0.05555, 0, 0, 0, 0.05555, 0, 0, 0, 0.10888];
                        GMM_Res{se, col}(counter).distMRFPrior = 1;
                    end
                end
            end

        end        
        
             
           function GMM_Res = matrix2gmm2(obj, ori_idxMatrix, ori_coordMatrix2X, ori_coordMatrix2Y, ori_SVidxMatrix)
            GMM_Res = cell(1, size(ori_idxMatrix, 2));
            for col = size(ori_idxMatrix, 2):-1:1
                counter = 0;
                for row = 1:size(ori_idxMatrix, 1)
                    if ori_idxMatrix(row, col) > 0
                        counter = counter + 1;
                        GMM_Res{1, col}(counter).m = [ori_coordMatrix2X(row, col) ori_coordMatrix2Y(row, col) 0];
                        GMM_Res{1, col}(counter).id = ori_idxMatrix(row, col)-1;
                        GMM_Res{1, col}(counter).lineage = ori_idxMatrix(row, col)-1;
                        GMM_Res{1, col}(counter).parent = ori_idxMatrix(row, max(1, col-1))-1;
                        GMM_Res{1, col}(counter).svIdx = ori_SVidxMatrix{row, col} - 1;
                        
                        GMM_Res{1, col}(counter).dims = 3;
                        GMM_Res{1, col}(counter).splitScore = 1;
                        GMM_Res{1, col}(counter).scale = [1 1 1.4];
                        GMM_Res{1, col}(counter).nu = 1;
                        GMM_Res{1, col}(counter).beta = 1;
                        GMM_Res{1, col}(counter).alpha = 1;
                        GMM_Res{1, col}(counter).W = [0.0230,0,0,0,0.0230,0,0,0,0.0452];
                        GMM_Res{1, col}(counter).nuPrior = 1;
                        GMM_Res{1, col}(counter).betaPrior = 1;
                        GMM_Res{1, col}(counter).alphaPrior = 1;
                        GMM_Res{1, col}(counter).mPrior = [ori_coordMatrix2X(row, col) ori_coordMatrix2Y(row, col) 0];
                        GMM_Res{1, col}(counter).WPrior = [0.05555, 0, 0, 0, 0.05555, 0, 0, 0, 0.10888];
                        GMM_Res{1, col}(counter).distMRFPrior = 1;
                    end
                end
            end

           end        
        
        
        function GMMs_tmp = getGMMsOfSpecifiedFrames(obj, GMMs_mamut, t1, t2) 
            % extracting GMMs_tmp from GMMs_mamut
            GMMs_tmp = cell(1, t2-t1+1);
            for t = t1:t2
                GMMs_tmp{1, t-t1+1} = GMMs_mamut{1, t};
            end
        end
        
        function GMMs_mamut = updateGMMsOfSpecificFrames(obj, GMMs_tmp, GMMs_mamut, t1, t2)
            % updating/writing GMMs_tmp to GMMs_mamut
            for t = t1:t2
                GMMs_mamut{1, t1+t-1} = GMMs_tmp{1, t};
            end
        end
        
        function GMMs_mamut = reOrderIdx(obj, GMMs_mamut, t1, t2)            
            for g = t2:-1:t1
                GMMs_latter = GMMs_mamut{1, g}';
                if g > 1
                    GMMs_former = GMMs_mamut{1, g-1}';
                    rowId_gmmId_former = zeros(size(GMMs_former, 1), 1);
                    for b = 1:size(GMMs_former, 1)
                        rowId_gmmId_former(b) = GMMs_former(b).id;
                    end
                    for a = 1:size(GMMs_latter, 1)
                        if GMMs_latter(a).id ~= a - 1
                            GMMs_latter(a).id = a - 1;
                            GMMs_latter(a).lineage = a - 1;
                        end
                        if GMMs_latter(a).parent > -1
                            tmp = GMMs_latter(a).parent;
                            pos = find(tmp == rowId_gmmId_former);
                            if ~isempty(pos)
                                GMMs_latter(a).parent = pos(1)-1;
                            end
                            %                             GMMs_latter(a).parent = find(tmp == row_gmmId_former)-1;
                        end
                    end
                end
                GMMs_mamut{1,g} = GMMs_latter';
            end
            for gg = 1:size(GMMs_mamut{1,t1}, 2)
                GMMs_mamut{1,t1}(gg).parent = -1;
                GMMs_mamut{1,t1}(gg).id = gg-1;
                GMMs_mamut{1,t1}(gg).lineage = gg-1;
            end

        end
        
        function [file3, GMMs] = linkHorizontalBackward(obj, GMMs, thisFrame, window, resWindowSize, minDist, file2, file3)
           [~, idx] = unique([file2.id].', 'rows', 'stable');  
           file2 = file2(idx);
            file2svIdx = [file2.svIdx];        
            orphanCells = [];
%             assign = [];
            for i = 1:size(file3, 2)
                pos = find(file3(i).svIdx == file2svIdx);
                if ~isempty(pos)
                    if length(pos) == 1
                        file3(i).parent = file2(pos).parent;
%                     else   
%                         if file2(pos(1)).parent ~= file2(pos(2)).parent                          
%                             assign = [assign; i pos];
%                         end
                    end
                else
                    file3(i).parent = -1;
                end           
%                 file3(i).parent
                if file3(i).parent == -1
                    orphanCells = [orphanCells; i];
                end
            end
%             %
%             if ~isempty(assign)
%                 uq = unique(assign(:, 1));
%                 rowNum = size(uq , 1);
%                 assign2 = cell(rowNum, 2);
%                 for rn = 1:rowNum
%                     assign2{rn, 1} = uq(rn);
%                 end
%                 for a = 1:size(assign, 1)
%                     a1 = assign(a, 1);
%                     a2 = assign(a, 2);
%                     posa = find(a1 == uq);
%                     assign2{posa, 2} = [assign2{posa, 2} a2];
%                 end
%                 order = [];
%                 for a = 1:size(assign2, 1)
%                     tmp1 = assign2{a, 2};
%                     flag = 0;
%                     for b = 1:size(assign2, 1)
%                         if a == b
%                             continue;
%                         end
%                         tmp2 = assign2{b, 2};
%                         if isequal(tmp1, tmp2)
%                             flag = b;
%                         end
%                     end
%                     order = [order; min(a, flag) max(a, flag)];
%                 end
%                 order = unique(order, 'rows');
%                 assign3 = cell(1, 2);
%                 for O = 1:size(order, 1)
%                     if order(O, 1) > 0
%                         assign3{O, 1} = [assign2{order(O, 1), 1} assign2{order(O, 2), 1}];
%                     else
%                         assign3{O, 1} = [assign2{order(O, 2), 1}];
%                     end
%                     assign3{O, 2} = [assign2{order(O, 2), 2}];
%                 end
%                 % assigning the daughters (of the same cell) to parents, this applies when
%                 % we correct cell lineages
%                 for aa = 1:size(assign3, 1)
%                     daughters = sort(assign3{aa, 1}, 'ascend');
%                     parents = sort(assign3{aa, 2}, 'ascend');
%                     if size(daughters) == size(parents)
%                         for d = 1:length(daughters)
%                             file3(daughters(d)).parent = file2(parents(d)).parent;
%                         end
% %                     else
% %                         file3(daughters(1)).parent = file2(parents(1)).parent;
%                     end
%                 end
%             end

%             orphanCells

            disp([num2str(length(orphanCells)) ' cells need synthetic positions.']);

            % Filling up the missing cell gap
            % considering the rescue part ... find out the GMM backwards based on window ...
            c = 0;
            GMMs_bw = [];
            for frame = max(1, thisFrame-1-window-1):thisFrame-2 % thisFrame-1-window-1
                c = c + 1;
                GMMs_bw{1,c} = GMMs{1, frame};
            end
            disp(['Frames being considered with fake cells are ' num2str(thisFrame-1-window) ' to ' num2str(thisFrame-2)]);
            xs = cell(1, window); ys = cell(1, window); %parents = cell(1, window); ids = cell(1, window);
            for k = size(GMMs_bw, 2):-1:1 % starting frome 2, as the first one is used as reference to generate an vector for calculating angles
                file4 = GMMs_bw{1,k}';
                xInFile4 = zeros(size(file4, 2), 1); yInFile4 = zeros(size(file4, 2), 1);
                %parent = zeros(size(file4, 2), 1); id = zeros(size(file4, 2), 1);
                parfor j = 1:size(file4, 2) % size(file4, 1)
                    xInFile4(j, 1) = file4(j).m(1);
                    yInFile4(j, 1) = file4(j).m(2);
                    %parent(j, 1) = file4(j).parent;
                    % id(j, 1) = file4(j).id;
                end
                xs{1, k} = xInFile4; ys{1, k} = yInFile4;
                %parents{1, k} = parent; ids{1, k} = id;
            end
            
            for orphan = 1:length(orphanCells)
                o = orphanCells(orphan);                
                if isempty(GMMs_bw)
%                     disp(['GMMs_bw No candidate found for orphancell with id ' num2str(o)]);
                    continue;
                end
                candidate = []; % saves frame id angle dist coordinates ..
                for k = size(GMMs_bw, 2):-1:2 % starting frome 2, as the first one is used as reference to generate an vector for calculating angles
                    file4 = GMMs_bw{1,k};
                    
                    xInFile4 = xs{1, k};
                    yInFile4 = ys{1, k};
                    
                    xInFile3 = ones(size(file4, 2), 1)*file3(o).m(1);
                    yInFile3 = ones(size(file4, 2), 1)*file3(o).m(2);
                    dist = sqrt((xInFile4-xInFile3).^2+(yInFile4-yInFile3).^2);
      
                    pos = find(dist < minDist);
                    if isempty(pos)                        
                        disp(['K = ' num2str(k) ' No candidate found for orphancell with id ' num2str(o) ' at (x, y) = ' num2str(file3(o).m(1)) ', ' num2str(file3(o).m(2))]);
                    continue;
                    end
                    
                    for po = 1:length(pos)
                        P0 = [xInFile4(pos(po), 1) yInFile4(pos(po), 1)]; P1 = [xInFile3(1,1) yInFile3(1,1)];
                        disp(['Found parent file4 with row ' num2str(pos(po)) ' thisFrame = ' num2str(thisFrame) ', K = ' num2str(k)]);
                        P2 = P0;
                        if file4(pos(po)).parent >= 0 
%                             &&  file4(pos(po)).parent+1 <=  size(GMMs_bw{1,k-1}, 2)
                            findrow = find(file4(pos(po)).parent == [GMMs_bw{1,k-1}.id]);
                           
                            if ~isempty(findrow)
                                 findrow = findrow(1);
                                P2 = [GMMs_bw{1,k-1}(findrow).m(1) GMMs_bw{1,k-1}(findrow).m(2)];
%                                 P2 = [GMMs_bw{1,k-1}(file4(pos(po)).parent+1).m(1) GMMs_bw{1,k-1}(file4(pos(po)).parent+1).m(2)];
                            end
                        end
                        disp(['K = ' num2str(k) '           P0=' num2str(P0) ',           P1=' num2str(P1) ',           P2=' num2str(P2)])
                        ang = atan2(abs(det([P2-P0;P1-P0])),dot(P2-P0,P1-P0))*180/pi;
                        if ang > 0
                            candidate = [candidate; thisFrame-window-3+k, file4(pos(po)).id  ang dist(pos(po)) P0 P1];
                        end                        
                    end                    
                end 
                
                if isempty(candidate)
                     disp(['No candidate found for orphancell with id ' num2str(o)]);
                    continue;
                else
                    % /* added on June 12
                    candidate = sortrows(candidate, 1, 'ascend');
                    frameOfGaps = unique(candidate(:, 1));
                    sortedCandidate = zeros(length(frameOfGaps), 8);
                    
                    for fog = 1:length(frameOfGaps)
                        seq = find(candidate(:, 1) == frameOfGaps(fog));
                        tmp = candidate(min(seq):max(seq), :);
                        tmp = sortrows(tmp, 3, 'descend');
                        sortedCandidate(fog, :) = tmp(1, :);
                    end
                    
                    sortedCandidate = sortrows(sortedCandidate, 1, 'descend');
                    
                    theCandidate = sortedCandidate(1,:);
                    fakeParent = theCandidate(1, 2);
                    
                    for frameToGenerateFakeCells = theCandidate(1)+1:thisFrame-1
%                                                     disp(['frameToGenerateFakeCells = ' num2str(frameToGenerateFakeCells)]);
                        xInFile4 = theCandidate(5); yInFile4 = theCandidate(6);
                        xInFile3 = theCandidate(7); yInFile3 = theCandidate(8);
                        %                         /* added on June 12
                        s = size(GMMs{1, frameToGenerateFakeCells}, 1); % size indicates the id goes sequentially, but if id = 1 3 5 7 9, ...
                        % while size is 5 ... so use s = max(all_ids);
                        maxId = getMaxIdPerFrame(obj, GMMs{1, frameToGenerateFakeCells});
%                                                    */
%                                                     disp(['size(GMMs{1, frameToGenerateFakeCells}, 1) = ' num2str(s)]);
%                                                     disp(['xInFile4=' num2str(xInFile4) ', yInFile4=' num2str(yInFile4)]);
%                                                     disp(['xInFile3=' num2str(xInFile3) ', yInFile3=' num2str(yInFile3)]);
%                                                     disp(['cell id = ' num2str(file3(i).id)]);
                        GMMs{1, frameToGenerateFakeCells}(s+1).m = [rand*(xInFile3-xInFile4)+xInFile4, rand*(yInFile3-yInFile4)+yInFile4, 0];
%                                                     disp(['m: ' num2str(GMMs{1, frameToGenerateFakeCells}(s+1).m(1)) ', ' num2str(GMMs{1, frameToGenerateFakeCells}(s+1).m(2))]);
                        GMMs{1, frameToGenerateFakeCells}(s+1).id = maxId + 1; % was s
                        GMMs{1, frameToGenerateFakeCells}(s+1).lineage = GMMs{1, frameToGenerateFakeCells}(s+1).id;
                        GMMs{1, frameToGenerateFakeCells}(s+1).parent = fakeParent;
                        fakeParent = GMMs{1, frameToGenerateFakeCells}(s+1).id;
                        GMMs{1, frameToGenerateFakeCells}(s+1).svIdx = [-1];
                        GMMs{1, frameToGenerateFakeCells}(s+1).dims = 3;
                        GMMs{1, frameToGenerateFakeCells}(s+1).splitScore = 1;
                        GMMs{1, frameToGenerateFakeCells}(s+1).scale = [1, 1, 1.4];
                        GMMs{1, frameToGenerateFakeCells}(s+1).nu = 1;
                        GMMs{1, frameToGenerateFakeCells}(s+1).beta = 1;
                        GMMs{1, frameToGenerateFakeCells}(s+1).alpha = 1;
                        GMMs{1, frameToGenerateFakeCells}(s+1).nuPrior = 1;
                        GMMs{1, frameToGenerateFakeCells}(s+1).betaPrior = 1;
                        GMMs{1, frameToGenerateFakeCells}(s+1).alphaPrior = 1;
                        GMMs{1, frameToGenerateFakeCells}(s+1).W = [0.0230,0,0,0,0.0230,0,0,0,0.0452];
                        GMMs{1, frameToGenerateFakeCells}(s+1).WPrior = [0.05555, 0, 0, 0, 0.05555, 0, 0, 0, 0.10888];
                        GMMs{1, frameToGenerateFakeCells}(s+1).mPrior =  GMMs{1, frameToGenerateFakeCells}(s+1).m;
                        GMMs{1, frameToGenerateFakeCells}(s+1).distMRFPrior = 1;
                        
                        str = ['Generating the fake cells of frame=' num2str(frameToGenerateFakeCells) ' at x=' ...
                            num2str(GMMs{1, frameToGenerateFakeCells}(s+1).m(1)) ', y=' num2str(GMMs{1, frameToGenerateFakeCells}(s+1).m(2)) ...
                            ' with id=' num2str(GMMs{1, frameToGenerateFakeCells}(s+1).id) ', its parent id=' num2str(GMMs{1, frameToGenerateFakeCells}(s+1).parent) ...
                            ' now cell count is ' num2str(size(GMMs{1, frameToGenerateFakeCells}, 1))];
                        disp(str);
                    end % */ added on June 12
                    if frameToGenerateFakeCells + 1 == thisFrame
                        file3(o).parent = GMMs{1, frameToGenerateFakeCells}(s+1).id;
                    else
                        file3(o).parent =  -1;
                    end
                    
                end               
            end         
            
        end
        
        
        
        
    end
end