classdef cellTracks
    
    properties
        tracks; % all tracks including the 3 types below
        fullTracks; % id(t) > 0 for every t in [0, T]
        deadTracks; % id(t) > 0 for every t in [t_0, t_T] and t_0 > 0 and t_T< T
        partialTracks; % id(t) > 0 for every t in [t_0, T] and t_0 > 0
        lineages; % group tracks of the same lineage together, saved in cell array
        lineagesTB; % representing lineages in the format of table
        mapLineage2Tracks; % a cell array, every entry e.g,{1, i} stores all the tracks belonging to the i-th lineage...
        tb_div; % this stores detected divisions, used for detecting false-positives.
                %col 1: dividing time; col 2: id of daughter 1 at dividing time;
                % col 3: id of daughter 2 at dividing time; col 4: the id
                % of lineage
        dT2fT;
        tripolar_lineage_id; % a vector saves (lineage-) ids of tripolar lineages
    end
    
    methods
        
        function  getStatus(obj)
            disp(['Generating ' num2str(size(obj.fullTracks, 2)) ' full tracks, ' num2str(size(obj.partialTracks, 2)) ' partial tracks, '...
                num2str(size(obj.deadTracks, 2)) ' dead tracks.']);
        end
        
        function fullTracks = getFullCellTracks(obj, startFrame, endFrame, GMMs)
            [gmms] = GMMs{1, endFrame};
            listSize = length(gmms); % size of gmms at the last frame ..
            fullTracks = cell(1, listSize);
            for i = 1:listSize
%                 disp(['Scanning ' num2str(i) '/' num2str(listSize) ' fTracks.'])
                fullTracks{1, i} = cellTrack;
                fullTracks{1, i} = fullTracks{1, i}.initializeCellTrack(startFrame, endFrame);
                blob = gmms(i);
                if isempty(blob.m)
                    continue;
                end
                x = blob.m(1); y = blob.m(2);
                idx = blob.lineage + 1;
                t = endFrame; svidx = blob.svIdx + 1;
                fullTracks{1, i} = fullTracks{1, i}.write2track(t, idx, x, y, svidx);
                parent = blob.parent + 1;
                if parent > 0
                    fullTracks{1, i}.indices(t-1) = parent;
                end
            end
            
            for frame = endFrame-1:-1:startFrame
                [gmms] = GMMs{1, frame};
                for j = 1:listSize
                    k = fullTracks{1, j}.indices(frame);
                    if k <= 0 %|| k > size(gmms, 1)                       
                        continue;
                    end
%                     disp(['The ' num2str(j) '-th row']);
%                     disp([', its parent id is ' num2str(fk)]);
%                     disp([', at ' num2str(k)]);
%                     disp(['-th row from frame ' num2str(frame)]);

                    blob = gmms(k);
                    if isempty(blob.m)
                        continue;
                    end
                    x = blob.m(1); y = blob.m(2);
                    idx = blob.lineage + 1;
                    t = frame; svidx = blob.svIdx + 1;
                    fullTracks{1, j} = fullTracks{1, j}.write2track(t, k, x, y, svidx);
                    parent = blob.parent + 1;
                    if parent > 0 && t > 1
                        fullTracks{1, j}.indices(t-1) = parent;
                    end
                end
            end
        end
        
        function partialTracks = getPartialCellTracks(obj, startFrame, endFrame)
            partialTracks = [];
            counter = 0;
            for i = 1:size(obj.fullTracks, 2)
%                 disp(['Scanning ' num2str(i) '/' num2str(size(obj.fullTracks, 2)) ' pTracks.'])
                if prod(obj.fullTracks{1, i}.indices) == 0
                    counter = counter + 1;
                    partialTracks{1, counter} = obj.fullTracks{1, i};
                end
            end
        end
        
        function [idxMatrix, coordMatrix2X, coordMatrix2Y] = getCoordMatrix(obj, startFrame, endFrame)
            listSize = size(obj.fullTracks, 2);
            idxMatrix = zeros(listSize, endFrame-startFrame+1);
            coordMatrix2X = zeros(listSize, endFrame-startFrame+1);
            coordMatrix2Y = zeros(listSize, endFrame-startFrame+1);
            for i = 1:listSize
                idxMatrix(i,:) = obj.fullTracks{1,i}.indices;
                coordMatrix2X(i,:) = obj.fullTracks{1,i}.xCoord;
                coordMatrix2Y(i,:) = obj.fullTracks{1,i}.yCoord;
            end
        end
        
        function [idxMatrix, coordMatrix2X, coordMatrix2Y, SVidxMatrix] = getMatrix(obj, startFrame, endFrame)
            listSize = size(obj.fullTracks, 2);
            idxMatrix = zeros(listSize, endFrame-startFrame+1);
            coordMatrix2X = zeros(listSize, endFrame-startFrame+1);
            coordMatrix2Y = zeros(listSize, endFrame-startFrame+1);
            SVidxMatrix = cell(listSize, endFrame-startFrame+1);
            for i = 1:listSize                
                idxMatrix(i,:) = obj.fullTracks{1,i}.indices;
                coordMatrix2X(i,:) = obj.fullTracks{1,i}.xCoord;
                coordMatrix2Y(i,:) = obj.fullTracks{1,i}.yCoord;
                for j = startFrame:endFrame
                    SVidxMatrix{i, j} = obj.fullTracks{1,i}.svIdx{1, j};
                end
            end
        end
        
        function [idxMatrix, coordMatrix2X, coordMatrix2Y, SVidxMatrix] = get_FT_DT_Matrix(obj, startFrame, endFrame)
            listSize = size(obj.fullTracks, 2); listSize2 = size(obj.deadTracks, 2);
            idxMatrix = zeros(listSize+listSize2, endFrame-startFrame+1);
            coordMatrix2X = zeros(listSize+listSize2, endFrame-startFrame+1);
            coordMatrix2Y = zeros(listSize+listSize2, endFrame-startFrame+1);
            SVidxMatrix = cell(listSize+listSize2, endFrame-startFrame+1);
            for i = 1:listSize                
                idxMatrix(i,:) = obj.fullTracks{1,i}.indices;
                coordMatrix2X(i,:) = obj.fullTracks{1,i}.xCoord;
                coordMatrix2Y(i,:) = obj.fullTracks{1,i}.yCoord;
%                 for j = startFrame:endFrame
                    SVidxMatrix(i, :) = obj.fullTracks{1,i}.svIdx;
%                 end
            end
            for k = 1:listSize2
                idxMatrix(listSize+k,:) = obj.deadTracks{1,k}.indices;
                coordMatrix2X(listSize+k,:) = obj.deadTracks{1,k}.xCoord;
                coordMatrix2Y(listSize+k,:) = obj.deadTracks{1,k}.yCoord;
%                 for j = startFrame:endFrame
                    SVidxMatrix(listSize+k, :) = obj.deadTracks{1,k}.svIdx;
%                 end
            end
        end
        
        
        function registeredDeadEvents = getRegisteredDeadEvents(obj, t) % registering all dead cells at time t, by recording cell idx
            registeredDeadEvents = [];
            if ~isempty(obj.deadTracks)
                for i = 1:size(obj.deadTracks, 2)
                    obj.deadTracks{1, i}.indices(t)
                    registeredDeadEvents = [registeredDeadEvents; obj.deadTracks{1, i}.indices(t)];
                end
%             else
%                  disp('deadTrack struct is empty ...');
            end
        end
        
        function deadTracks = identifyDeadCellTracks(obj, startFrame, endFrame, GMMs)
            deadTracks = [];
            counter = 0;
            for t = endFrame:-1:startFrame+1
                disp(['Scanning frame ... ' num2str(t-1)]);
                gmms1 = GMMs{1, t}; 
                gmms2 = GMMs{1, t-1};
                
                if size(gmms1, 1) == 1
                    gmms1 = gmms1';
                end
                if size(gmms2, 1) == 1
                    gmms2 = gmms2';
                end
                
                pIdx = [];                
%                 disp(['size(gmms1, 1) = ' num2str(size(gmms1, 1)) ', size(gmms2, 1) = ' num2str(size(gmms2, 1))]);
                for i = 1:size(gmms1, 1)
                    if ~isempty(gmms1(i).W)
                        pIdx = [pIdx gmms1(i).parent]; % parent idx of cells at t
                    end
                end
                cIdx = [];
                for j = 1:size(gmms2, 1)
                    cIdx = [cIdx gmms2(j).lineage]; % cell idx of cells at t-1
                end

                deadOnes = setdiff(cIdx, pIdx);                
                if isempty(deadOnes)
%                     disp(['No dead tracks found ... at time t = ' num2str(t) ' (time counting from 1).']);
                    continue;
                end
                
                for k = 1:length(deadOnes)
                    dIdx = deadOnes(k)+1; % dead cell idx at time t-1, counting from 1                 
                    
                    registeredDeadEvents = obj.getRegisteredDeadEvents(t-1);
                    
                    if ismember(dIdx, registeredDeadEvents)
                        continue;
                    end
                    
                    counter = counter + 1;
                    deadTracks{1, counter} = cellTrack;
                    deadTracks{1, counter} = deadTracks{1, counter}.initializeDeadCellTrack(startFrame, endFrame);
                    x = gmms2(dIdx).m(1); 
                    y = gmms2(dIdx).m(2); 
                    svidx = gmms2(dIdx).svIdx+1;
                    deadTracks{1, counter} = deadTracks{1, counter}.write2track(t-1, dIdx, x, y, svidx); % a track ends at t-1, so at time t-1, we need to write in its own idx
%                     disp(['Cell at x: ' num2str(num2str(x)) ' y: ' num2str(y) ' disappears at frame (mamut) ' num2str(t-2) ' of id = ' num2str(deadOnes(k))]);
                    idx2add = deadTracks{1, counter}.search4parent(t-1, dIdx, GMMs);
                    
                    tao = t-2;
                    for id = 1:length(idx2add)
%                         disp(['Scanning ' num2str(id) '/' num2str(length(idx2add)) ' dTracks.'])
                        currentId = idx2add(id);
                        deadTracks{1, counter}.framesPresented(tao) = tao;
                        deadTracks{1, counter}.indices(tao) = GMMs{1, tao}(currentId).id + 1; % here might be wrong ...
                        deadTracks{1, counter}.xCoord(tao) = GMMs{1, tao}(currentId).m(1);
                        deadTracks{1, counter}.yCoord(tao) = GMMs{1, tao}(currentId).m(2);
                        deadTracks{1, counter}.svIdx{1, tao} = GMMs{1, tao}(currentId).svIdx + 1;
                        tao = tao - 1;
                    end
                end
            end
        end
        
        function [deadTracks, cancelOff] = canceloffDTracksWithFTracks(obj, startFrame, endFrame)
            fTracks = obj.fullTracks; dTracks = obj.deadTracks;
            FTidxMatrix = zeros(size(fTracks, 2), endFrame - startFrame + 1);            
            parfor i = 1:size(fTracks, 2)
                FTidxMatrix(i,:) = fTracks{1, i}.indices;
            end            
            cancelOff = [];
            parfor i = 1:size(dTracks, 2)
                i
                c = (ones(size(fTracks, 2), 1)*dTracks{1, i}.indices).*FTidxMatrix;
                result = find(sum(c') > 0.9*(endFrame - startFrame + 1));
                if ~isempty(result)
                    cancelOff = [cancelOff result];
                end
            end
            cancelOff = fliplr(unique(cancelOff));
            if ~isempty(cancelOff) 
                disp([num2str(length(cancelOff)) ' dead tracks are cancelled off ..' ]);
                for j = 1:length(cancelOff)
                    obj.deadTracks(:, cancelOff(j)) = [];
                end               
            end
            deadTracks = obj.deadTracks;
            disp('Done');
        end
        
        
        function e = isEmptyTracks(obj, pPool)
            e = 'true';  a = [];
            pPoolSize = size(pPool, 2);
            for i = 1:pPoolSize
                a = [a; pPool{1, i}.isNotEmptyTrack()];
            end
            if sum(a) > 0
                e = 'false';
            end            
            pos = find(a > 0); 
%             str = '';
%             for i = 1:length(pos)
%                 str = [str ' ' num2str(pos(i)) '-th'];
%             end
%             disp(['The partial tracks remain to be wired are the ' str ' patial tracks.']);           
        end
        
        function remain = deadTracksRemain(obj)
            tracks = obj.deadTracks;
            remain = 'false';
            for i = 1:size(tracks, 2)
                if sum(tracks{1, i}.indices) > 0
                    remain = 'true';
                    disp(['The dTracks remain to be checked are ' num2str(i) '-th dTrack.']);
%                     break;
                end
            end
        end
        
        function [obj, numRescued] = rescue2(obj, slidingWindow, distance, maxLoop)
            pPool = obj.partialTracks; pPoolSize = size(pPool, 2);
            dPool = obj.deadTracks; dPoolSize = size(dPool, 2);
            numRescued = 0;
            if pPoolSize > 0 && dPoolSize > 0                 
                numLoop = 0;
                numRescued = zeros(pPoolSize, 1);
                
                pPoolCoordMatrix2X = zeros(pPoolSize, length(pPool{1, 1}.indices));
                pPoolCoordMatrix2Y = zeros(pPoolSize, length(pPool{1, 1}.indices));
                dPoolCoordMatrix2X = zeros(dPoolSize, length(dPool{1, 1}.indices));
                dPoolCoordMatrix2Y = zeros(dPoolSize, length(dPool{1, 1}.indices));
                vanish_dx = zeros(dPoolSize, 1); vanish_dy = zeros(dPoolSize, 1);
                
                for j = 1:pPoolSize
                    pPoolCoordMatrix2X(j, :) = pPool{1, j}.xCoord;
                    pPoolCoordMatrix2Y(j, :) = pPool{1, j}.yCoord;
                end
                for j = 1:dPoolSize
                    dPoolCoordMatrix2X(j, :) = dPool{1, j}.xCoord;
                    dPoolCoordMatrix2Y(j, :) = dPool{1, j}.yCoord;
                    dead_t = max(find(dPoolCoordMatrix2X(j, :)>0));
                    if isempty(dead_t)
                        dead_t = max(find(dPoolCoordMatrix2Y(j, :)>0));
                    end
                    vanish_dx(j, 1) = dPool{1, j}.xCoord(dead_t);
                    vanish_dy(j, 1) = dPool{1, j}.yCoord(dead_t);
                end
                %             checkStas = isequal(obj.isEmptyTracks(pPool), 'false');
                while numLoop < maxLoop
                    numLoop = numLoop + 1;
                    disp([' Rescuing Loop ... ' num2str(numLoop)]);
                    for p = 1:pPoolSize
                        disp(['Rescuing Loop ... ' num2str(numLoop) '   working on ' num2str(p) '/' num2str(pPoolSize) ' pTracks ...'])
                        px = ones(dPoolSize, 1) * pPoolCoordMatrix2X(p, :);
                        py = ones(dPoolSize, 1) * pPoolCoordMatrix2Y(p, :);
                        t = min(find(px(1,:) > 0)); % the time point when the pTrack first pops up
                        sliding_t = max(1, t-1-slidingWindow):t-1; %1:max(1, t-1-slidingWindow)
                        tmp_x = px.*dPoolCoordMatrix2X > 0; 
%                         aaa = ones(size(tmp_x, 1), 1); bbb = [1:1:size(tmp_x, 2)]; ccc = aaa*bbb;
%                         tmp_x = tmp_x.*ccc;                  
                        checker1 = sum(tmp_x')'; 
                        checker2 = sum(dPoolCoordMatrix2X(:, sliding_t)')';
                        if length(sliding_t) == 1
                            checker2 = dPoolCoordMatrix2X(:, sliding_t);
                        end

                        pos1 = unique([find(checker1 == 0); find(tmp_x(:, (t)) > 0)]); pos2 = find(checker2 > 0);
                        pos = intersect(pos1, pos2);
                        if isempty(pos)
%                                                     disp('isempty(pos)')
                            continue;
                        end
                        px = px(1,:); py = py(1,:);
                        myX = px(t);                        myY = py(t); % coordinate of the partial track, when it pops up firstly...
                        myXnext = px(min(t+1, length(px)));                        myYnext = py(min(t+1, length(py)));
                        myX2 = ones(length(vanish_dx), 1)*myX;                        myY2 = ones(length(vanish_dy), 1)*myY;
                        dist = sqrt((myX2 - vanish_dx).^2 + (myY2 - vanish_dy).^2);
                        candidate = [];
                        pos3 = intersect(find(dist <= distance), pos);
                        
                        for po = 1:length(pos3)
                            P0 = [myXnext myYnext]; P1 = [myX myY]; P2 = [vanish_dx(pos3(po)) vanish_dx(pos3(po))];
                            ang = atan2(abs(det([P2-P0;P1-P0])),dot(P2-P0,P1-P0))*180/pi;
                            candidate = [candidate; vanish_dx(pos3(po)) vanish_dy(pos3(po)) pos3(po) ang dist(pos3(po))];
                        end
                        if isempty(candidate)
%                                                     disp('Empty candidate');
                        end
                        
                        if size(candidate, 1) > 0
                            candidate = sortrows(candidate, 4, 'ascend');
                            toChoose = candidate(1, 3);
                            disp(['Loop ' num2str(numLoop) '... connecting ... the ' num2str(p) '-th pTrack with ' num2str(toChoose) '-th dTrack!']);
                            %                         disp('before');
                            %                         pPoolCoordMatrix2X(p, :)
                            pPool{1, p} = pPool{1, p}.updatePartialTrack(dPool{1, toChoose});
                            obj.partialTracks{1, p} = obj.partialTracks{1, p}.updatePartialTrack(dPool{1, toChoose});
                            %                         dPool{1, toChoose} = dPool{1, toChoose}.deleteTrack();
                            obj.deadTracks{1, toChoose} = obj.deadTracks{1, toChoose}.deleteTrack();
                            %                         dPoolCoordMatrix2X(toChoose, :) = zeros(1, length(dPoolCoordMatrix2X(toChoose, :)));
                            %                         dPoolCoordMatrix2Y(toChoose, :) = zeros(1, length(dPoolCoordMatrix2Y(toChoose, :)));
                            pPoolCoordMatrix2X(p, :) = pPool{1, p}.xCoord;
                            pPoolCoordMatrix2Y(p, :) = pPool{1, p}.yCoord;
                            %                         disp('after');
                            %                         pPoolCoordMatrix2X(p, :)
                            %                         progress = [progress; p toChoose];
                        end
                        
                        ftc = pPool{1, p}.isFullyTracked();
                        if ftc == 1
                            obj = obj.addingPT2FTs(pPool{1, p}, pPool{1, p}.indices(end)); % adding corrected pTracks to obj.fTracks
                            numRescued(p, 1) = 1;
                            %                     pPool{1, j} = pPool{1, j}.deleteTrack();
                            %                         disp(['Loop ' num2str(numLoop) ':This pTrack is now fully tracked = ' num2str(ftc)]);
                        end
                    end
                    
                end
                disp([num2str(sum(numRescued)) ' incomplete tracks are rescued!']);
                
            end
        end
        
        function [obj, numRescued] = rescueFT(obj, slidingWindow, distance, maxLoop)
            pPool = obj.fullTracks; pPoolSize = size(pPool, 2);  
            dPool = obj.deadTracks; dPoolSize = size(dPool, 2);          
            numLoop = 0;         
            numRescued = zeros(pPoolSize, 1); 
%             progress = [];
            pPoolCoordMatrix2X = zeros(pPoolSize, length(pPool{1, 1}.indices));
            pPoolCoordMatrix2Y = zeros(pPoolSize, length(pPool{1, 1}.indices));
            dPoolCoordMatrix2X = zeros(dPoolSize, length(dPool{1, 1}.indices));
            dPoolCoordMatrix2Y = zeros(dPoolSize, length(dPool{1, 1}.indices));
            vanish_dx = zeros(dPoolSize, 1); vanish_dy = zeros(dPoolSize, 1);            
            for j = 1:pPoolSize
                pPoolCoordMatrix2X(j, :) = pPool{1, j}.xCoord;
                pPoolCoordMatrix2Y(j, :) = pPool{1, j}.yCoord;
            end
            for j = 1:dPoolSize
                dPoolCoordMatrix2X(j, :) = dPool{1, j}.xCoord;
                dPoolCoordMatrix2Y(j, :) = dPool{1, j}.yCoord;
                dead_t = max(find(dPoolCoordMatrix2X(j, :)>0));
                vanish_dx(j, 1) = dPool{1, j}.xCoord(dead_t); 
                vanish_dy(j, 1) = dPool{1, j}.yCoord(dead_t); 
            end
%             checkStas = isequal(obj.isEmptyTracks(pPool), 'false');      
            while numLoop < maxLoop
                numLoop = numLoop + 1;   
                disp([' Rescuing Loop ... ' num2str(numLoop)]);
                for p = 1:pPoolSize
                    disp(['   working on ' num2str(p) '/' num2str(pPoolSize) ' pTracks ...'])
                    px = ones(dPoolSize, 1) * pPoolCoordMatrix2X(p, :);
                    py = ones(dPoolSize, 1) * pPoolCoordMatrix2Y(p, :);
                    t = min(find(px(1,:) > 0)); % the time point when the pTrack first pops up
                    sliding_t = max(1, t-1-slidingWindow):t-1; %1:max(1, t-1-slidingWindow)
                    tmp_x = px.*dPoolCoordMatrix2X;
                    checker1 = sum(tmp_x')'; checker2 = sum(dPoolCoordMatrix2X(:, sliding_t)')';
                    if length(sliding_t) == 1
                        checker2 = dPoolCoordMatrix2X(:, sliding_t);
                    end
                    pos1 = find(checker1 == 0); pos2 = find(checker2 > 0);
                    pos = intersect(pos1, pos2);
                    if isempty(pos)
                        continue;
                    end
                    px = px(1,:); py = py(1,:);
                    myX = px(t); myY = py(t); % coordinate of the partial track, when it pops up firstly...
                    myXnext = px(min(t+1, length(px))); myYnext = py(min(t+1, length(py)));
                    myX2 = ones(length(vanish_dx), 1)*myX; myY2 = ones(length(vanish_dy), 1)*myY;
                    dist = sqrt((myX2 - vanish_dx).^2 + (myY2 - vanish_dy).^2);
                    candidate = [];
                    pos3 = intersect(find(dist <= distance), pos);
                    for po = 1:length(pos3)
                        P0 = [myXnext myYnext]; P1 = [myX myY]; P2 = [vanish_dx(pos3(po)) vanish_dy(pos3(po))];
                        ang = atan2(abs(det([P2-P0;P1-P0])),dot(P2-P0,P1-P0))*180/pi;
                        candidate = [candidate; vanish_dx(pos3(po)) vanish_dy(pos3(po)) pos3(po) ang dist(pos3(po))];
%                         disp(['P0 = ' num2str(P0) ', P1 = ' num2str(P1) ', P2 = ' num2str(P2)]);
                    end     
                    
%                     candidate
                    if size(candidate, 1) > 0
                        candidate = sortrows(candidate, 4, 'ascend');
                        toChoose = candidate(1, 3);
                        disp(['Loop ' num2str(numLoop) '... connecting ... the ' num2str(p) '-th pTrack with ' num2str(toChoose) '-th dTrack!']);
%                         disp('before'); 
%                         pPoolCoordMatrix2X(p, :)
                        pPool{1, p} = pPool{1, p}.updatePartialTrack(dPool{1, toChoose});
                        obj.fullTracks{1, p} = obj.fullTracks{1, p}.updatePartialTrack(dPool{1, toChoose});
                        dPool{1, toChoose} = dPool{1, toChoose}.deleteTrack();
                        obj.deadTracks{1, toChoose} = obj.deadTracks{1, toChoose}.deleteTrack();
                        dPoolCoordMatrix2X(toChoose, :) = zeros(1, length(dPoolCoordMatrix2X(toChoose, :)));
                        dPoolCoordMatrix2Y(toChoose, :) = zeros(1, length(dPoolCoordMatrix2Y(toChoose, :)));
                        pPoolCoordMatrix2X(p, :) = pPool{1, p}.xCoord;
                        pPoolCoordMatrix2Y(p, :) = pPool{1, p}.yCoord;
%                         disp('after');
%                         pPoolCoordMatrix2X(p, :)
%                         progress = [progress; p toChoose]; 
                    end
                    
                    ftc = pPool{1, p}.isFullyTracked();
                    if ftc == 1
                        obj = obj.addingPT2FTs(pPool{1, p}, pPool{1, p}.indices(end)); % adding corrected pTracks to obj.fTracks
                        numRescued(p, 1) = 1;
                        %                     pPool{1, j} = pPool{1, j}.deleteTrack();
                        %                         disp(['Loop ' num2str(numLoop) ':This pTrack is now fully tracked = ' num2str(ftc)]);
                    end                    
                end       
                
            end
            disp([num2str(sum(numRescued)) ' incomplete tracks are rescued!']);
        end
        
         function [obj, numRescued] = rescue3(obj, slidingWindow, distance, maxLoop)
            pPool = obj.partialTracks; pPoolSize = size(pPool, 2);
            dPool = obj.deadTracks; dPoolSize = size(dPool, 2);
            numLoop = 0;
            numRescued = zeros(pPoolSize, 1);
            %             progress = [];
            pPoolCoordMatrix2X = zeros(pPoolSize, length(pPool{1, 1}.indices));
            pPoolCoordMatrix2Y = zeros(pPoolSize, length(pPool{1, 1}.indices));
            dPoolCoordMatrix2X = zeros(dPoolSize, length(dPool{1, 1}.indices));
            dPoolCoordMatrix2Y = zeros(dPoolSize, length(dPool{1, 1}.indices));
            vanish_dx = zeros(dPoolSize, 1); vanish_dy = zeros(dPoolSize, 1);
            for j = 1:pPoolSize
                pPoolCoordMatrix2X(j, :) = pPool{1, j}.xCoord;
                pPoolCoordMatrix2Y(j, :) = pPool{1, j}.yCoord;
            end
            for j = 1:dPoolSize
                dPoolCoordMatrix2X(j, :) = dPool{1, j}.xCoord;
                dPoolCoordMatrix2Y(j, :) = dPool{1, j}.yCoord;
                dead_t = max(find(dPoolCoordMatrix2X(j, :)>0));
                if isempty(dead_t)
                    continue;
                end
                vanish_dx(j, 1) = dPool{1, j}.xCoord(dead_t);
                vanish_dy(j, 1) = dPool{1, j}.yCoord(dead_t);
            end
            %             checkStas = isequal(obj.isEmptyTracks(pPool), 'false');
                        
            updateTB = zeros(pPoolSize, maxLoop);
            while numLoop < maxLoop
                numLoop = numLoop + 1;
                disp([' 3 Rescuing Loop ... ' num2str(numLoop)]);   
                parfor p = 1:pPoolSize
                    tmpx = pPoolCoordMatrix2X; tmpy = pPoolCoordMatrix2Y;
                    px = ones(dPoolSize, 1) * tmpx(p, :);
                    py = ones(dPoolSize, 1) * tmpy(p, :);
                    t = min(find(px(1,:) > 0)); % the time point when the pTrack first pops up
                    sliding_t = max(1, t-1-slidingWindow):t-1; %1:max(1, t-1-slidingWindow)
                    tmp_x = px.*dPoolCoordMatrix2X;
                    checker1 = sum(tmp_x')'; checker2 = sum(tmpx(:, sliding_t)')';
                    if length(sliding_t) == 1
                        checker2 = tmpx(:, sliding_t);
                    end
                    pos1 = find(checker1 == 0); pos2 = find(checker2 > 0);
                    pos = intersect(pos1, pos2);
                    if isempty(pos)
                        continue;
                    end
                    px = px(1,:); py = py(1,:);
                    myX = px(t); myY = py(t); % coordinate of the partial track, when it pops up firstly...
                    myXnext = px(min(t+1, length(px))); myYnext = py(min(t+1, length(py)));
                    myX2 = ones(length(vanish_dx), 1)*myX; myY2 = ones(length(vanish_dy), 1)*myY;
                    dist = sqrt((myX2 - vanish_dx).^2 + (myY2 - vanish_dy).^2);
                    candidate = [];
                    pos3 = intersect(find(dist <= distance), pos);
%                     disp(['Length candidate: ' num2str(length(pos3))]);
                    for po = 1:length(pos3)
                        P0 = [myXnext myYnext]; P1 = [myX myY]; P2 = [myX2(pos3(po)) myY2(pos3(po))];
                        ang = atan2(abs(det([P2-P0;P1-P0])),dot(P2-P0,P1-P0))*180/pi;
                        candidate = [candidate; vanish_dx(pos3(po)) vanish_dy(pos3(po)) pos3(po) ang dist(pos3(po))];
                    end
                    if isempty(candidate)
                        continue;
                    end
                    
                    candidate = sortrows(candidate, 4, 'ascend');
%                     toChoose = candidate(1, 3);
                    %                         disp(['Loop ' num2str(numLoop) '... connecting ... the ' num2str(p) '-th pTrack with ' num2str(toChoose) '-th dTrack!']);
                    %                         disp('before');
                    %                         pPoolCoordMatrix2X(p, :)
%                     tmpTB(p, :) = candidate(1, 3);
                    updateTB(p, numLoop) = candidate(1, 3);
                end
            end
            for p = 1:pPoolSize
                toChoose = updateTB(p, 1);
                if toChoose > 0
                    pPool{1, p} = pPool{1, p}.updatePartialTrack(dPool{1, toChoose});
                    obj.partialTracks{1, p} = obj.partialTracks{1, p}.updatePartialTrack(dPool{1, toChoose});
                    dPool{1, toChoose} = dPool{1, toChoose}.deleteTrack();
                    obj.deadTracks{1, toChoose} = obj.deadTracks{1, toChoose}.deleteTrack();
                    dPoolCoordMatrix2X(toChoose, :) = zeros(1, length(dPoolCoordMatrix2X(toChoose, :)));
                    dPoolCoordMatrix2Y(toChoose, :) = zeros(1, length(dPoolCoordMatrix2Y(toChoose, :)));
                    pPoolCoordMatrix2X(p, :) = pPool{1, p}.xCoord;
                    pPoolCoordMatrix2Y(p, :) = pPool{1, p}.yCoord;
                    %                         disp('after');
                    %                         pPoolCoordMatrix2X(p, :)
                    %                         progress = [progress; p toChoose];
                    
                    ftc = pPool{1, p}.isFullyTracked();
                    if ftc == 1
                        obj = obj.addingPT2FTs(pPool{1, p}, pPool{1, p}.indices(end)); % adding corrected pTracks to obj.fTracks
                        numRescued(p, 1) = 1;
                        %                     pPool{1, j} = pPool{1, j}.deleteTrack();
                        %                         disp(['Loop ' num2str(numLoop) ':This pTrack is now fully tracked = ' num2str(ftc)]);
                    end
                end
            end                
            end

        
        function [obj, numRescued] = rescue(obj, slidingWindow, distance, maxLoop)
            pPool = obj.partialTracks; 
            dPool = obj.deadTracks;
            pPoolSize = size(pPool, 2);            
            numLoop = 0;  
%             progress2  = [];            
            checkStas = isequal(obj.isEmptyTracks(pPool), 'false');      
            numRescued = 0;
            while checkStas && numLoop < maxLoop
                numLoop = numLoop + 1;                
                for j = 1:pPoolSize % for debugging, change here the index of pTracks to specific number ...
                    pTrack = pPool{1, j};                    
                    if pTrack.isNotEmptyTrack()==0
%                         disp(['Loop ' num2str(numLoop) ' the ' num2str(j) '-th pTrack is empty.']);
                        continue;
                    end
                    t = pTrack.getShowtime4partialTrack();
                    myX = pTrack.xCoord(t); myY = pTrack.yCoord(t); % coordinate of the partial track, when it pops up firstly...
                    myXnext = pTrack.xCoord(min(t+1, length(pTrack.xCoord))); myYnext = pTrack.yCoord(min(t+1, length(pTrack.yCoord)));
                    t = t - 1;
                    xys = [];
                  
                    for i = 1:size(dPool, 2)
                        if prod(dPool{1, i}.indices(t+1:min(size(dPool{1, i}.indices, 2), t+slidingWindow))) > 0
%                             disp(['Loop ' num2str(numLoop) ' the ' num2str(i) '-th dTrack is none 0.']);
                            continue;
                        end
                        if sum(dPool{1, i}.indices(max(1, t-slidingWindow):t)) == 0
%                             disp(['Loop ' num2str(numLoop) ' the ' num2str(i) '-th dTrack is empty.']);
                            continue;
                        end
                        for slidingt = t:-1:max(1, t-slidingWindow)
                            x = dPool{1, i}.xCoord(slidingt);
                            y = dPool{1, i}.yCoord(slidingt);
                            if x*y ~= 0
                                xys = [xys; x y i];               
                                break;
                            end
                        end
                        if isempty(xys)
                            disp(['Loop ' num2str(numLoop) ' the ' num2str(i) '-th dTrack is empty with coordinates.']);
                        end
                    end
                    candidate = [];
                    %/* adjusted on June 19
                    myX2 = ones(size(xys, 1), 1)*myX; myY2 = ones(size(xys, 1), 1)*myY;
                    dist = sqrt((myX2 - xys(:, 1)).^2 + (myY2 - xys(:, 2)).^2);
                    pos = find(dist <= distance);
                    for po = 1:length(pos)
                        P0 = [myXnext myYnext]; P1 = [myX myY]; P2 = [myX2(pos(po)) myY2(pos(po))];
                        ang = atan2(abs(det([P2-P0;P1-P0])),dot(P2-P0,P1-P0))*180/pi;
                        candidate = [candidate; xys(pos(po), :) ang dist(pos(po))];
                    end
                    %                     for k = 1:size(xys, 1)
%                         thisX = xys(k, 1); thisY = xys(k, 2);
%                         dist = sqrt((myX - thisX)^2 + (myY - thisY)^2);
% %                         disp(['XYS is not empty, dist = ' num2str(dist)]);
%                         if dist <= distance
%                             P0 = [myXnext myYnext]; P1 = [myX myY]; P2 = [thisX thisY]; 
%                             ang = atan2(abs(det([P2-P0;P1-P0])),dot(P2-P0,P1-P0))*180/pi;
%                             candidate = [candidate; xys(k, :) ang dist];
% %                             disp(['candidate(of x=' num2str(myX) ',y=' num2str(myY) ') info: x=' num2str(thisX) ' y=' num2str(thisY) ' ang=' num2str(ang) ' dist=' num2str(dist)])
% %                             disp(['Loop ' num2str(numLoop) '... connecting ... the ' num2str(j) '-th pTrack with ' num2str(xys(k, 3)) '-th dTrack!']);
% %                             pPool{1, j} = pPool{1, j}.updatePartialTrack(dPool{1, xys(k, 3)});
% %                             obj.partialTracks{1, j} = obj.partialTracks{1, j}.updatePartialTrack(dPool{1, xys(k, 3)});
% %                             dPool{1, xys(k, 3)} = dPool{1, xys(k, 3)}.deleteTrack();
% %                             obj.deadTracks{1, xys(k, 3)} = obj.deadTracks{1, xys(k, 3)}.deleteTrack();
% %                             break; % here we assume ONLY one cell will be close to ...
%                         end
%                     end
                          %*/
                    
                    if size(candidate, 1) > 0
                        % /* adjusted on June 19
                        candidate = sortrows(candidate, 4, 'ascend');
                        toChoose = candidate(1, 3);
%                         pos = find(candidate(:, 4) == min(candidate(:, 4))); % first find the minimul turning angle
%                         toChoose = candidate(pos(1), 3);
%                         for pp = 1:length(pos)
%                             if candidate(pos(pp), 5) < candidate(pos(1), 5) % if multiple minimal turning angles found, use distance
%                                 toChoose = candidate(pos(pp), 3);
%                             end
%                         end
                         % */
                        disp(['Loop ' num2str(numLoop) '... connecting ... the ' num2str(j) '-th pTrack with ' num2str(toChoose) '-th dTrack!']);
                        
                        pPool{1, j} = pPool{1, j}.updatePartialTrack(dPool{1, toChoose});
                        obj.partialTracks{1, j} = obj.partialTracks{1, j}.updatePartialTrack(dPool{1, toChoose});
                        dPool{1, toChoose} = dPool{1, toChoose}.deleteTrack();
                        obj.deadTracks{1, toChoose} = obj.deadTracks{1, toChoose}.deleteTrack();
%                         progress2 = [progress2; j toChoose];
                    end
                   
                    ftc = pPool{1, j}.isFullyTracked();
                    if ftc == 1
                        obj = obj.addingPT2FTs(pPool{1, j}, pPool{1, j}.indices(end)); % adding corrected pTracks to obj.fTracks
                        numRescued = numRescued + 1;
                        pPool{1, j} = pPool{1, j}.deleteTrack();
%                         disp(['Loop ' num2str(numLoop) ':This pTrack is now fully tracked = ' num2str(ftc)]);
                    end
                end
%                 disp(['After LOOP ... ' num2str(numLoop)]);
                checkStas = isequal(obj.isEmptyTracks(pPool), 'false');
%                 disp(['Loop ' num2str(numLoop) ': obj.isEmptyTrackS(pPool) = ' num2str(checkStas)]);                
            end     
            disp([num2str(numRescued) ' cell tracks are rescued!']);
        end % examines if a pTrack could be rescued
        
                
        function obj = addingPT2FTs(obj, pT, i) % adding a corrected pTrack to the list of full track, or to the i-th row of full tracks
            if i <= size(obj.fullTracks, 2)
                obj.fullTracks{1, i} = pT;
            else
            obj.fullTracks{1, size(obj.fullTracks, 2) + 1} = pT;
            end
% %             obj.getStatus();
        end
        
        
        function [ alpharad ] = anglevec(obj, veca, vecb )  % Calculate angle between two vectors
            alpharad = acos(dot(veca, vecb) / sqrt( dot(veca, veca) * dot(vecb, vecb)));
        end
    
        function obj = identifyLineage(obj, startFrame, endFrame)
            obj.tracks = [obj.fullTracks obj.deadTracks obj.partialTracks];
            listSize = size( obj.tracks, 2);
            idxMatrix = zeros(listSize, endFrame-startFrame+1);           
            for i = 1:listSize
%                 disp(['i=' num2str(i) '/' num2str(listSize) '  ' num2str(length(obj.tracks{1,i}.indices))]);
                idxMatrix(i,:) =  obj.tracks{1,i}.indices;                
            end
            colCurrentFrame = idxMatrix(:, 1);
            [a,b] = hist(colCurrentFrame,unique(colCurrentFrame));
            a = a(2:end); b = b(2:end); % show a, b'
            order = find(a >= 2 );
            obj.lineage = cell(1, length(order)); 
            obj.mapLineage2Tracks = cell(1, length(order));
            for i = 1:length(order)
                idRepetive = b(order(i));
                pos = find(colCurrentFrame == idRepetive);
                flag = 0; % check if tracks within one lineage are equal
                for j = 1:length(pos) - 1
                    if isequal( obj.tracks{1, pos(j)},  obj.tracks{1, pos(j + 1)})
                        flag = 1; disp(['Lineage ' num2str(i) ' is found equal']);
                        break;                        
                    end
                end
                if flag == 0
                    tt = []; mapping = [];
                    for j = 1:length(pos)
                        tt = [tt obj.tracks{1, pos(j)}];
                        mapping = [mapping pos(j)];
                    end
                    obj.lineage{1, i} = tt;
                    obj.mapLineage2Tracks{1, i} = mapping;
%                     obj.lineage{1, i} = cell(1, length(pos));
%                     for j = 1:length(pos)
%                      obj.lineage{1, i}{1, j}  = tracks{1, pos(j)};%% pos(j) is added as index of fullTracks, easy to retrieve ...
%                     end
                end
            end
        end
        
        function obj = identifyLineageSimple(obj, startFrame, endFrame, includeDT)
            fTracks = obj.fullTracks;
            pTracks = obj.partialTracks;
            dTracks = obj.deadTracks;
            
            tracks = fTracks;

            if includeDT == 1
                tracks = [fTracks dTracks];
            end
            listSize = size(tracks, 2);
            listSize2 = size(fTracks, 2);
            
            idxMatrix = zeros(listSize, endFrame-startFrame+1);
            
            for i = 1:listSize
                %                 disp(['i=' num2str(i) '/' num2str(listSize) '  ' num2str(length(obj.tracks{1,i}.indices))]);
                idxMatrix(i,:) =  tracks{1,i}.indices;
            end
            
            colCurrentFrame = idxMatrix(:, 1); % WAS 1, LOOKING AT THE FIRST COLUMN
            [a,b] = hist(colCurrentFrame,unique(colCurrentFrame));
            a = a(2:end); b = b(2:end); % show a, b'
            order = find(a >= 2 );
            lineageCounter = 0;
            totalLineages = [];
            
            for k = 1:length(order)
                idRepetive = b(order(k));
                pos = find(colCurrentFrame == idRepetive);
                idMatrixOfLineage = []; pos2tracks = [];
                for p = 1:length(pos)
                    idMatrixOfLineage = [idMatrixOfLineage; idxMatrix(pos(p), :)];
                    pos2tracks = [pos2tracks; pos(p)];
                end
                
                divCounter = 0;
                lineage_tb = cell(1, 4);
                disp(['Scanning ' num2str(k) '-th out of ' num2str(length(order)) ' lineages.']);
                
                % To determine div time
                idxMinus = abs(idMatrixOfLineage(1, :) - idMatrixOfLineage(2, :));
                divTime = min(find(idxMinus > 0));
                % To get parent Id
                pId = idMatrixOfLineage(1, (divTime - 1));
                tmp = []; dist2parent = [];
                px = fTracks{1, pos(1)}.xCoord(divTime - 1);
                py = fTracks{1, pos(1)}.yCoord(divTime - 1);
                
                for p = 1:length(pos)
                    tmp = [tmp idMatrixOfLineage(p, (divTime))];
                    x = fTracks{1, pos(p)}.xCoord(divTime);
                    y = fTracks{1, pos(p)}.yCoord(divTime);
                    dist2parent= [dist2parent sqrt((px-x)^2 + (py-y)^2)];
                end
                
                divCounter = divCounter + 1;
                lineage_tb{divCounter, 1} = idRepetive;
                lineage_tb{divCounter, 2} = divTime;
                lineage_tb{divCounter, 3} = [pId, tmp]; % pId is the id of cell at frame divTime-1, tmp saves the ids of daughters at frame f+1
                
                dp = find(dist2parent == min(dist2parent));
                lineage_tb{divCounter, 4} = [pos(dp(1)), pos'];
                lineage_tb{divCounter, 5} = [px py];
                
                if ~isempty(lineage_tb)
                    lineageCounter = lineageCounter + 1;
                    totalLineages{1, lineageCounter} = lineage_tb;
                end                
            end
            obj.lineages = totalLineages;
            
            coordDivision = [];
            counter = 0;
            for p = 1:size(totalLineages, 2)
                theLineage = totalLineages{1, p};
                coordDivision = [coordDivision; p theLineage{5} theLineage{2}];
            end
            obj.lineagesTB = coordDivision;
        end
        
        function [obj, a, b, colCurrentFrame] = identifyLineage2(obj, startFrame, endFrame, includeDT) %colCurrentFrame
            obj.tracks = [obj.fullTracks];
            if includeDT == 1
                obj.tracks = [obj.fullTracks obj.deadTracks];
            end
            listSize = size( obj.tracks, 2);
            listSize2 = size( [obj.fullTracks ], 2);
          
            idxMatrix = zeros(listSize, endFrame-startFrame+1);
            
            for i = 1:listSize
%                 disp(['i=' num2str(i) '/' num2str(listSize) '  ' num2str(length(obj.tracks{1,i}.indices))]);
                idxMatrix(i,:) =  obj.tracks{1,i}.indices;
            end
%             idxMatrix
            colCurrentFrame = idxMatrix(:, 2); % WAS 1, LOOKING AT THE FIRST COLUMN
            [a,b] = hist(colCurrentFrame,unique(colCurrentFrame));
%             a = a(2:end); b = b(2:end); % show a, b'
            order = find(a >= 2 );
            lineageCounter = 0;
            totalLineages = [];
            
            for k = 2:length(order)               
                idRepetive = b(order(k));
                pos = find(colCurrentFrame == idRepetive);
                idMatrixOfLineage = []; pos2tracks = [];
                for p = 1:length(pos)
                    idMatrixOfLineage = [idMatrixOfLineage; idxMatrix(pos(p), :)];
                    pos2tracks = [pos2tracks; pos(p)];
                end
%                 idMatrixOfLineage
                divCounter = 0;
                lineage_tb = cell(1, 4);
                disp(['Scanning ' num2str(k) '-th out of ' num2str(length(order)) ' lineages.']);
                
                for f = size(idMatrixOfLineage, 2)-1:-1:2
                    if length(idMatrixOfLineage(:,f))-length(unique(idMatrixOfLineage(:,f))) == 0
                        %                         disp('length(idMatrixOfLineage(:,f))-length(unique(idMatrixOfLineage(:,f))) == 0')
                        continue;
                    else
                        if length(idMatrixOfLineage(:,f-1))-length(unique(idMatrixOfLineage(:,f-1))) == 0
                            continue;
                        else                        
                        [aa, bb] = hist(idMatrixOfLineage(:,f), unique(idMatrixOfLineage(:,f)));
                        div = find(aa>=2);
                        %                         disp(['f = ' num2str(f) ' ' num2str(length(div))]);
                        %                         div
                        for d = 1:length(div)
                            pId = bb(div(d));
                            
                            pos1 = find(idMatrixOfLineage(:,f) == pId);
                            if length(pos1) < 2 % determining number of daughters, 1 means no div
                                %                                 disp('length(pos1)');
                                continue;
                            end
                            tmp = zeros(1, length(pos1));
                            tmp1 = tmp;
                            for q = 1:length(pos1)
                                tmp(q) = idMatrixOfLineage(pos1(q), f+1); % check next frame if id is consistant ...
                                tmp1(q) = idMatrixOfLineage(pos1(q), end); % the idx at the last frame will be used as final cell Ids ..
                            end
                            if length(unique(tmp)) == 1
                                %                                 disp('length(unique(tmp)) == 1')
                                continue;
                            end
                            dist2parent = zeros(size(tmp1));
                            lifespan = 100;
                            for q = 1:length(pos1)
                                if tmp1(q) <= 0
                                    tmp1(q) = -(pos2tracks(pos1(q)) - listSize2);
                                    px = obj.deadTracks{1, -tmp1(q)}.xCoord(f); py = obj.deadTracks{1, -tmp1(q)}.yCoord(f);
                                    x = obj.deadTracks{1, -tmp1(q)}.xCoord(f + 1); y = obj.deadTracks{1, -tmp1(q)}.yCoord(f+1);
                                    lifespan = max(find(obj.deadTracks{1, -tmp1(q)}.indices > 0)) - f;
                                else
                                    px = obj.fullTracks{1, tmp1(q)}.xCoord(f); py = obj.fullTracks{1, tmp1(q)}.yCoord(f);
                                    x = obj.fullTracks{1, tmp1(q)}.xCoord(f + 1); y = obj.fullTracks{1, tmp1(q)}.yCoord(f+1);
                                end
                                dist2parent(q) = sqrt((px-x)^2 + (py-y)^2);
                            end
                            %                             dist2parent
                            if max(tmp1) < 0
                                %                                 disp('max(tmp1) < 0')
                                continue;
                            end
                            if lifespan <= 1 && f + 1 < endFrame
                                continue;
                            end
                            
                            if f + 2 > endFrame
                                continue;
                            end
                            
                            divCounter = divCounter + 1;
                            lineage_tb{divCounter, 1} = idRepetive;
                            lineage_tb{divCounter, 2} = f;
                            lineage_tb{divCounter, 3} = [pId, tmp]; % pId is the id of cell at frame f, tmp saves the ids of daughters at frame f+1
                            dp = find(dist2parent == min(dist2parent));
                            lineage_tb{divCounter, 4} = [tmp1(dp(1)), tmp1];
                            %                             disp('before: '); tmp1
                            %                             disp(['divCounter = ' num2str(divCounter)]);
                            %                             pos1
                            if divCounter > 1
                                for deltaDivCounter = 1:divCounter-1
                                    tmp2 = lineage_tb{divCounter - deltaDivCounter, 4};
                                    replacement = tmp2(1);
                                    tmp4 = tmp2(2:end);
                                    if ismember(tmp4, tmp1)
                                        tmp3 =  [setdiff(tmp1, tmp2) replacement];
                                        lineage_tb{divCounter, 4} = [tmp1(1), tmp3]; %
                                        %                                         disp(['@ the ' num2str(divCounter) '-th division, checking confirmed daughters of division ' ...
                                        %                                             num2str(divCounter - deltaDivCounter) ', with {1, 4} = ']);
                                        %                                         tmp2
                                        %                                         disp('Found intersection ... ');
                                        %                                         setdiff(tmp1, tmp2)
                                        %                                         disp('Going to replace daughters ');
                                    end
                                end
                            end
                            %                               disp('after: '); tmp1
                            
                        end
                        
                    end
                    end
                end
            
            if ~isempty(lineage_tb{1,1})
                lineage_tb = sortrows(lineage_tb, 2);
                if divCounter <= (endFrame - startFrame)*0.25
                    if size(lineage_tb, 1) == 1 && max(lineage_tb{1, 4}) < 0 % to filter out the merging case when tree height = 1
%                         disp('size(lineage_tb, 1) == 1 && max(lineage_tb{1, 4}) < 0')
                        continue;
                        else
                            lineageCounter = lineageCounter + 1;
                            totalLineages{1, lineageCounter} = lineage_tb;
                        end
                    end
                end
            end
            % totalLineages need "quality control" ...
            totalLineages2 = cell(1,3);
            counter = 0;
             for L = 1:size(totalLineages, 2)
                thisTree = totalLineages{1, L};
                if size(thisTree, 1) > (endFrame - startFrame)*1%4/60
%                     disp('size(thisTree, 1) > (endFrame - startFrame)*1')
                    continue;
                end
                % control, if the depth of one lineage is larger than the number of
                % hrs, this lineage must be tracked WRONG ....
                for depth = 1:size(thisTree, 1)
                    if endFrame - thisTree{depth, 2} < 1 % do not trust divisions at the last 2 frames ...
                        disp([' L = ' num2str(L) ' endFrame - thisTree{depth, 2} < 5'])
                        continue;
                    end
                    if thisTree{depth, 2} < 1 % do not trust divisions at the first 2 frames ...
                        continue;
                    end
                    counter = counter + 1;
                    totalLineages2(counter, 1) = thisTree(depth, 1);
                    totalLineages2(counter, 2) = thisTree(depth, 2); 
                    totalLineages2(counter, 3) = thisTree(depth, 4);
%                     totalLineages2(counter, 4) = thisTree(depth, 4);
                end
             end
            obj.lineages = totalLineages;
            obj.lineagesTB = totalLineages2;
        end
        
        function obj = identifyLineageOnTheGo(obj, startFrame, endFrame, includeDT) %colCurrentFrame
            obj.tracks = [obj.fullTracks];
            startingK = 1;
            if includeDT == 1
                obj.tracks = [obj.fullTracks obj.deadTracks];
                startingK = 2;
            end
            listSize = size( obj.tracks, 2);
            listSize2 = size( [obj.fullTracks ], 2);
          
            idxMatrix = zeros(listSize, endFrame-startFrame+1);
            
            for i = 1:listSize
%                 disp(['i=' num2str(i) '/' num2str(listSize) '  ' num2str(length(obj.tracks{1,i}.indices))]);
                idxMatrix(i,:) =  obj.tracks{1,i}.indices;
            end
            
            colCurrentFrame = idxMatrix(:, 1); % WAS 1, LOOKING AT THE FIRST COLUMN
            [a,b] = hist(colCurrentFrame,unique(colCurrentFrame));
%             a = a(2:end); b = b(2:end); % show a, b'
            order = find(a >= 2 );
            lineageCounter = 0;
            totalLineages = [];            

            for k = startingK:length(order)               
                idRepetive = b(order(k));
                pos = find(colCurrentFrame == idRepetive);
                idMatrixOfLineage = []; pos2tracks = [];
                for p = 1:length(pos)
                    idMatrixOfLineage = [idMatrixOfLineage; idxMatrix(pos(p), :)];
                    pos2tracks = [pos2tracks; pos(p)];
                end
                
%                 idMatrixOfLineage
                
%                 idMatrixOfLineage
                divCounter = 0;
                lineage_tb = cell(1, 4);
                disp(['Scanning ' num2str(k) '-th out of ' num2str(length(order)) ' lineages.']);
                
                for f = size(idMatrixOfLineage, 2)-1
                    if length(idMatrixOfLineage(:,f))-length(unique(idMatrixOfLineage(:,f))) == 0
                        disp('length(idMatrixOfLineage(:,f))-length(unique(idMatrixOfLineage(:,f))) == 0')
                        continue;
                    else
                        [aa, bb] = hist(idMatrixOfLineage(:,f), unique(idMatrixOfLineage(:,f)));
                        div = find(aa>=2);
                        disp(['f = ' num2str(f) ' numDivs = ' num2str(length(div))]);
                        for d = 1:length(div)
                            pId = bb(div(d));                            
                            pos1 = find(idMatrixOfLineage(:,f) == pId);
                            if length(pos1) < 2 % determining number of daughters, 1 means no div
%                                 disp('length(pos1)');
                                continue;
                            end
                            tmp = zeros(1, length(pos1));
                            tmp1 = tmp;
                            for q = 1:length(pos1)
                                tmp(q) = idMatrixOfLineage(pos1(q), f+1); % check next frame if id is consistant ...
                                tmp1(q) = idMatrixOfLineage(pos1(q), end); % the idx at the last frame will be used as final cell Ids ..
                            end
   
                            if length(unique(tmp)) == 1
%                                 disp('length(unique(tmp)) == 1')
                                continue;
                            end
                            dist2parent = zeros(size(tmp1)); 
                            for q = 1:length(pos1)
                                if tmp1(q) <= 0
                                    tmp1(q) = -(pos2tracks(pos1(q)) - listSize2);
                                    px = obj.deadTracks{1, -tmp1(q)}.xCoord(f); py = obj.deadTracks{1, -tmp1(q)}.yCoord(f);
                                    x = obj.deadTracks{1, -tmp1(q)}.xCoord(f + 1); y = obj.deadTracks{1, -tmp1(q)}.yCoord(f+1);
                                else
                                    px = obj.fullTracks{1, tmp1(q)}.xCoord(f); py = obj.fullTracks{1, tmp1(q)}.yCoord(f);
                                    x = obj.fullTracks{1, tmp1(q)}.xCoord(f + 1); y = obj.fullTracks{1, tmp1(q)}.yCoord(f+1);
                                end
                                dist2parent(q) = sqrt((px-x)^2 + (py-y)^2);
                            end 
                            dist2parent = round(dist2parent*10)/10;

                            divCounter = divCounter + 1;
                            lineage_tb{divCounter, 1} = idRepetive;
                            lineage_tb{divCounter, 2} = f;
                            lineage_tb{divCounter, 3} = [pId, tmp]; % pId is the id of cell at frame f, tmp saves the ids of daughters at frame f+1
                            
                            dp = find(dist2parent == min(dist2parent));

                            lineage_tb{divCounter, 4} = [tmp1(dp(1)), tmp1];                            
                        end
                    end
                end
                
                if divCounter > 0
                    lineageCounter = lineageCounter + 1;
                    totalLineages{1, lineageCounter} = lineage_tb;
                end
            end
            
              obj.lineages = totalLineages;
        end

        function mamut2 = GMMs_Res2Mamut(obj, res, mamut)
            mamut2 = res';
            if ~isempty(mamut)
                svIdx_res = [res.svIdx];
                svIdx_mamut = [mamut.svIdx];
                toAdd = setdiff(svIdx_mamut, svIdx_res);
                if ~isempty(toAdd)
                    counter = size(mamut2, 1);
                    for i = 1:size(toAdd, 2)
                        svId = toAdd(i);
                        pos = find(svId == svIdx_mamut);
                        for j = 1:length(pos)
                            counter = counter + 1;
                            mamut2(counter) = mamut(pos(j));
                        end
                    end
                end
            end
        end
        
        function [obj, case1, case2] = getTripolarDiv(obj, tripolar_input) 
            % case 1 saves lineage ids that are grouped as tripolar; 
            % case 2 saves lineage ids that cannot be grouped 
            case1 = []; case2 = [];
            obj.tripolar_lineage_id = [];
            for tri = 1:size(tripolar_input, 1)
                event_time = tripolar_input(tri, 1);
                event_coordinate = tripolar_input(tri, 2:3);
                thrDist = 20; % first detect the lineage of tripolar, using a large distance
                tripolar = 0;
                for L = 1:size(obj.lineages, 2)
                    lineage = obj.lineages{1, L};
                    candidate = lineage{end, 4}(1);
                    if candidate < 0 
                        candidate = max(lineage{end, 4}); % specify one cell which belongs to the lineage
                    end
                    if candidate < 0
%                         disp('candidate < 0');
                        continue;
                    end
                    candidate_xy = [obj.fullTracks{1, candidate}.xCoord(event_time) obj.fullTracks{1, candidate}.yCoord(event_time)];
                    if sqrt(sum((event_coordinate-candidate_xy)'.^2)) < thrDist
                        tripolar = L;
                        disp(tripolar);
                        break;
                    end
                end
                dt = obj.deadTracks;
                if tripolar == 0
                    case2 = [case2; tri L];
                    disp('This tripolar division cannot be wired ...');
                else
                    cells = [];
                    tripolarLineage = obj.lineages{1, tripolar};
                    for i = 1:size(tripolarLineage, 1)
                        if abs(tripolarLineage{i, 2} - event_time) > 5                           
                            continue;
                        end
                        cells = [cells tripolarLineage{i, 4}];
                    end
                    if isempty(cells)
                         case2 = [case2; tri tripolar];
                        continue;
                    end
                    pId = cells(1);
                    cells = unique(cells);
                    deadOnes = [];liveOnes = [];
                    for c = 1:length(cells)
                        if cells(c) > 0
                            liveOnes = [liveOnes; cells(c)];
                        else
                            Len = max(find(dt{1, -cells(c)}.indices > 0)) - min(find(dt{1, -cells(c)}.indices > 0));
                            deadOnes = [deadOnes; cells(c) Len+1];
                        end
                    end
                    deadOnes = sort(deadOnes, 2);
%                     deadOnes

                    cells = cells(min(find(cells>0)):max(find(cells>0)));
                    if ~isempty(deadOnes)
                        while size(cells, 2) < 3 && deadOnes(end, 2) > 100
                            cells = [cells deadOnes(end, 1)];
                            deadOnes(end, :) = [];
                        end
%                         if size(cells, 2) < 3 && deadOnes(end, 2) > 100 % we check if a long dTrack can compensate
%                             cells = [cells deadOnes(end, 1)];
%                         end
                    end
%                     cells
                    if cells(1) == 0
                        cells(1) = [];
                    end
%                     if ~ismember(pId, cells)
                        cells = [pId cells];
%                     end
                    triLineage{1,1} = pId;
                    triLineage{1,2} = event_time;
                    triLineage{1,4} = cells;
                    
                    disp(['sum(cells >0) = ' num2str(sum(cells >0))]);
%                     if sum(cells >0) == 4
                        obj.lineages{1, tripolar} = triLineage;
                        obj.tripolar_lineage_id = [obj.tripolar_lineage_id; tripolar];
%                     end
case1 = [case1; tri tripolar];
                end
            end
        end
        
        function obj = correctLineage(obj, startFrame, endFrame, slidingWindow, distance)
            listSize = size(obj.lineage, 2);
            for i = 1:listSize
                ddTracks = []; fTracks = [];
                currentLineage = obj.lineage{1, i};
                if isempty(currentLineage)
                    continue;
                end                
                m = [];
                for j = 1:size(currentLineage, 2)
                    if prod(currentLineage(j).indices) == 0
                        disp(['Merging case detected ... of lineage ' num2str(i) ', the ' num2str(j) '-th track.']);
                        ddTracks = [ddTracks currentLineage(j)];
                        obj.mapLineage2Tracks{1, i}(j) = 0;
                    else
                        m = [m; currentLineage(j).indices];
                        fTracks = [fTracks currentLineage(j)];
                    end
                end      
                % check if fTracks generate false positive ...
                if ~isempty(m) && size(m, 1)>1                                                                                                                 
                    div_t = max(find(m(1,:) == m(2, :)));           
                    dTracks = obj.deadTracks;                    
                    isFalsePositive = 'false';
                    for p = 1:size(dTracks, 2)
                        if sum(dTracks{1, p}.indices) == 0
                            continue;
                        end
%                         disp(['Processing dTrack ' num2str(p) ' ... with guessed false-positives ... ']);
                        dt = max(find(dTracks{1, p}.indices > 0)); % at time dt + 1, track is disappearing, meaning dt is the detected div time 
                        if dt > div_t || dt + slidingWindow < div_t
                            continue;
                        end
                        disp(['Detected div time is ' num2str(div_t) ', and the dTrack ends at ' num2str(dt) '.']);
                        deadX = dTracks{1, p}.xCoord(dt); deadY = dTracks{1, p}.yCoord(dt);
%                         disp(['deadX = ' num2str(deadX) ', thisY = ' num2str(deadY) ' at time ' num2str(dt)]);
%                         disp(['Processing lineage ' num2str(i) ' ... ']);

                        for sliding_t = dt + 1:div_t
                            thisX = fTracks(1).xCoord(sliding_t);
                            thisY = fTracks(1).yCoord(sliding_t);
%                             disp(['thisX = ' num2str(thisX) ', thisY = ' num2str(thisY) ' at time ' num2str(sliding_t)]);
                            dist = sqrt((deadX - thisX)^2 + (deadY - thisY)^2);
                            if dist < distance
                                disp('false-positive is detected!');
                                isFalsePositive = 'true';
                                break;
                            end
                        end  
                        
                        if isequal(isFalsePositive, 'true')
%                             disp('isFalsePositive');
                            break;
                        end
                        
                    end
                    %  ... DO SOMETHING ...to correct the false
                    %  positive
                    if isequal(isFalsePositive, 'true') && div_t < endFrame
                        d1XY = [fTracks(1).xCoord(div_t + 1) fTracks(1).yCoord(div_t+1)];
                        d2XY = [fTracks(2).xCoord(div_t + 1) fTracks(2).yCoord(div_t+1)];
                        dist1 = sqrt((deadX - d1XY(1))^2 + (deadY - d1XY(2))^2);
                        dist2 = sqrt((deadX - d2XY(1))^2 + (deadY - d2XY(2))^2);
                        toCorrect = 2;
                        if dist1 < dist2
                            toCorrect = 1;
                        end
                        fTracks(toCorrect) = fTracks(toCorrect).mergeDT2FT(dTracks{1, p});                       
                        dTracks{1,p} = dTracks{1,p}.deleteTrack();
                        dt2 = max(find(fTracks(1).indices(1:dt) == fTracks(2).indices(1:dt))); % double check if the corrected tracks generate false positive divisions
%                         disp(['dt2 = ' num2str(dt2) ', dt = ' num2str(dt)]);
                        if dt2 < dt 
                            for pp = 1:size(dTracks, 2)
                                if prod(dTracks{1, pp}.indices == 0)
                                    continue;
                                end                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
                                myX = fTracks(1).xCoord(dt2);
                                myY = fTracks(1).yCoord(dt2);
                                ddt2 = max(find(dTracks{1, pp}.indices > 0));
                                ddX = dTracks{1, pp}.xCoord(ddt2); ddY = dTracks{1, pp}.yCoord(ddt2);
                                if  sqrt((ddX - myX)^2 + (ddY - myY)^2)< distance
                                    disp(['Found with dTrack ' num2str(pp) ' with ' num2str(pp) '-th dTrack, ends at ' num2str(ddt2)'.']);
                                    d1XY = [fTracks(1).xCoord(dt2 + 1) fTracks(1).yCoord(dt2+1)];
                                    d2XY = [fTracks(2).xCoord(dt2 + 1) fTracks(2).yCoord(dt2+1)];
                                    dist1 = sqrt((ddX - d1XY(1))^2 + (ddY - d1XY(2))^2);
                                    dist2 = sqrt((ddX - d2XY(1))^2 + (ddY - d2XY(2))^2);
                                    toCorrect = 2;
                                    if dist1 < dist2
                                        toCorrect = 1;
                                    end                                    
                                    fTracks(toCorrect) = fTracks(toCorrect).mergeDT2FT(dTracks{1, pp});
                                    dTracks{1, pp} = dTracks{1, pp}.deleteTrack();                                    
                                end
                            end
                        end                        
                        mapping = obj.mapLineage2Tracks{1, i}; % these are the cell tracks (their row numbers) belonging to the i-th lineage
                        if ~isempty(mapping)
                            for mm = 1:min(2, size(mapping, 2))  % to fix in the future, so far we ONLY handle two tracks per lineage ...                               
                                rId = mapping(mm);
                                if rId > 0 && rId <= size(obj.fullTracks, 2)
                                    obj.fullTracks{1, rId} = fTracks(mm);
%                                     changed on March 9%                                     
%                                     obj.tracks{1, rId} = fTracks(mm);
                                end
                            end
                        end
                          obj.lineage{1, i} = [];
                    end
                end
            end
        end
        
        function [mergingCases, obj] = correctLineage2(obj, startFrame, endFrame, slidingWindow, thrDist)
            mergingCases = [];
            dTracks = obj.deadTracks;
            for i = 1:size(obj.lineages, 1)
                thisLineage = obj.lineages{1, i};
                for j = 1:size(thisLineage, 1)
                    div_t = thisLineage{j, 2};
                    div_x = obj.fullTracks{1, thisLineage{j, 3}(1)}.xCoord(div_t);
                    div_y = obj.fullTracks{1, thisLineage{j, 3}(1)}.yCoord(div_t);
                    %                 isFalsePositive = 'false';
                    for d = 1:size(dTracks, 2)
                        dead = dTracks{1,d};
                        dead_t = max(find(dead.indices > 0));
                        %                     disp(['dead_t = ' num2str(dead_t) ', div_t = ' num2str(div_t)]);
                        if isempty(dead_t)
                            continue;
                        end
                        if abs(dead_t - div_t) > slidingWindow % dead_t > div_t
                            continue;
                        end
                        if div_t - dead_t > slidingWindow
                            continue;
                        end
                        deadX = 0; deadY = 0; disappering_t = 0;
                        for sliding_t = max(1, div_t - slidingWindow):max(div_t, dead_t)
                            dead_x = dead.xCoord(sliding_t);
                            dead_y = dead.yCoord(sliding_t);
                            distance = sqrt((div_x-dead_x)^2+(div_y-dead_y)^2);
                            if distance < thrDist
                                disp(['Merging case detected ... of lineage ' num2str(obj.lineagesTB{i, 1}) ', the ' num2str(i) '-th row on the lineage table.']);
                                %                             isFalsePositive = 'true';
                                deadX = dead_x; deadY = dead_y;
                                disappering_t = sliding_t;
                                break;
                            end
                        end
                        if deadX*deadY~=0 %isequal(isFalsePositive,'true')
                            mergingCases = [mergingCases; i];
                            thisLineage{j, 1} = -1;
                            break;
                        end
                    end
                end
            end
        end
        
        function obj = correctLineage3(obj, startFrame, endFrame, slidingWindow, thrDist)
            dTracks = obj.deadTracks;           
            for i = 116%1:size(obj.lineages, 2)
                theLineage = obj.lineages{1, i};
                wrongDaughter = [];
                correctDaughter = cell(size(theLineage, 1), 2);
                for j = size(theLineage, 1):-1:1                 
                    div_t = theLineage{j, 2};
                    div_x = obj.fullTracks{1, theLineage{j, 4}(1)}.xCoord(div_t); 
                    div_y = obj.fullTracks{1, theLineage{j, 4}(1)}.yCoord(div_t);
                    correctDaughter{j, 2} = theLineage{j, 4}(1);
                    isFalsePositive = 'false';
                    daughters = [];
                    for dd = 2:length(theLineage{j, 4})
                        if ismember(theLineage{j, 4}(dd), wrongDaughter)
                            theLineage{j, 4}(dd) = -1;
                        end
                        daughters = [daughters theLineage{j, 4}(dd)];
                    end
%                     disp(['before: with j = ' num2str(j)]);
%                     daughters
                    for delta_j = 0:size(theLineage, 1)-j
                        confirmedDaughters = correctDaughter{j+delta_j, 1};
                        if ~isempty(confirmedDaughters)
                            inter = intersect(daughters, confirmedDaughters);
                            if ~isempty(inter)
                                for ii = 1:length(inter)
                                    pos = find(daughters == inter(ii));
                                    daughters(pos) = [];
                                end
                                daughters = [daughters correctDaughter{j+delta_j, 2}];
                            end
                        end
                    end
                    disp(['after: with j = ' num2str(j)]);
                    daughters
                    numDaughters = sum(daughters > 0);
                    if numDaughters <= 1
                        theLineage{j, 1} = -1;
                        continue;
                    end                    
                    for d = 1:size(dTracks, 2)    
                        dead = dTracks{1,d};
                        dead_t = max(find(dead.indices > 0));
                        dead_t
                        div_t
                        %                     disp(['dead_t = ' num2str(dead_t) ', div_t = ' num2str(div_t)]);
                        if isempty(dead_t)
                            disp('isempty(dead_t)')
                            continue;
                        end
                        %                         if abs(dead_t - div_t) > slidingWindow % dead_t > div_t
                        %                             continue;
                        %                         end
                        if dead_t > div_t
                            continue;
                        end
                        if div_t - dead_t > slidingWindow
%                             disp('div_t - dead_t > slidingWindow')
                            continue;
                        end
                        deadX = 0; deadY = 0;
                        for sliding_t = max(div_t, dead_t):-1:max(1, div_t - slidingWindow)
                            dead_x = dead.xCoord(sliding_t);
                            dead_y = dead.yCoord(sliding_t);
                            if dead_x*dead_y == 0
                                continue;
                            end
                            distance = sqrt((div_x-dead_x)^2+(div_y-dead_y)^2);
                            if distance < thrDist
%                                 disp(['Merging case detected ... of lineage ' num2str(theLineage{j, 1}) ', the ' num2str(j) '-th row on the lineage tree.']);
%                                 isFalsePositive = 'true';
                                deadX = dead_x; deadY = dead_y;
                                disp(['Found possible dTrack with dist = ' num2str(distance) ' at time ' num2str(sliding_t)]);  
                                break;
                            end
                        end  
                        if deadX*deadY~=0 %isequal(isFalsePositive,'true')
                            % if a mergine case is identified, we re-write the
                            % tracks below:
                           
                            deadX = dead.xCoord(max(sliding_t, dead_t)); deadY = dead.yCoord(max(sliding_t, dead_t));
%                             disp(['deadX = ' num2str(deadX) ' deadY = ' num2str(deadY) ' div_x = ' num2str(div_x) ' div_y = ' num2str(div_y) ...
%                                 ' at time = ' num2str(max(sliding_t, dead_t))]);
                            % now look at daughters
                            daughter_dist = zeros(numDaughters, 1); 
                            for num = 1:numDaughters 
                                if daughters(num) == -1
                                    daughter_dist(num, 1) = 100000;
                                else                                    
                                    xx = obj.fullTracks{1, theLineage{j, 4}(num + 1)}.xCoord(div_t + 1); % +1 because the first id is for parent
                                    yy = obj.fullTracks{1, theLineage{j, 4}(num + 1)}.yCoord(div_t + 1);
                                    dist = sqrt((xx-deadX)^2+(yy-deadY)^2);
                                    daughter_dist(num, 1) = dist;
                                end
                            end
                            minDaughter = find(daughter_dist == min(daughter_dist));
                            sameIds = intersect(dead.indices, obj.fullTracks{1, daughters(minDaughter)}.indices);
                           
                            if length(sameIds) < 0.2 * div_t
                                wrongDaughter = [wrongDaughter; daughters(minDaughter)];
%                                 obj.fullTracks{1, daughters(minDaughter)} = mergeDT2FT(obj.fullTracks{1, daughters(minDaughter)}, dead);
                                dTracks{1,d}.indices = zeros(1, endFrame-startFrame+1); % erase the dTrack that has been used ...
                                theLineage{j, 4}(minDaughter + 1) = -1;
                                theLineage{j, 1} = -1; % set lineage Id = 0, as merging case ..
                                isFalsePositive ='true';
                                disp(['A merging case detected, we are going to work on the ' num2str(d) '-th dTrack ...']);
                                disp(['Inserting the ' num2str(d) '-th dTrack ... to the ' num2str(daughters(minDaughter)) '-th fTrack.']);
                                break; % the loop for dTrack
                            end
                                                     
                        end    
                    end
                    % no merging cases ...
                    if isequal(isFalsePositive, 'false')
                        disp(['No merging case for the ' num2str(j) '-th division']);
                        correctDaughter{j, 1} = setdiff(daughters, -1);
                    end
                end
%                 obj.lineages{1, i} = theLineage; % updating the lineage
            end            
        end
        
        function obj = daughter2daughter_dist(obj, startFrame, endFrame, thrDist)
             for i = 29%1:size(obj.lineages, 2)
                 theLineage = obj.lineages{1, i};
                 for j = size(theLineage, 1):-1:1
                     daughter2parent{j, 2} = theLineage{j, 4}(1);
                     daughter2parent{j, 1} = theLineage{j, 4}(2:end);
                 end
                for j = size(theLineage, 1):-1:1                 
                    div_t = theLineage{j, 2};
                    isFalsePositive = 'false';
                    daughters = theLineage{j, 4}(2:end);
                     for delta_j = 0:size(theLineage, 1)-j
                        confirmedDaughters = daughter2parent{j+delta_j, 1};
                        if ~isempty(confirmedDaughters)
                            inter = intersect(daughters, confirmedDaughters);
                            if ~isempty(inter)
                                for ii = 1:length(inter)
                                    pos = find(daughters == inter(ii));
                                    daughters(pos) = [];
                                end
                                daughters = [daughters daughter2parent{j+delta_j, 2}];
                            end
                        end
                     end
                     numDaughters = length(daughters);
                     daughter_xy = zeros(numDaughters, 2);
                     for num = 1:numDaughters
                         xx = obj.fullTracks{1, theLineage{j, 4}(daughters(num) + 1)}.xCoord(div_t + 1); % +1 because the first id is for parent
                         yy = obj.fullTracks{1, theLineage{j, 4}(daughters(num) + 1)}.yCoord(div_t + 1);
                         dist = sqrt((xx-deadX)^2+(yy-deadY)^2);
                         daughter_xy(num, :) = [xx yy];
                     end
                     d2d_dist = sqrt(sum((daughter_xy(:, 1) - daughter_xy(:, 2)).^2)); % now only taking into account 2 daughtrs
                     if d2d_dist < thrDist
                     end
                end
             end
        end
        
        function obj = breakLineageOnTheGo(obj)
            dTracks = obj.deadTracks;
            for i =1:size(obj.lineages, 2)
                theLineage = obj.lineages{1, i};
                div_t = theLineage{1, 2} + 1;
                pId = theLineage{1, 4}(1);
                dId = setdiff(theLineage{1, 4}, pId);    
                dId = dId(1);
                if pId == 0
                    theLineage{1, 1} = -1;
                    continue;
                end
                if pId < 0
                    div_x = obj.deadTracks{1, -pId}.xCoord(div_t);
                    div_y = obj.deadTracks{1, -pId}.yCoord(div_t);
                else
                    div_x = obj.fullTracks{1, pId}.xCoord(div_t);
                    div_y = obj.fullTracks{1, pId}.yCoord(div_t);
                end
                dTracks2choose = [];
                for d = 1:size(dTracks, 2)
                    dead = dTracks{1,d};
                    dead_t = max(find(dead.indices > 0)); 
%                          disp(['dead_t = ' num2str(dead_t) ', div_t = ' num2str(div_t)]);
                    if isempty(dead_t)
                        continue;
                    end
                    if dead_t >= div_t
                        continue;
                    end

                    dead_x = dead.xCoord(dead_t); dead_y = dead.yCoord(dead_t);
                    deadDist = sqrt( (div_x-dead_x)^2 + (div_y-dead_y)^2 );           

                    if deadDist < 10 && dId > 0
                        dTracks2choose = [dTracks2choose; d dead_t];
                    end
                end
                candidate = zeros(size(dTracks2choose, 1), 3);
                for c = 1:size(dTracks2choose, 1)
                    dead_x = dTracks{1, dTracks2choose(c, 1)}.xCoord(max(1, dTracks2choose(c, 2)));
                    dead_y = dTracks{1, dTracks2choose(c, 1)}.yCoord(max(1, dTracks2choose(c, 2))); % middle
                    x0 = dTracks{1, dTracks2choose(c, 1)}.xCoord(max(1, dTracks2choose(c, 2)-1));
                    y0 = dTracks{1, dTracks2choose(c, 1)}.yCoord(max(1, dTracks2choose(c, 2)-1)); % early
                    P1 = [div_x div_y]; P0 = [dead_x dead_y]; P2 = [x0 y0];
                    ang = atan2(abs(det([P2-P0;P1-P0])),dot(P2-P0,P1-P0))*180/pi;
                    candidate(c, :) = [dTracks2choose(c, :) ang];
                end
                candidate = sortrows(candidate, 3, 'descend'); 
                if ~isempty(candidate) 
%                     dd = candidate(1, 1); dead_t = candidate(1, 2);
                    disp(['Breaking the ' num2str(dId) '-th fTrack ...']);
                    obj.fullTracks{1, dId}.indices(1:end-1) = 0;
                    obj.fullTracks{1, dId}.xCoord(1:end-1) = 0;
                    obj.fullTracks{1, dId}.yCoord(1:end-1) = 0;
%                     obj.fullTracks{1, dId}.svIdx(1:end-1) = 0;
                end             

            end
        end
        
        function [obj, numCorrected, numDiv, info] = correctLineageOnTheGo(obj, death2birth_dist)
            numCorrected = 0;
            numDiv = 0;
            info = [];
            dTracks = obj.deadTracks;
            for i = 1:size(obj.lineages, 2)
                              theLineage = obj.lineages{1, i};
                div_t = theLineage{1, 2} + 1;
                pId = theLineage{1, 4}(1);
                
                % commented
                dId = setdiff(theLineage{1, 4}, pId);    
%                 dId = dId(1);
%                 dId = unique(theLineage{1, 4});
                
                if pId == 0 || length(dId) > 2
                    theLineage{1, 1} = -1;
                    continue;
                end
                numDiv  = numDiv + 1;
                if pId < 0
                    div_x = obj.deadTracks{1, -pId}.xCoord(div_t);
                    div_y = obj.deadTracks{1, -pId}.yCoord(div_t);
                else
                    div_x = obj.fullTracks{1, pId}.xCoord(div_t);
                    div_y = obj.fullTracks{1, pId}.yCoord(div_t);
                end

                dTracks2choose = [];
                for d = 1:size(dTracks, 2)
                    dead = dTracks{1,d};
                    dead_t = max(find(dead.indices > 0)); 
%                          disp(['dead_t = ' num2str(dead_t) ', div_t = ' num2str(div_t)]);
                    if isempty(dead_t)
%                         disp('dead_t is empty')
                        continue;
                    end
                    if dead_t >= div_t
                        continue;
                    end

                    dead_x = dead.xCoord(dead_t); 
                    dead_y = dead.yCoord(dead_t);
                    deadDist = sqrt( (div_x-dead_x)^2 + (div_y-dead_y)^2 );       

                    if round(deadDist) <= death2birth_dist && min(dId) > 0
%                         div_x
%                         div_y
%                         dead_x
%                         dead_y
%                         deadDist
                        dTracks2choose = [dTracks2choose; d dead_t];
                    end
                end
   
                
                candidate = zeros(size(dTracks2choose, 1), 3);
                for c = 1:size(dTracks2choose, 1)
                    dead_x = dTracks{1, dTracks2choose(c, 1)}.xCoord(max(1, dTracks2choose(c, 2)));
                    dead_y = dTracks{1, dTracks2choose(c, 1)}.yCoord(max(1, dTracks2choose(c, 2))); % middle
                    x0 = dTracks{1, dTracks2choose(c, 1)}.xCoord(max(1, dTracks2choose(c, 2)-1));
                    y0 = dTracks{1, dTracks2choose(c, 1)}.yCoord(max(1, dTracks2choose(c, 2)-1)); % early
                    P1 = [div_x div_y]; P0 = [dead_x dead_y]; P2 = [x0 y0];
                    ang = atan2(abs(det([P2-P0;P1-P0])),dot(P2-P0,P1-P0))*180/pi;
                    candidate(c, :) = [dTracks2choose(c, :) ang];
                end
                candidate = sortrows(candidate, 3, 'descend'); 

                if ~isempty(candidate)
                    
                    dd = candidate(1,1);  dead_t = candidate(1, 2);
                    closeId = 1;
%                     dId
                    if length(dId) > 1
%                         added

                        dead_x = obj.deadTracks{1, dd}.xCoord(dead_t);
                        dead_y = obj.deadTracks{1, dd}.yCoord(dead_t);
                        x1 = obj.fullTracks{1, dId(1)}.xCoord(div_t);
                        y1 = obj.fullTracks{1, dId(1)}.yCoord(div_t);
                        x2 = obj.fullTracks{1, dId(2)}.xCoord(div_t);
                        y2 = obj.fullTracks{1, dId(2)}.yCoord(div_t);
                        myDist = pdist2([dead_x dead_y], [x1 y1; x2 y2]);
%                                                                     [dead_x dead_y]
%                                                                     [x1 y1; x2 y2]
                        closeId = find(myDist == min(myDist));
                    end
                
                    disp(['Updating the ' num2str(dId(closeId)) '-th fTrack with ' num2str(dd) '-th dTrack...']);

                    obj.fullTracks{1, dId(closeId)}.indices(1:dead_t) = obj.deadTracks{1, dd}.indices(1:dead_t);
                    obj.fullTracks{1, dId(closeId)}.xCoord(1:dead_t) = obj.deadTracks{1, dd}.xCoord(1:dead_t);
                    obj.fullTracks{1, dId(closeId)}.yCoord(1:dead_t) = obj.deadTracks{1, dd}.yCoord(1:dead_t);
                    obj.fullTracks{1, dId(closeId)}.svIdx(1:dead_t) = obj.deadTracks{1, dd}.svIdx(1:dead_t);
                    info = [info; obj.deadTracks{1, dd}.xCoord(dead_t) obj.deadTracks{1, dd}.yCoord(dead_t)];
                    % added 
                    if dead_t < length(obj.fullTracks{1, dId(closeId)}.indices) - 1
                        obj.fullTracks{1, dId(closeId)}.indices(end - 1) = 0;
                    end
                    numCorrected = numCorrected + 1;
                    % end
                end             

            end
        end
        

        function [obj, corrected] = correctLineage4(obj, startFrame, endFrame, slidingWindow, thrDist, update)
            dTracks = obj.deadTracks;
            obj.dT2fT = zeros(size(dTracks, 2), 1);
            corrected = [];
            for i = 1:size(obj.lineages, 2)
                disp(['Scanning the ' num2str(i) '-th lineage out of ' num2str(size(obj.lineages, 2)) ' lineages ...'])
                theLineage = obj.lineages{1, i};
                wrongDaughter = [];
                daughter2parent = cell(size(theLineage, 1), 2);
                for j = size(theLineage, 1):-1:1       
                    if theLineage{j, 1} < 0
                        continue;
                    end
                
%                   disp(['----------- examining division j = ' num2str(j) ' ----------------']);
                    div_t = theLineage{j, 2}-1;
                    pId = theLineage{j, 4}(1);
                    if pId == 0
                        theLineage{j, 1} = -1;
                        continue;
                    end
                    if pId < 0
                        div_x = obj.deadTracks{1, -pId}.xCoord;
                        div_y = obj.deadTracks{1, -pId}.yCoord;
                    else
                        div_x = obj.fullTracks{1, pId}.xCoord;
                        div_y = obj.fullTracks{1, pId}.yCoord;
                    end

                    daughter2parent{j, 2} = theLineage{j, 4}(1);
                    if daughter2parent{j, 2} < 0
                        daughter2parent{j, 2} =  obj.dT2fT(-daughter2parent{j, 2});
                    end
                    isFalsePositive = false;
                    daughters = [];
                    for dd = 2:length(theLineage{j, 4})
                        if theLineage{j, 4}(dd) < 0 && obj.dT2fT(-theLineage{j, 4}(dd)) > 0
                            theLineage{j, 4}(dd) =  obj.dT2fT(-theLineage{j, 4}(dd));
                        end
                        daughters = [daughters theLineage{j, 4}(dd)];
                    end
%                     disp(['before: with j = ' num2str(j)]);  daughters
                    for delta_j = size(theLineage, 1)-j:-1:0
                        confirmedDaughters = daughter2parent{j+delta_j, 1};
                        if ~isempty(confirmedDaughters)
                            inter = intersect(daughters, confirmedDaughters);
                            if ~isempty(inter) && length(inter) > 1
                                for ii = 1:length(inter)
                                    pos = find(daughters == inter(ii));
                                    pos = pos(1);
                                    daughters(pos) = [];
                                end
                                daughters = [daughters daughter2parent{j+delta_j, 2}];
                            end
                        end
                    end
                    daughters = unique(daughters);
                    
%                     disp(['after: with j = ' num2str(j)]);  daughters
%                     numDaughters = sum(daughters > 0);
                    numDaughters = length(daughters);
                    
                    % Adding one more condition here: if two daughters are
                    % always next to each other for consistantly n (=5)
                    % frames, highly likely it is caused by
                    % over-segmentation                    
                    oversegmentation = 0;
                    daughters_xy = cell(numDaughters, 1);
                    for num = 1:numDaughters
                        if daughters(num) < 0
                            xx = obj.deadTracks{1, -daughters(num)}.xCoord(div_t + 1: min(div_t + 10, endFrame-startFrame+1)); % +1 because the first id is for parent
                            yy = obj.deadTracks{1, -daughters(num)}.yCoord(div_t + 1: min(div_t + 10, endFrame-startFrame+1));
                        else     
                            xx = obj.fullTracks{1, daughters(num)}.xCoord(div_t + 1: min(div_t + 10, endFrame-startFrame+1)); % +1 because the first id is for parent
                            yy = obj.fullTracks{1, daughters(num)}.yCoord(div_t + 1: min(div_t + 10, endFrame-startFrame+1));
                        end
                        daughters_xy{num, 1} = [xx; yy]';
                    end

                    if numDaughters == 2
                        data1 = daughters_xy{1, 1};
                        data2 = daughters_xy{2, 1};
                        dist2each = sqrt(sum((data1-data2).^2, 2));
                        if dist2each(end) < 5
                            oversegmentation = 1;
                        end
                        disp(['daughters: ' num2str(daughters(1)) ' ' num2str(daughters(2))]);
                    end

                    disp(['over-segmentation = ' num2str(oversegmentation)]);             
                    
                    if  prod(daughters) == 0 || numDaughters <= 1 || oversegmentation == 1 || numDaughters > 50 % was <= 1
                        theLineage{j, 1} = -1;
                        continue;
                    end  
    
                    for d = 1:size(dTracks, 2)    
                        dead = dTracks{1,d};
                        dead_t = max(find(dead.xCoord > 0));
%                         disp(['dead_t = ' num2str(dead_t) ', div_t = ' num2str(div_t)]);
                        if isempty(dead_t)
%                             disp('isempty(dead_t)')
                            continue;
                        end
                        if length(find(dead.xCoord(1:dead_t) > 0)) <= 2 % dead track also needs a certain length
%                             disp('length(find(dead.xCoord(1:dead_t) > 0)) <= 2')
                            continue;
                        end
                        if dead_t <= 1
                            continue;
                        end
                        
                        if dead_t > div_t + 2 % here 2 frames considered after fake cells generated, might have 2 frames extra (for dead tracks)
%                             disp('dead_t > div_t + 2');
                            continue;
                        end
                        if div_t - dead_t > slidingWindow
%                             disp('div_t - dead_t > slidingWindow')
                            continue;
                        end
                        deadX = 0; deadY = 0;
%                         disp(['dead_t = ' num2str(dead_t) ', div_t = ' num2str(div_t)]);
                        
                        for sliding_t = max(div_t, dead_t):-1:max(1, div_t - slidingWindow)
                            dead_x = dead.xCoord(sliding_t);
                            dead_y = dead.yCoord(sliding_t);
                            if dead_x*dead_y == 0
%                                 disp('dead_x*dead_y == 0');
                                continue;
                            end
                            distance = floor(sqrt((div_x(sliding_t)-dead_x)^2+(div_y(sliding_t)-dead_y)^2));
%                             disp(['distance = ' num2str(distance) ' with deadx = ' num2str(dead_x) ', deady = ' num2str(dead_y) ', div_x = ' num2str(div_x(div_t)) ', div_y = ' num2str(div_y(div_t))]);
                            if distance <= thrDist
%                                 disp(['Merging case detected ... of lineage ' num2str(theLineage{j, 1}) ', the ' num2str(j) '-th row on the lineage tree.']);
%                                 isFalsePositive = 'true';
                                deadX = dead_x; deadY = dead_y;
%                                 disp(['Found possible dTrack with dist = ' num2str(distance) ' at time ' num2str(sliding_t)]);  
                                break;
                            end
                        end  
                        if deadX*deadY~=0 %isequal(isFalsePositive,'true')
                            % if a mergine case is identified, we re-write the
                            % tracks below:
                           
                            deadX = dead.xCoord(max(sliding_t, dead_t)); deadY = dead.yCoord(max(sliding_t, dead_t));
%                             disp(['deadX = ' num2str(deadX) ' deadY = ' num2str(deadY) ' div_x = ' num2str(div_x) ' div_y = ' num2str(div_y) ...
%                                 ' at time = ' num2str(max(sliding_t, dead_t))]);
%                             now look at daughters
                            daughter_dist = zeros(numDaughters, 1); 
                            daughter_xy = zeros(numDaughters, 2);
                            daughter_angle = zeros(numDaughters, 1);
                            disp(['j = ' num2str(j) ' with numDaughters = ' num2str(numDaughters)]);
                            for num = 1:numDaughters  
%                                 if daughters(num) < 0
%                                     xx = obj.deadTracks{1, -daughters(num)}.xCoord(div_t + 1); % +1 because the first id is for parent
%                                     yy = obj.deadTracks{1, -daughters(num)}.yCoord(div_t + 1);   
%                                 else
%                                     xx = obj.fullTracks{1, daughters(num)}.xCoord(div_t + 1); % +1 because the first id is for parent
%                                     yy = obj.fullTracks{1, daughters(num)}.yCoord(div_t + 1);
%                                 end
                                data = daughters_xy{num, 1};
                                xx = data(1,1); yy = data(1, 2);
                                dist = sqrt((xx-deadX)^2+(yy-deadY)^2);
                                daughter_dist(num, 1) = dist;
                                daughter_xy(num, :) = [xx yy];    
%                                 max(sliding_t, dead_t)
                                P0 = [deadX deadY]; P2 = [xx yy]; P1 = [dead.xCoord(max(sliding_t, dead_t)-1) dead.yCoord(max(sliding_t, dead_t)-1)];
                                daughter_angle(num) = atan2(abs(det([P2-P0;P1-P0])),dot(P2-P0,P1-P0))*180/pi;
                            end
                            disp(['daughters: ' num2str(daughters(1)) ' ' num2str(daughters(2))]);
%                             disp(['theLineage{j, 4}(num + 1): ' num2str(theLineage{j, 4}(2)) '  ' num2str(theLineage{j, 4}(3))]);
%                             daughter_dist
                            
                            minDaughter = find(daughter_dist == min(daughter_dist));
                            minDaughter = minDaughter(1);
                            if max(daughter_dist) - min(daughter_dist) < 5 % if both daughters are too close to the dTrack, we use angle
%                                 disp('both daughters are too close to the dTrack... using angle to deternine');
                                 minDaughter = find(daughter_angle == max(daughter_angle));
                                 minDaughter = minDaughter(1);
                            end
%                             disp(['length(minDaughtre) = ' num2str(minDaughter)]);
%                             disp(['daughter_xy(1, :) = ' num2str(daughter_xy(1, 1)) ' ' num2str(daughter_xy(1, 2))]);
%                             disp(['daughter_xy(2, :) = ' num2str(daughter_xy(2, 1)) ' ' num2str(daughter_xy(2, 2))]);
%                             disp(['daughter_dist = ' num2str(daughter_dist(1)) ' ' num2str(daughter_dist(2))]);
%                             disp(['daughter_angle(1, :) = ' num2str(daughter_angle(1)) ]);
%                             disp(['daughter_angle(2, :) = ' num2str(daughter_angle(2)) ]);
%                             disp(['First found minDaughter = ' num2str(minDaughter) '  and ' num2str(daughters(minDaughter)) ]);
                              
                             if daughters(minDaughter) > 0
                                 sameIds = intersect(dead.indices, obj.fullTracks{1, daughters(minDaughter)}.indices);
                             else
                                 sameIds = intersect(dead.indices, obj.deadTracks{1, -daughters(minDaughter)}.indices);
                             end
%                              sameIds
                            disp(['Found dTrack with id = ' num2str(d) ' at j = ' num2str(j)]);
                            
%                             daughter_xy = daughter_xy(1:2,:)';
%                             if sqrt(sum((daughter_xy(:, 1) - daughter_xy(:, 2)).^2)) > 20 % just after division, ... daughter-to-daughter distance... 
% %                                 disp(['The ' num2str(j) '-th row, two daughters are too far .. with dist = ' num2str(sqrt(sum((daughter_xy(:, 1) - daughter_xy(:, 2)).^2)))]);                                
%                                 %if the two daughter cells are far away
%                                 %from each other (dist > 10 pixels)
%                                 sameIds = 0;
%                             end                            
                            
                            if size(sameIds, 2) < 0.95 * div_t && ~ismember(-d, daughters)
                                disp(['length(sameIds) = ' num2str(size(sameIds, 2)), ' div_t = ' num2str(div_t)]);
%                                 disp(['daughter_dist = ' num2str(daughter_dist(1)) ' ' num2str(daughter_dist(2))]);
                                wrongDaughter = [wrongDaughter; daughters(minDaughter)];
                                if update == 1
                                    if daughters(minDaughter) > 0
                                        obj.fullTracks{1, daughters(minDaughter)} = mergeDT2FT(obj.fullTracks{1, daughters(minDaughter)}, dead);
                                    else
                                        obj.deadTracks{1, -daughters(minDaughter)} = mergeDT2DT(obj.deadTracks{1, -daughters(minDaughter)}, dead);
                                    end
                                end
                                dTracks{1,d}.indices = zeros(1, endFrame-startFrame+1); % erase the dTrack that has been used ...
%                                 theLineage{j, 4}(minDaughter + 1) = -1;
                                theLineage{j, 1} = -1; % set lineage Id = 0, as merging case ..
                                isFalsePositive = true;
%                                 disp(['A merging case detected, we are going to work on the ' num2str(d) '-th dTrack ...']);
                                 obj.dT2fT(d) = daughters(minDaughter);
%                                 dT2fT(d)
                                disp(['Inserting the ' num2str(d) '-th dTrack ... to the ' num2str(daughters(minDaughter)) '-th fTrack.']);
                                break; % the loop for dTrack
                            end
                                                     
                        end    
                    end
                    % no merging cases ...
                    if isFalsePositive
                        corrected = [corrected; i];
                    else
                        disp(['No merging case for the ' num2str(j) '-th division']);
                        if theLineage{j, 4}(1) < 0
                            theLineage{j, 4}(1) = obj.dT2fT(-theLineage{j, 4}(1));
                        end
%                         daughter2parent{j, 1} = setdiff(daughters, -1);
                    end
                     daughter2parent{j, 1} = setdiff(daughters, -1);
                end
                if update == 1
                    obj.lineages{1, i} = theLineage; % updating the lineage
                end
                  
            end     
            obj.deadTracks = dTracks;
        end
        
        function obj = rearrangeLineage(obj)
            for i = 1:size(obj.lineages, 2)                      
                disp(['Re-arranging lineage ' num2str(i) '/' num2str(size(obj.lineages, 2) )]);
                target = obj.lineages{1, i};
                parent2daughters = cell(size(target, 1), 2);                
                for j = size(target, 1):-1:2
                    parent2daughters{j,1} = target{j,4}(1);
                    parent2daughters{j,2} = target{j,4}(2:end);                    
                    for delta_j = j-1:-1:1
                        inds = target{delta_j,4}(2:end);
                        inter = intersect(inds, parent2daughters{j,2});
                        if ~isempty(inter)&& size(inter, 2) > 1
                            s = [parent2daughters{j,1} setdiff(inds, parent2daughters{j,2})];
                            target{delta_j,4} = [target{delta_j,4}(1) s];
                        end
                    end
                end
                obj.lineages{1, i} = target;
            end
            
        end
        
        function lineages = getLineageFromLR(obj, lR, lineageCorrWindowSize, thr)
            lineages = cell(1, 1);%size(lR{1,1}, 1)
            size(lineages, 2)               
            
            for f = 1:size(lR, 2)
                firstLineages = lR{1,f};
                if ~isempty(firstLineages{1, 1})
                    disp('not empty')
                    break;
                end                   
            end
            firstLineageTimer = f;
            
            counter = 0;
            for L = 1:size(firstLineages, 1)
                ids = firstLineages{L, 3};
                if  prod(ids <= thr) > 0
                    counter = counter + 1;
                    lineages{1, counter}(1) = firstLineages(L, 1);
                    lineages{1, counter}(2) = firstLineages(L, 2);
                    lineages{1, counter}(3) = firstLineages(L, 3);
                end
            end
            
            size(lineages, 2)
            if size(lR, 2) > 1
                for i = firstLineageTimer+1:size(lR, 2)
                    disp(['lR i = ' num2str(i)]);
                    curLineages = lR{1,i};
                    if isempty(curLineages{1,1})
                        continue;
                    end
                    for j = 1:size(curLineages, 1)
                        curLineages{j, 2} = curLineages{j, 2} + (i-1)*lineageCorrWindowSize - 1;
                        parent = curLineages{j, 3}(1);
                        flag = 0;
                        
                        for k = 1:size(lineages, 2)
                            ll = lineages{1, k};
                            if ismember(parent, ll{1, 3}) && parent <= thr
                                s = size(ll, 1);
                                ll(s + 1, 1) = curLineages(j, 1);
                                ll(s + 1, 2) = curLineages(j, 2);
                                ll(s + 1, 3) = curLineages(j, 3);
                                flag = 1;
                                %                                 disp(['adding k = ' num2str(k) ', adding j = ' num2str(j)]);
                                %                                 size(ll, 1)
                                lineages{1, k} = ll;
                                break;
                            end
                        end
                        
                        if flag == 0
                            adding = cell(1,3);
                            adding(1) = curLineages(j, 1);
                            adding(2) = curLineages(j, 2);
                            adding(3) = curLineages(j, 3);
                            lineages{1, size(lineages, 2) + 1} = adding;
                        end
                        
                        
                    end
                end
            end
            % to generate lineage tree structure
            for i = 1:size(lineages, 2)
                curLineage = lineages{1, i};
                if size(curLineage, 1) < 2
                    continue;
                end
                parent2daughter = cell(size(curLineage, 1), 2);
                
                for depth = size(curLineage, 1):-1:1
                    parent2daughter{depth, 1} = curLineage{depth, 3}(1);
                    parent2daughter{depth, 2} = curLineage{depth, 3}(2:end);
                    for delta_depth = depth - 1:-1:1
                        inter = intersect(curLineage{delta_depth, 3}(2:end),  parent2daughter{depth, 2});
                        if isempty(inter)
                            continue;
                        end
                        %                         setdiff(curLineage{depth, 3}(2:end), inter);
                        curLineage{delta_depth, 3} = [curLineage{delta_depth, 3}(1) parent2daughter{depth, 1} setdiff(curLineage{delta_depth, 3}(2:end), inter)];
                    end
                end
                lineages{1, i} = curLineage;
            end
            
        end
        
        function obj = getLineageTB(obj, startFrame, endFrame, initial, final)
             counter = 0; lineagesTB = cell(1, 3);
             for L = 1:size(obj.lineages, 2)
                thisTree = obj.lineages{1, L};
                for depth = 1:size(thisTree, 1)
                    if final == 1 && endFrame - thisTree{depth, 2} < 3 % do not trust divisions at the last 2 frames ...
                        continue;
                    end
                    if initial == 1 && thisTree{depth, 2} < 6 % do not trust divisions at the first 2 frames ...
                        continue;
                    end
                    if thisTree{depth, 1} < 0 
                        continue;
                    end

                    if length(thisTree{depth, 4}) > 3 % was 4
                        continue;
                    end 
                    
                    counter = counter + 1; 
                    lineagesTB{counter, 1} = thisTree{depth, 4}(1);
                    lineagesTB(counter, 2) = thisTree(depth, 2); 
                    tmp = thisTree{depth, 4};                   
                    tmp = [tmp(1) unique(tmp(2:end))];
                    % */ Li added this on June 24
%                     pos = find(tmp > 0);
%                     if length(pos) < 4
%                         continue;
%                     end
                    %*/ till here
                    
                    lineagesTB(counter, 3) = {tmp};

                end
             end
            obj.lineagesTB = lineagesTB;
        end
        
        function obj = identifyDiv(obj, startFrame, endFrame)
            obj.tb_div = [];
            listSize = size(obj.lineage, 2);
            for i = 1:listSize
                currentLineage = obj.lineage{1, i};
                if isempty(currentLineage)
                    continue;
                end
                m = []; 
                for j = 1:size(currentLineage, 2)
                    m = [m; currentLineage(j).indices];
                end               
                div_t = max(find(m(1,:) == m(2, :))) + 1;
                obj.tb_div
                
                if div_t <= endFrame
                    m(:, div_t)'
                    obj.tb_div = [obj.tb_div; div_t m(:, div_t)' i]; % last col saves which lineage the div belongs to
                end
            end
            obj.tb_div = sortrows(obj.tb_div);
        end 
        
        function lR = getPositiveLineage2(obj, endFrame, initial, final)
            counter = 0; 
            lR = cell(1,1);
             for L = 1:size(obj.lineages, 2)
                thisTree = obj.lineages{1, L};
                for depth = 1:size(thisTree, 1)
                    if final == 1 && endFrame - thisTree{depth, 2} <= 2 % do not trust divisions at the last 2 frames ...
                        continue;
                    end
                    if initial == 1 && thisTree{depth, 2} <= 2 % do not trust divisions at the first 2 frames ...
                        continue;
                    end
                    if thisTree{depth, 1} < 0 
                        continue;
                    end
                    
                    if thisTree{depth, 4}(1) <=0 
                        continue;
                    end

                    if length(thisTree{depth, 4}) > 3 % was 4
                        continue;
                    end                     
                              
                    tmp = thisTree{depth, 4};     
                    
                    if length(unique(tmp)) == length(tmp)
                        continue;
                    end
                    counter = counter + 1;
                    disp(['counter = ' num2str(counter)]);
                    tmp = [tmp(1) unique(tmp(2:end))];
                    flag = 1;
                    lR{counter, 1} = flag;
                    lR{counter, 2} = thisTree{depth, 2};
%                     lR{counter, 3} = useful;
                    lR{counter, 3} = tmp;                    
                    
                end
             end
                
        end
        
        function lR = getPositiveLineage(obj, endFrame, initial, final, validIdx, validDTidxMatrix, ctr)
            counter = 0; 
            lR = cell(1,1);
             for L = 1:size(obj.lineages, 2)
                thisTree = obj.lineages{1, L};
                for depth = 1:size(thisTree, 1)
                    if final == 1 && endFrame - thisTree{depth, 2} <= 10 % do not trust divisions at the last 2 frames ...
                        continue;
                    end
                    if initial == 1 && thisTree{depth, 2} <= 6 % do not trust divisions at the first 2 frames ...
                        continue;
                    end
                    if thisTree{depth, 1} < 0 
                        continue;
                    end
                    
                    if thisTree{depth, 4}(1) <=0 
                        continue;
                    end

                    if length(thisTree{depth, 4}) > 4 % was 4
                        continue;
                    end 
                    
                              
                    tmp = thisTree{depth, 4};     
           
                    useful = [];
                    if prod(ismember(unique(tmp(2:end)), validIdx)) == 0
                        %                         continue;
                       useful = [];
                       disp(['Depth = ' num2str(depth) ', prod(ismember(unique(tmp(2:end)), validIdx)) == 0']);
                       tmp2 = [tmp(1) unique(tmp(2:end))];
                       for u = 1:size(tmp2, 2)
                           if ismember(tmp2(u), validIdx)
                               useful = [useful tmp2(u)];
                           else
                               if ismember(tmp2(u), validDTidxMatrix)
                                   disp(['u = ' num2str(u) ' tmp2(u) = ' num2str(tmp2(u))]);
                                   dRows = find(tmp2(u) == validDTidxMatrix);
                                   if ~isempty(dRows)   
                                       m = zeros(1, length(dRows));
                                        for dr = 1:length(dRows)
                                            m(dr) = ctr(dRows(dr));
                                        end
                                        dd = find(m == max(m));
                                        useful = [useful -dRows(dd(1))];
                                   end
                               end
                           end                          
                       end         
                    else
                        alpha = tmp;
                        for uni = 1:length(alpha)
                            pos = find(alpha(uni) == validIdx);
                            if length(pos)> 0
                                useful = [useful pos(1)];
                            end
                        end
                        
                    end
                    
                    if isempty(useful)
                        disp('isempty(useful)')
                        continue;
                    end
                    
                    
                end
             end
        end
        
        
    end
end

