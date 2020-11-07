% Amalgamation.m
% Developed by Alex Koch, alexrkoch@gmail.com

% Measure amalgamation length and compare it to bottom lengths of 
% channel belts over moving window range. To use this, you need two 
% sets of polylines; the amalgamation contacts, and the channel belt 
% bases. Shapefiles of amalgamation contacts must be meter or cm 
% interpolated (see VerticeInterpolator.m)

clear;

% Load the channel belt base lines, and amalgamation lines
load('CC1Amalg_mInterp.mat');
amalg = CC1Amalg;
load('CC1AmalgBases_mInterp.mat');
bases = CC1AmalgBases;
load('CC_Centroids.mat'); 
points = CC1_Centroids_Flattened;

plotAllLines = 0 ;
finalOutput = []; %stores [iii, ii, avgAmalgIndex];

% Define moving mindows
movingWindowIncrement = 50; % this is the amount that the moving window
    % jumps laterally after each run. 
movingWindowLength = 1000 ; %input Length of window must be divisible 
    % by cell width (cw)
movingWindowThickness = 100 ;
increment = 5 ; % vertical window movement increment
xminimum = floor(min(points(:, 1))) ; 
xmaximum = ceil(max(points(:, 1))) ;
R = ceil(((xmaximum - xminimum) - movingWindowLength)/ ...
    movingWindowIncrement);

saveVars = 1 ;
saveName = ['cc1-ai-',num2str(movingWindowLength),'wide-',...
  num2str(movingWindowThickness),'thick-', ...
  num2str(movingWindowIncrement)] ;

% lateral windows
for iii = 0:R
disp('Lateral Window =')
disp(iii)
winStart = xminimum + (movingWindowIncrement*iii) ;
winEnd = (xminimum + movingWindowLength) + (movingWindowIncrement*iii);
windowPoints = [];
for i=1:size(points,1) %determines if points are in window
    if points(i,1)>winStart && points(i,1)<winEnd
       newRow = [points(i,1) points(i,2)] ;
       windowPoints = vertcat(windowPoints, newRow) ;
    end
end

% vertical windows
base = floor(min(windowPoints(:,2))) ;
top = (ceil(max(windowPoints(:,2))))-(movingWindowThickness) ;
stratThickness = max(windowPoints(:,2)) - min(windowPoints(:,2)) ;
limit = top-base ; 
amalgIncides = [];

for ii = base:increment:top %loop through vertical windows
partials = [] ; 
completes = [] ; 
winTop = ii + movingWindowThickness ; 
winBase = ii ; 
amalgLength = 0;
basesLength = 0.0001;
amalgIndex = 0;

%% Add total amalgamation length
for a = 1:length(amalg) %evaluate if polyline is in window fully or 
        % partially
    if (amalg(a).BoundingBox(1,1)) >= winStart && ...
            (amalg(a).BoundingBox(2,1)) <= winEnd && ...
            (amalg(a).BoundingBox(1,2)) >= winBase && ...
            (amalg(a).BoundingBox(2,2)) <= winTop % if its fully in
        completes = horzcat(completes, a); 
        % Calculate the length of the polyline
        for h = 1:((length(amalg(a).X))-1) % for each vertice except 
                % the last one
            amalgLength = amalgLength + (sqrt(((amalg(a).X(h) - ...
                amalg(a).X(h+1))^2) + ((amalg(a).Y(h) - ...
                amalg(a).Y(h+1))^2)));
                %add the current segment length to the total
                %amalgLength for this shape
        end
    else % see if its all out, if not, then add to partials list
        allOut = 0; %this will tally the number of vertices that are 
            % now in the window.  
        for d = 1:length(amalg(a).X) % omit those mischevious polylines 
                % which have all vertices outside the window
            if amalg(a).X(d) >= winEnd || ...
                amalg(a).X(d) <= winStart || ...
                amalg(a).Y(d) >= winTop || ...
                amalg(a).Y(d) <= winBase 
                allOut = allOut + 1;
            end
        end
        if allOut ~= length(amalg(a).X)
            partials = horzcat(partials, a);
        end
    end % end if poly is in or not
end % end for 'a' polygon

for b = 1:length(partials) %crops each polygon that has any vertices 
        % outside of the window, adds its area to net
    partialPoly = [];
    for c = 1:length(amalg(partials(b)).X) %each of the X vals in b^^
        if (amalg(partials(b)).X(c)) >= winStart && ...
                (amalg(partials(b)).X(c)) <= winEnd && ...
                (amalg(partials(b)).Y(c)) >= winBase && ...
                (amalg(partials(b)).Y(c)) <= winTop
            verticeRow = [(amalg(partials(b)).X(c)) ...
                (amalg(partials(b)).Y(c))]; %if the vertice is entirely 
                %within window, add it to the new XY list which will be 
                % used to create a clipped polygon.
            partialPoly = vertcat(partialPoly, verticeRow);
        end
    end % end c vertice search
    % Calculate the length of the partial polyline
    if length(partialPoly(:,1)) > 1
        for h = 1:(length(partialPoly)-1) % for each vertice except the 
                % last one
            amalgLength = amalgLength + (sqrt(((partialPoly(h,1) - ...
                partialPoly((h+1),1))^2) + ((partialPoly(h,2) - ...
                partialPoly((h+1),2))^2)));
                %add the current segment length to the total
                %amalgLength for this shape
        end
    end

end % end b partial poly

% Add total channel belt base length
partials = [];
completes = [];
for a = 1:length(bases) %evaluate if polygon is in window fully or 
        % partially
    if (bases(a).BoundingBox(1,1)) >= winStart && ...
            (bases(a).BoundingBox(2,1)) <= winEnd && ...
            (bases(a).BoundingBox(1,2)) >= winBase && ...
            (bases(a).BoundingBox(2,2)) <= winTop % if its fully in
        completes = horzcat(completes, a); 
        % Calculate the length of the polyline
        for h = 1:((length(bases(a).X))-1) % for each vertice except 
                % the last
            basesLength = basesLength + (sqrt(((bases(a).X(h) - ...
                bases(a).X(h+1))^2) + ((bases(a).Y(h) - ...
                bases(a).Y(h+1))^2)));
                %add the current segment length to the total
                %basesLength for this shape
        end
    else % see if its all out, if not, then add to partials list
        allOut = 0; %this will tally the number of vertices that are 
            % now in the window.  
        for d = 1:length(bases(a).X) % omit those mischevious polylines 
                % which have all vertices outside the window
            if bases(a).X(d) >= winEnd || ...
                bases(a).X(d) <= winStart || ...
                bases(a).Y(d) >= winTop || ...
                bases(a).Y(d) <= winBase 
                allOut = allOut + 1;
            end
        end
        if allOut ~= length(bases(a).X)
            partials = horzcat(partials, a);
        end
    end % end if poly is in or not
end % end for 'a' polygon

for b = 1:length(partials) %crops each polygon that has any vertices 
        % outside of the window, adds its area to net
    partialPoly = [];
    verticeRow = [];
    for c = 1:length(bases(partials(b)).X) %each of the X vals in b^^
        if (bases(partials(b)).X(c)) >= winStart && ...
                (bases(partials(b)).X(c)) <= winEnd && ...
                (bases(partials(b)).Y(c)) >= winBase && ...
                (bases(partials(b)).Y(c)) <= winTop
            verticeRow = [(bases(partials(b)).X(c)) ...
                (bases(partials(b)).Y(c))]; %if the vertice is entirely 
                %within window, add it to the new XY list which will be 
                % used to create a clipped polygon.
        end
        partialPoly = vertcat(partialPoly, verticeRow);
    end % end c vertice 
    % Calculate the length of the partial polyline
    if length(partialPoly(:,1)) > 1
        for h = 1:(length(partialPoly)-1) % for each vertice except the 
                % last one
            basesLength = basesLength + (sqrt(((partialPoly(h,1) - ...
                partialPoly((h+1),1))^2) + ((partialPoly(h,2) - ...
                partialPoly((h+1),2))^2)));
                %add the current segment length to the total
                %basesLength for this shape
        end
    end
end % end b partial poly

amalgIndex = amalgLength / basesLength ;
newRow = [iii, ii, amalgIndex];
finalOutput = vertcat(finalOutput, newRow);
    
end % end vertical window

end % end lateral window


if saveVars == 1
    save(saveName);
end

%%
if plotAllLines == 1
    figure;
    mapshow(bases, 'Color', 'k')
    mapshow(amalg, 'Color', 'g');
%     rectangle('Position', [winStart, ii, movingWindowLength, ...
        % movingWindowThickness]);   %use this for the full vert window
end












