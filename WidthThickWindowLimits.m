% WidthThickWindowLimits.m
% Developed by Alex Koch, contact at alexrkoch@gmail.com

% Calculates the maximum width and thickness of polygons within 
% lateral and vertical moving windows.


clear;

%% Load pointsets
load('CC_Centroids.mat'); %This is a variable saved with all the pointsets 
    %of the project.
points = CC5_Centroids_Flattened; % set desired pointset to variable

load('CC5_Shapefile_mInterp_noNaN.mat');
shp = CC5shp; %load the corresponding shapefile to the pointset above



%% define moving window
movingWindowIncrement = 50; % this is the amount that the moving window
    % jumps laterally after each run. 
movingWindowLength = 775; %input Length of window must be divisible 
    % by cell width (cw)
movingWindowThickness = 100;
xminimum = floor(min(points(:, 1))); 
xmaximum = ceil(max(points(:, 1)));
R = ceil(((xmaximum - xminimum) - movingWindowLength)/ ... 
    movingWindowIncrement);
widths = [];
widthsNoNorm = [];
thicknesses = [];
thicknessesNoNorm = [];
widthToThicks = [];
widthToThicksNoNorm = [];
printFig = 0;
wtVstrat = 0;
printWT = 0;
saveVars = 1;
outputWidths = [];
outputThicknesses = [];
clipWidths = 0 ;

for iii=0%:R %loop through lateral windows
disp('Lateral Window =')
disp(iii)
winStart = xminimum + (movingWindowIncrement*iii);
winEnd = (xminimum + movingWindowLength) + (movingWindowIncrement*iii);
windowPoints = [];  %will be filled with only the points in the window
for i=1:size(points,1) %determines if points are in window
    if points(i,1)>winStart && points(i,1)<winEnd
       newRow = [points(i,1) points(i,2)];
       windowPoints = vertcat(windowPoints, newRow);    
    end
end


%% Calculate WT of vertical window position

base = floor(min(windowPoints(:,2)));
increment = 5;  %INPUT vertical moving window increment
top = (ceil(max(windowPoints(:,2))))-(movingWindowThickness);
stratThickness = max(windowPoints(:,2)) - min(windowPoints(:,2));
limit = top-base ;

for ii = base:increment:top %loop through vertical windows
    disp(['ii = ',num2str(ii),':',num2str(increment),':',num2str(top)])
    widths = [] ;
    thicknesses = [] ;
    widthToThick = 0 ;
    n = 0 ; % number of belts measured 
    partials = [] ; 
    completes = [] ; 
    winTop = ii + movingWindowThickness ; 
    winBase = ii ; 
    e = 0 ; %this will be the index for windowPolys

    for a = 1:length(shp) %evaluate if polygon is in window fully or 
        %partially
        if (shp(a).BoundingBox(1,1)) >= winStart && ...
                (shp(a).BoundingBox(2,1)) <= winEnd && ...
                (shp(a).BoundingBox(1,2)) >= winBase && ...
                (shp(a).BoundingBox(2,2)) <= winTop
            completes = horzcat(completes, a);
            e = e+1;
            windowPolys(e).Index = a ;
            windowPolys(e).X = shp(a).X ;
            windowPolys(e).Y = shp(a).Y ;
            windowPolys(e).MaxThickness = shp(a).MaxThickness;
            thicknesses = vertcat(thicknesses, (shp(a).MaxThickness));
            windowPolys(e).Width = shp(a).Width;
            widths = vertcat(widths, (shp(a).Width));
        else
            allOut = 0; %starts at 1 if each polygon has a NaN 
                %cell at the end of their vertices. I've removed these
                %from some files though.
            for d = 1:length(shp(a).X) % omit those mischevious 
                % polygons which have all their vertices outside the 
                % bounding box.
                if shp(a).X(d) >= winEnd || ...
                    shp(a).X(d) <= winStart || ...
                    shp(a).Y(d) >= winTop || ...
                    shp(a).Y(d) <= winBase 
                    allOut = allOut + 1;
                end
            end
            if allOut ~= length(shp(a).X)
                partials = horzcat(partials, a);
            end
        end
    end
if clipWidths == 1
    for b = 1:length(partials)
    partialPoly = [];
        for c = 1:length(shp(partials(b)).X) %each of the X values 
                % within b^^
            % below: if vertice is fully in window, count it in. 
            if (shp(partials(b)).X(c)) >= winStart && ...
                    (shp(partials(b)).X(c)) <= winEnd && ...
                    (shp(partials(b)).Y(c)) >= winBase && ...
                    (shp(partials(b)).Y(c)) <= winTop
                verticeRow = [(shp(partials(b)).X(c)), ...
                    (shp(partials(b)).Y(c))]; %if the vertice is 
                    % entirely within window, add it to the new XY list 
                    % which will be used to create a clipped polygon.
                partialPoly = vertcat(partialPoly, verticeRow); 
            end % end going through vertices of 'b' partial polygon
        end
        if sum(partialPoly)>0
                partialPolyNew = [partialPoly(1,1),partialPoly(1,2)];
                for p = 2:(length(partialPoly)-1)
                if partialPoly(p,1) ~= partialPoly((p-1),1) || ...
                        partialPoly(p,2) ~= partialPoly((p-1),2)
                    newRow = [partialPoly(p,1),partialPoly(p,2)];
                    partialPolyNew = vertcat(partialPolyNew, newRow);
                end
                end
                partialPolyNew = vertcat(partialPolyNew, ...
                    [partialPoly(1,1),partialPoly(1,2)]);
                e = e+1 ;
                windowPolys(e).X = partialPolyNew(:,1);
                windowPolys(e).Y = partialPolyNew(:,2);
                windowPolys(e).Index = partials(b);
                windowPolys(e).MaxThickness = ...
                    shp(partials(b)).MaxThickness;    
                thicknesses = vertcat(thicknesses, ...
                    (shp(partials(b)).MaxThickness));
                windowPolys(e).Width = max(partialPolyNew(:,1)) - ...
                    min(partialPolyNew(:,1));
                widths = vertcat(widths, (max(partialPolyNew(:,1))...
                    - min(partialPolyNew(:,1))));
                n = n + 1 ;
        end
        end % end for b partials
        if printFig == 1
            fig = gcf;
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 3 3];
            print(figSaveTitle,'-dtiff','-r1200')
        end
        aveWidth = mean(widths);
        aveThickness = mean(thicknesses);
        newRowWidth = [iii, ii, aveWidth];
        newRowThickness = [iii, ii, aveThickness];
        outputWidths = vertcat(outputWidths, newRowWidth);
        outputThicknesses = vertcat(outputThicknesses,newRowThickness);
        
else %don't clip, just use whole widths.
     for b = 1:length(partials)
        e = e+1;
        windowPolys(e).Index = a ;
        windowPolys(e).X = shp(partials(b)).X ;
        windowPolys(e).Y = shp(partials(b)).Y ;
        windowPolys(e).MaxThickness = shp(partials(b)).MaxThickness;
        thicknesses = vertcat(thicknesses, ...
            (shp(partials(b)).MaxThickness));
        windowPolys(e).Width = shp(partials(b)).Width;
        widths = vertcat(widths, (shp(partials(b)).Width));
        n = n + 1 ;
    end % end for b partials
        if printFig == 1
            fig = gcf;
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 3 3];
            print(figSaveTitle,'-dtiff','-r1200')
        end
        aveWidth = mean(widths);
        aveThickness = mean(thicknesses);
        newRowWidth = [iii, ii, aveWidth];
        newRowThickness = [iii, ii, aveThickness];
        outputWidths = vertcat(outputWidths, newRowWidth);
        outputThicknesses = vertcat(outputThicknesses,newRowThickness);

end %end clipWidths
    end %end vertical window    
   
    
% Plot current window through stratigraphy
    if wtVstrat == 1
        figure;
        hold on;
        for e = 1:length(widths)
            if widths(e,1) == iii
                scatter(widths(e,2), widths(e,3), 50, 'k.');
            end
        end
        xlabel('WT');
        ylabel('Strat. Level (m)');
        title(figTitle)
        if printWT == 1
            fig = gcf;
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 3 3];
            print(figWTSave,'-dtiff','-r1200')
        end
        hold off;
    end

end %end lateral window
%%
saveName = ['cc5-wt-',num2str(clipWidths),'-',...
    num2str(movingWindowLength), 'wide-',...
    num2str(movingWindowThickness),'thick-',...
    num2str(movingWindowIncrement), 'latinc-',num2str(increment),...
    'vertinc'] ;
if saveVars == 1
    save(saveName);
end



    
    
    
    
    
    
    
