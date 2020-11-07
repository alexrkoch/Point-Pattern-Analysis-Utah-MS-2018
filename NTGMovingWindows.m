% NTGMovingWindows.m
% Developed by Alex Koch, contact at alexrkoch@gmail.com
% Repository at https://github.com/alexrkoch

% This script utilizes shapefiles of channel belt polygons, and
% calculates the net-to-gross (NTG), which is net polygon to gross data
% area (representing net sand to gross rock area). It calculates NTG
% over lateral and vertical moving windows. It also can calculate the 
% 1D NTG of a vertical line in the center of any given window. This is 
% to simulate the NTG that a single well would encounter, and compare
% it to the 2D NTG in the surrounding subsurface. 

clear;

%% Load pointsets
load('CC_Centroids.mat') ; %This is a variable saved with all the 
    % pointsets of the project.

points = CC5_Centroids_Flattened ; %set desired pointset to variable 

load('CC5_Shapefile_mInterp_noNaN.mat');
shp = CC5shp; %load the corresponding shapefile to the pointset above


plot2Dvs1D = 0 ;
plotCorrelations = 0 ;
plotHist = 0;
correlations = [] ;

% figure of all the polygons
% figure;
% hold on;
% xticks([347:1:348])
% for a = 1:length(shp)
%     fill(shp(a).X, shp(a).Y, [0 0 0]) ;
% end

%% define moving window

movingWindowIncrement = 50; % this is the amount that the moving window
    % jumps laterally after each run. 
movingWindowLength = 600 ;  %The way the code is written currently is  
% meant to handle even values for window length only. 
movingWindowThickness = 15 ; %match this with the frame thickness of  
    % the quadrat (if you intend to do quadrat-NTG correlation)
xminimum = floor(min(points(:, 1))) ; 
xmaximum = ceil(max(points(:, 1))) ;
R = ceil(((xmaximum - xminimum) - movingWindowLength)/ ...
    movingWindowIncrement);
gross = movingWindowThickness * movingWindowLength ;
increment = 15;  %INPUT vertical moving window increment
saveName = ['cc5-1D-2Dntg-',num2str(movingWindowLength), 'wide-',...
    num2str(movingWindowThickness),'thick-',...
    num2str(movingWindowIncrement)] ;
figTitle = ['cc5-1D-2Dntg-',num2str(movingWindowLength), 'wide-',...
    num2str(movingWindowThickness),'thick-',...
    num2str(movingWindowIncrement),'Correlation_Figure'] ;

% Control Panel
ntgVstrat = 0 ;
printNTG = 0 ;
saveVars = 1 ;
ntgs2DAll = [];
   
for iii = 0:R
disp('Lateral Window =')
disp(iii)
lineX = xminimum + (iii*movingWindowIncrement) +(movingWindowLength/2);
winStart = xminimum + (movingWindowIncrement*iii) ;
winEnd = (xminimum + movingWindowLength) + (movingWindowIncrement*iii);
windowPoints = [] ;  
ntgs2D = [] ;
ntgs1D = [] ;
for i=1:size(points,1) %determines if points are in window
    if points(i,1)>winStart && points(i,1)<winEnd
       newRow = [points(i,1) points(i,2)] ;
       windowPoints = vertcat(windowPoints, newRow) ; 
    end
end


%% Calculate NTG of vertical window position
base = floor(min(windowPoints(:,2))) ;
top = (ceil(max(windowPoints(:,2))))-(movingWindowThickness) ;
stratThickness = max(windowPoints(:,2)) - min(windowPoints(:,2)) ;
limit = top-base ; 

for ii = base:increment:top %loop through vertical windows
nets = 0 ; 
partials = [] ; 
completes = [] ; 
winTop = ii + movingWindowThickness ; 
winBase = ii ; 


for a = 1:length(shp) %evaluate if polygon is in window fully or 
        % partially
if (shp(a).BoundingBox(1,1)) >= winStart && ...
        (shp(a).BoundingBox(2,1)) <= winEnd && ...
        (shp(a).BoundingBox(1,2)) >= winBase && ...
        (shp(a).BoundingBox(2,2)) <= winTop;
    nets = nets + shp(a).SHAPE_Area; %if bounding box is entirely 
        % in the window, add it to the net.
    completes = horzcat(completes, a);

else
    allOut = 1; %starts at 1 because each polygon has a NaN cell at 
        % the end of their vertices. 
    for d = 1:length(shp(a).X) % omit those mischevious polygons 
            % which have all vertices out side the bounding box.
        if ((shp(a).X(d) < winStart && shp(a).Y(d) > winTop) ... 
                || (shp(a).X(d) > winStart && shp(a).X(d) < ...
                    winEnd && shp(a).Y(d) > winTop) ...
                || (shp(a).X(d) > winEnd && shp(a).Y(d) > winTop) ... 
                || (shp(a).X(d) > winEnd && shp(a).Y(d) < winTop &&...
                    shp(a).Y(d) > winBase) ...
                || (shp(a).X(d) > winEnd && shp(a).Y(d) < winBase) ... 
                || (shp(a).X(d) > winStart && shp(a).X(d) < winEnd ...
                    && shp(a).Y(d) < winBase) ...
                || (shp(a).X(d) < winStart && shp(a).Y(d) < winBase)... 
                || (shp(a).X(d) < winStart && shp(a).Y(d) < winTop...
                    && shp(a).Y(d) > winBase)); 
            allOut = allOut + 1;
        end
    end
    if allOut ~= length(shp(a).X)
        partials = horzcat(partials, a);
    end
end
end

for b = 1:length(partials) %crops each polygon that has any vertices 
        % outside of the window, adds its area to net
partialPoly = [];
for c = 1:length(shp(partials(b)).X) %each of the X values w/in b^^
if (shp(partials(b)).X(c)) >= winStart && (shp(partials(b)).X(c)) ...
        <= winEnd && (shp(partials(b)).Y(c)) >= winBase && ...
        (shp(partials(b)).Y(c)) <= winTop
    verticeRow = [(shp(partials(b)).X(c)) (shp(partials(b)).Y(c))]; 
        %if the vertice is entirely within window, add it to the new 
        %XY list which will be used to create a clipped polygon.

elseif (shp(partials(b)).X(c)) < winStart && (shp(partials(b)).Y(c))...
        >= winBase && (shp(partials(b)).Y(c)) <= winTop
    verticeRow = [winStart, (shp(partials(b)).Y(c))]; %if out on North

elseif (shp(partials(b)).X(c)) > winEnd && (shp(partials(b)).Y(c))...
        >= winBase && (shp(partials(b)).Y(c)) <= winTop
    verticeRow = [winEnd, (shp(partials(b)).Y(c))]; %if out on South

elseif (shp(partials(b)).X(c)) >= winStart && ...
        (shp(partials(b)).X(c)) <= winEnd && (shp(partials(b)).Y(c))...
        > winTop
    verticeRow = [(shp(partials(b)).X(c)), winTop]; %if out on Top

elseif (shp(partials(b)).X(c)) >= winStart && ...
        (shp(partials(b)).X(c)) <= winEnd && (shp(partials(b)).Y(c))...
        < winBase
    verticeRow = [(shp(partials(b)).X(c)), winBase]; %if out on Base

elseif (shp(partials(b)).X(c)) < winStart && (shp(partials(b)).Y(c))...
        > winTop
    verticeRow = [winStart, winTop]; %if out on North and Top

elseif (shp(partials(b)).X(c)) > winEnd && (shp(partials(b)).Y(c))...
        > winTop
    verticeRow = [winEnd, winTop]; %if out on South and Top

elseif (shp(partials(b)).X(c)) > winEnd && (shp(partials(b)).Y(c))...
        < winBase
    verticeRow = [winEnd, winBase]; %if out on South and Base         

elseif (shp(partials(b)).X(c)) < winStart && (shp(partials(b)).Y(c))...
        < winBase
    verticeRow = [winStart, winBase]; %if out on North and Base
end

partialPoly = vertcat(partialPoly, verticeRow);
end


partialArea = polyarea(partialPoly(:,1), partialPoly(:,2)) ;
nets = nets + partialArea ;
end

ntg = nets/gross ;
iiNorm = (ii-base)*(stratThickness/limit) ; %normalize window position 
    % to stratigraphic level
newRow = [iii, ii, ntg] ;
ntgs2D = vertcat(ntgs2D, newRow);
end %end vertical window - 2D area
 
% Keep running list of all lateral window results
ntgs2DAll = vertcat(ntgs2DAll, ntgs2D);

%% Loop through each polygon and perform the workflow

shpIntersect = []; %store list of polys that intersect our 1D line
for a = 1:length(shp) % search each polygon
    if (shp(a).BoundingBox(1,1)) < lineX &&(shp(a).BoundingBox(2,1))...
            > lineX
        shpIntersect = vertcat(shpIntersect, a);
    end
end

% If Index exceeds matrix dimensions, and no intersecting vertices are
% found for a shape that does intersect the lineX, check the value of
% lineX. The way the code is written currently is meant to handle even
% values for window length only. 
for a = 1:length(shpIntersect) % search each polygon
    yTemp = [];
    for b = 1:length(shp(shpIntersect(a)).X)
        if shp(shpIntersect(a)).X(b) == lineX
            yTemp = vertcat(yTemp, shp(shpIntersect(a)).Y(b));
        end
    end
    yTempSort = sort(yTemp);
    struct1D(a).base = yTempSort(1);
    struct1D(a).top = yTempSort(2);
end
    
% now go through the vert windows and tabulate net.

for ii = base:increment:top
    nets = [];
    for a = 1:length(struct1D)
        if struct1D(a).base >= ii && struct1D(a).top <= (ii + ...
                movingWindowThickness)
            nets = vertcat(nets, (struct1D(a).top - struct1D(a).base));
                % if poly is fully in window
        elseif struct1D(a).base < ii && struct1D(a).top <= (ii + ...
                movingWindowThickness) && struct1D(a).top > ii
            nets = vertcat(nets, (struct1D(a).top - ii));
                % if top is in but base is below window
        elseif struct1D(a).base >= ii && struct1D(a).base < (ii + ...
                movingWindowThickness) && struct1D(a).top > (ii +  ...
                movingWindowThickness)
            nets = vertcat(nets, ((ii + movingWindowThickness) - ...
                struct1D(a).base));
                % if base is in but top is above window
        elseif struct1D(a).base < ii &&  struct1D(a).top > (ii + ...
                movingWindowThickness)
            nets = vertcat(nets, movingWindowThickness);
                % if base is in but top is above window
        end
    end
    net = sum(nets);
    ntg = net/movingWindowThickness ;
    iiNorm = (ii-base)*(stratThickness/limit) ; %normalize window 
        % position to stratigraphic level
    newRow = [iii, iiNorm, ntg] ;
    ntgs1D = vertcat(ntgs1D, newRow);
end

%% Correlation Coefficient

Correlations1D2D = corrcoef(ntgs2D(:,3), ntgs1D(:,3));
Corr1D2D = Correlations1D2D(1,2);
newRow = [iii,Corr1D2D];
correlations = vertcat(correlations, newRow);

%% Plots
if plot2Dvs1D == 1
    figure;
    hold on;
    plot(ntgs2D(:,3), ntgs2D(:,2));
    plot(ntgs1D(:,3), ntgs1D(:,2), '--')
    title([num2str(winStart), ' to ', num2str(winEnd)]);
    xlabel('NTG')
    ylabel('Stratigraphic Level (m)')
end

end % end lateral window


%% increase iii by 1
% for e = 1:length(correlations)
%     correlations(e,1) = correlations(e,1)+1;
% end


%%

% if plotCorrelations == 1
    figure;
    scatter(correlations(:,1),correlations(:,2), 75, 'k.');
    ylim([0 1])
    xlim([0 37])
%     title('Correlation of 1D to 2D NTG','FontSize',10)
%     ylabel('Correlation Coefficient','FontSize',10)
%     xlabel('Lateral Window','FontSize',10)
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 2 2];
    print('CC5-1D-2Dntg-600wide-25thick-50inc-sctr','-dtiff','-r1200');  
        %print saves this out to file
% end
%     if plotHist == 1

    figure;
    histogram(correlations(:,2), 'FaceColor', 'k');
%         title([OC,' Cluster Index to NTG']);
%         xlabel('r^2 value');
%         ylabel('Frequency');
    xlim([0 1])
%         if printHist == 1
        fig = gcf;
        fig.PaperUnits = 'inches';
        fig.PaperPosition = [0 0 2 2];
        print('CC5-1D-2Dntg-600wide-25thick-50inc-hst',...
            '-dtiff','-r1200');
%         end
%     end
aveCorrCoeff = mean(correlations(:,2));
%%
saveName = ['cc5-1D-2Dntg-',num2str(movingWindowLength), 'wide-',...
    num2str(movingWindowThickness),'thick-',...
    num2str(movingWindowIncrement)] ;
if saveVars == 1
    save(saveName); % 'correlations','ntgs1D','ntgs2D');
end

    
    
    
