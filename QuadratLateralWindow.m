% QuadratLateralWindow.m
% Developed by Alex Koch, contact at alexrkoch@gmail.com
% Modified from Benhallam, W. 2015, "Spatial analysis of channel-belt
% stacking patterns: Metrics to discriminate between local and regional
% controls on deposition in the fluvial John Henry Member of the 
% Straight Cliffs Formation, southern Utah" Appendix C.

% Performs a vertical moving window quadrat point pattern analysis, 
% over any set of lateral moving windows.

% Custom functions called by this script;
%   Cell.m

clear;

%% Load pointsets
load('CC_Centroids.mat'); %This is a variable saved with all the 
    % pointsets of the project.
points = CC5_Centroids_Flattened; % here I specify which to use
load('CC5_Shapefile_mInterp_noNaN.mat');
shp = CC5shp; % pick the shapefile that matches your point set

%input all the moving window lengths you want to test.
movingWindowLengths = [1610]; 
for b = 1:length(movingWindowLengths)

%% Define moving window, set up empty arrays
movingWindowIncrement = 50; % this is the amount that the moving window
    % jumps laterally after each run. 
movingWindowLength = movingWindowLengths(b); %input Length of window 
    % must be divisible by cell width (cw)
xminimum = floor(min(points(:, 1)));
xmaximum = ceil(max(points(:, 1)));
R = ceil(((xmaximum - xminimum) - movingWindowLength)/...
    movingWindowIncrement);

aveWidthRatios = [];

chs = [50];
cwGuide =[805];
cws = [];

%% determine cell widths that will be divisors of the current 
% movingWindowLength
for c = 1:length(cwGuide)
    a = movingWindowLength / cwGuide(c);
    d = floor((movingWindowLength / round(a))*100)/100;
    if (movingWindowLength / d) >=2 % ensure that for small windows 
            % lengths, at least 2 cells will fit widthwise.
        cws(c) = d;
    end
end

for ix = 1:length(cws) %loop through different cell widths
% for ixx = 1:length(chs) %loop through different cell thicknesses

ratiosChecks = [];
allRatios = [];
allRatiosNoNorm = [];
allP90s = [];
allP10s = [];
maxminBreech = [];
allStats = [];
    
for iii= 0:1:R
%Adjust the limit number based on window size and increment

%% Pull the window subset from points.
winStart = xminimum + (movingWindowIncrement*iii);  %Gets minimum X 
    % value of window
winEnd = (xminimum + movingWindowLength) + (movingWindowIncrement*iii);  
    % Gets maximum X value of window
windowPoints = [];  %this will be filled with only the points in the 
    % window
for i=1:size(points,1) %determines if points are in window
    if points(i,1)>winStart && points(i,1)<winEnd
       newRow = [points(i,1) points(i,2)];
       windowPoints = vertcat(windowPoints, newRow);%adds points to 
            % windowPoints.
    end
end

%% inputs and variables

clear scatter ;
%INPUTS 1 is on, 0 is off
runData = 1 ;   %runs actual dataset
runSims = 1 ;   %runs and plots random simulations
plotFig = 0 ;   %creates Quadrat figure for individual window
printFig = 0 ;  %saves Quadrat figure to file
movingAverage = 0 ; %calculates moving average of the cluster index 
    % for individual window plot
plotGrid = 0 ;  %plots poinset with last cell overlay
printGrid = 0 ; %saves pointset and cell overlay to file
% The code has options to plot high and low percentiles of simulation
% results, in addition to max and min (as in the K-function). Pick 
% below which you want to plot (can plot multiple).
plotP60and40 = 0 ; 
plotP70and30 = 0 ;
plotP80and20 = 0 ;
plotP90and10 = 1 ;
plotMaxAndMin = 1 ;
calculateBreeches = 1; %calculates how far the cluster index lies 
    % outside the simulation value
plotBreeches = 0 ;
printBreeches = 0 ;
averageWidths = 0 ; %averages all cell widths for each cell 
    % thickness value
saveVars = 1 ; % do you want to save the variables?
windowLength = winEnd - winStart;
cw = cws(ix) ; %INPUT cell width must be a divisor of windowLength
ch = chs(ix) ; %INPUT cell height must be < total strat thickness / 2. 

% display progress of the code during run
disp(['Window Length =  ',num2str(movingWindowLength)])
disp(['Cell Dimensions = ',num2str(cw),' x ', num2str(ch)])
disp(['Lateral Window ',num2str(iii)])

rows = 2;  %INPUT number of cell rows (can be 1, 2, or 3)
increment = 5;  %INPUT increment of vertical moving window change.
base = floor(min(windowPoints(:,2)));
top = (ceil(max(windowPoints(:,2))))-(ch*rows);
ratios = [];
ratiosNoNorm = []; % holds results without normalization.
ratiosWLatWind = [];
MonteCarloSimulations = 99 ; %INPUT number of random simulations
limit = top - base ;
stratThickness = max(windowPoints(:,2)) - min(windowPoints(:,2));

%INPUT name the title of the Quadrat and Grid figures to save as.
figTitle = ['CC5F-fullQuad-',num2str(cw*100),'x',...
    num2str(ch*100),'-plot-', num2str(winStart)];
gridTitle = ['CC5F-fullQuad-',num2str(cw*100),'x',...
    num2str(ch*100),'-grid-', num2str(winStart)];
totalFigTitle = ['CC51-fullQuad-',num2str(cw*100),'x',...
    num2str(ch*100),'-bar-', num2str(winStart)];

%% Run quadrat for the real data
if runData == 1
%%
for k=base:increment:top
    cells = Cell();
    %Constructing the grid
    %Cell xy coordinates define the lower left corner
    walker = 1;
    for i=winStart:cw:(winEnd-cw)
        cells(walker) = Cell(i, k, cw, ch);  %constructs one cell
        if rows == 1 
            walker = walker + 1;    %to use only one row of cells
        elseif rows == 2 
            cells(walker + 1) = Cell(i, k+ch, cw, ch);  %constructs a 
                % second cell on top of the first.
            walker = walker + 2; %sets up walker for the next pair of 
                % cells.  
        elseif rows == 3 
            cells(walker + 1) = Cell(i, k+ch, cw, ch);  %constructs a 
                % second cell on top of the first.
            cells(walker + 2) = Cell(i, k+(ch*2), cw, ch);  %constructs 
                % a second cell on top of the first.
            walker = walker + 3;  %sets up walker for the next pair of 
                % cells.  
        end
    end

    %Counting the points in the grid
    for i=1:length(windowPoints)
       ni = floor(windowPoints(i,1)-winStart);  %for small values, 
            % large windows x = 0...
       if ni == windowLength %When a point is directly on the edge of 
            % the window, this loop puts it inside.
           ni = ni - 1 ;
       end
       n = floor(ni/cw);
       m = (n*rows)+1 ;     %needs to result in a positive integer.  
       %keep the following, to find if in Y range
       if rows == 1 
         if(windowPoints(i,2) < k+ch && windowPoints(i,2) >= k) %if 
                % Yvalue is within cell then...
           index = m;  %index needs to be the cell index that the point 
                % is within. 
           cells(index) = cells(index).incrementCount();  %add 1 to the 
                % count of that cell.
         end
       end
       if rows == 2 
         if(windowPoints(i,2) < k+ch && windowPoints(i,2) >= k) %if 
                % Yvalue is within cell then...
           index = m;  %index needs to be the cell index that the point 
                % is within. 
           cells(index) = cells(index).incrementCount(); %add 1 to the 
                % count of that cell.
         elseif (windowPoints(i,2) < k+(2*ch) && windowPoints(i,2) ...
                 >= k+ch) %if Yvalue is within 2nd story cell then...
           index = m+1;
           cells(index) = cells(index).incrementCount(); %add 1 to the 
                % count for that cell.
         end
       end
       if rows == 3 
         if(windowPoints(i,2) < k+ch && windowPoints(i,2) >= k) %if 
                % Yvalue is within cell then...
           index = m; %index needs to be the cell index that the point 
                % is within. 
           cells(index) = cells(index).incrementCount(); %add 1 to the 
                % count of that cell.
         elseif (windowPoints(i,2) < k+(2*ch) && windowPoints(i,2) ...
                 >= k+ch) %if Yvalue is within 2nd story cell then...
           index = m+1;    
           cells(index) = cells(index).incrementCount(); %add 1 to the 
                % count for that cell.  
         elseif (windowPoints(i,2) < k+(3*ch) && windowPoints(i,2) ...
                 >= k+ch) %if Yvalue is within 3rd story cell then...
           index = m+2;   
           cells(index) = cells(index).incrementCount();%add 1 to the 
                % count for that cell.
         end
       end
    end
    %Counting the total number of points in the grid
    counts = {cells.count};
    pointsIn = 0; 
    for i=1:length(counts)
        pointsIn = pointsIn + counts{i};    % counts the total number 
            % of points within the grid
    end


    %Computing the mean, variance, and ratio of the grid
    mean = pointsIn/length(cells);
    sum = 0;
    for i=1:length(cells)
        sum = sum + (cells(i).count - mean)^2;
    end

    variance = sum/(length(cells) - 1);
    ratio = variance/mean;
    kNorm = (k-base)*(stratThickness/limit) ; %normalize the vertical 
        % moving window positions to strat thickness.
    row = [iii, kNorm, ratio];
    ratiosRow = [ratio, kNorm];
    rowNoNorm = [iii, k, ratio];
    ratiosRowNoNorm = [ratio, k];
    ratios = vertcat(ratios, ratiosRow);
    ratiosNoNorm = vertcat(ratiosNoNorm, ratiosRowNoNorm);
    ratiosWLatWind = vertcat(ratiosWLatWind, row);
    allRatiosNoNorm = vertcat(allRatiosNoNorm, rowNoNorm);
    mean = 1; 
    %calculate all stats: max, min, range, mean, std for each quad
    %window
    widths = [];
    areas = [];
    clear mean; 
    for f = 1:length(points)
        if (points(f,2) >= k) && (points(f,2) <= (k+increment))...
                && (points(f,1) >= winStart) && (points(f,1) <= winEnd)
            widths = vertcat(widths, shp(f).Width);
%                     areas = []; vertcat(areas, shp(f).SHAPE_Area);
        end    
    end  
    Wmin = min(widths);
    Wmean = mean(widths);
    Wmax = max(widths);
    Wrng = Wmax - Wmin;
    Wstd = std(widths);
    TF = isnan(Wmean);
    if TF == 0
        newRow = [iii, k, Wmin, Wmean, Wmax, Wrng, Wstd];
        allStats = vertcat(allStats, newRow); 
    end

end % end vertical window
allRatios = vertcat(allRatios, ratiosWLatWind); %Master array to 
    % compare against GeologicStats

%% Perform moving average
if movingAverage == 1
movingAverage = [];
aveSize = 5;
for i=aveSize+1:aveSize:length(ratios)
    %i
    average = 0;
    height = 0;
    for j=i-aveSize:i-1
        average = ratios(j,1) + average;
        height = ratios(j,2) + height;
    end
    average = average/aveSize;
    height = height/aveSize;
    newRow = [average height];
    movingAverage = vertcat(movingAverage, newRow);
end
end
end

%% Create simulations, run quadrat
if runSims == 1
allResults = [];
N = size(windowPoints, 1);

for ii=1:MonteCarloSimulations
    randPoints = rand(N,2); %makes an Nx2 martix of random numbers 
        % between 0 and 1
    randPoints(:,1) = winStart + windowLength.*randPoints(:,1); 
        % multiplies width by each of the x Rands (from 0 to 1)
    randPoints(:,2) = base + ((ceil(max(windowPoints(:,2))))...
        -base).*randPoints(:,2); %same as ^^ but for y's.
    randRatios = [];

    for k=base:increment:top
        cells = Cell();
        %Constructing the grid
        %Cell xy coordinates define the lower left corner
        walker = 1;
        for i=winStart:cw:(winEnd-cw)
            cells(walker) = Cell(i, k, cw, ch);  %constructs one cell
            if rows == 1 ;
                walker = walker + 1;    %to use only one row of cells
            elseif rows == 2 ;
                cells(walker + 1) = Cell(i, k+ch, cw, ch);  %constructs
                    % a second cell on top of the first.
                walker = walker + 2;    %sets up walker for the next 
                    % pair of cells.  
            elseif rows == 3 ;
                cells(walker + 1) = Cell(i, k+ch, cw, ch);  %constructs 
                    % a second cell on top of the first.
                cells(walker + 2) = Cell(i, k+(ch*2), cw, ch); 
                    % constructs a second cell on top of the first.
                walker = walker + 3; %sets up walker for the next pair 
                    % of cells.  
            end
        end

        %Counting the points in the grid
        %If you get an error like 'exceeds index' in this counting
        %process, and you find that the index and m values are greater
        %than the greatest cell index, it might be because your cells
        %are too small, and there is a gap at the southern end, where
        %there is window space (a point may be within the window) but
        %not inside any cell. this causes this error, so then tighten
        %up your cells so they are just right against the window.
        for i=1:length(randPoints)
           ni = floor(randPoints(i,1)-winStart);  %for small values, 
                % large windows x = 0...
           if ni == windowLength %When a point is directly on the edge 
                % of the window, this loop puts it inside.
               ni = ni - 1 ;
           end
           n = floor(ni/cw);
           m = (n*rows)+1 ;     %needs to result in a positive integer
           %keep the following, to find if in Y range
           if rows == 1 
             if(randPoints(i,2) < k+ch && randPoints(i,2) >= k) 
                    %if Yvalue is within cell then...
               index = m;
               cells(index) = cells(index).incrementCount();  
                    %...add 1 to the count of that cell.
             end
           end
           if rows == 2 
             if(randPoints(i,2) < k+ch && randPoints(i,2) >= k) 
                    %if Yvalue is within cell then...
               index = m;
               cells(index) = cells(index).incrementCount();  
                    %...add 1 to the count of that cell.
             elseif(randPoints(i,2) < k+(rows*ch) && randPoints(i,2)...
                     >= k+ch) %if Yvalue is within 2nd story cell then
               index = m+1;
               cells(index) = cells(index).incrementCount();  
                    %...add 1 to the count for that cell.
             end
           end
           if rows == 3 
             if(randPoints(i,2) < k+ch && randPoints(i,2) >= k) 
                    %if Yvalue is within cell then...
               index = m;   %index needs to be the cell index that the 
                    % point is within. 
               cells(index) = cells(index).incrementCount(); %...add 
                    % 1 to the count of that cell.
             elseif (randPoints(i,2) < k+(2*ch) && randPoints(i,2) ...
                     >= k+ch) %if Yvalue is within 2nd story cell then
               index = m+1;    
               cells(index) = cells(index).incrementCount();
                    %...add 1 to the count for that cell.  
             elseif (randPoints(i,2) < k+(3*ch) && randPoints(i,2) ...
                     >= k+ch) %if Yvalue is within 3rd story cell then
               index = m+2;   
               cells(index) = cells(index).incrementCount();                 
                    %...add 1 to the count for that cell.
             end
           end
         end

        % Counting the total number of points in the grid
        counts = {cells.count};
        pointsIn = 0; 
        for i=1:length(counts)
            pointsIn = pointsIn + counts{i};
        end

        % Computing the mean, variance, and ratio of the grid
        mean = pointsIn/length(cells);
        sum = 0;
        for i=1:length(cells)
            sum = sum + (cells(i).count - mean)^2;
        end

        variance = sum/(length(cells) - 1);
        randRatio = variance/mean;
        randRatios = [randRatios ; randRatio] ;

    end

    allResults = horzcat(allResults, randRatios);

end
simulations = allResults;

%% Get max and min and percentile values of simulated cluster indices
    maxs = [];
    mins = [];
    for i=1:size(simulations, 1) %goes through each 'window' and takes 
            % the max and min out of any simulation.
        maxs = [maxs; max(simulations(i,:))];
        mins = [mins; min(simulations(i,:))];
    end

    % Get percentiles of simulations
    if plotP60and40 == 1
        P60s = quantile(simulations, 0.6, 2);
        P40s = quantile(simulations, 0.4, 2);
    end
    if plotP70and30 == 1
        P70s = quantile(simulations, 0.7, 2);
        P30s = quantile(simulations, 0.3, 2);
    end
    if plotP80and20 == 1
        P80s = quantile(simulations, 0.8, 2);
        P20s = quantile(simulations, 0.2, 2);
    end
    if plotP90and10 == 1
        P90s = quantile(simulations, 0.9, 2);
        P10s = quantile(simulations, 0.1, 2);
    end        

end

% Find where the pattern is greater/less than p90/p10 of simulations.
if calculateBreeches == 1
for t=1:size(ratiosNoNorm, 1)
    ratiosCheck = []; %will hold results for just this lat window.
    if ratiosNoNorm(t,1) > P90s(t,1)
        ratiosCheck = [(iii), (ratiosNoNorm(t,1)-P90s(t,1)), ...
            ratiosNoNorm(t,2), ratiosNoNorm(t,1)];
        ratiosChecks = vertcat(ratiosChecks, ratiosCheck);
    elseif ratios(t,1) < P10s(t,1)
        ratiosCheck = [(iii), (ratiosNoNorm(t,1)-P10s(t,1)), ...
            ratiosNoNorm(t,2),ratiosNoNorm(t,1)];
        ratiosChecks = vertcat(ratiosChecks, ratiosCheck);
    else
        ratiosCheck = [(iii), 0, ratiosNoNorm(t,2),ratiosNoNorm(t,1)];
        ratiosChecks = vertcat(ratiosChecks, ratiosCheck);
    end
end
end

if calculateBreeches == 1
for t=1:size(ratiosNoNorm, 1)
    if ratiosNoNorm(t,1) > maxs(t,1)
        maxminBreech = vertcat(maxminBreech, ...
            (ratiosNoNorm(t,1)-maxs(t,1)));
    elseif ratios(t,1) < mins(t,1)
        maxminBreech = vertcat(maxminBreech, ...
            (ratiosNoNorm(t,1)-mins(t,1)));
    else
        maxminBreech = vertcat(maxminBreech, 0);
    end
end
end

%% Plot figures
figTitle = ['CC5Quadrat Cluster Indices, Window ',num2str(winStart),...
    ' m to ', num2str(winEnd), ' m'];
figTitle = ['CC5 Quadrat Window-',num2str(iii)];
gridTitle = 'CC5 Points_1';
totalFigTitle = 'cc1 full quadrat 90 x 60 m cells';

%Plot quadrat results and simulations
if plotFig == 1
    figure;
    scatter(ratios(:,1), ratios(:,2), 50, 'k.'); %All cluster indices

    %Plot simulated cluster indices.
    hold on;
%         plot(movingAverage(:,1), movingAverage(:,2),'k', ...
          % 'LineWidth', 2);   %real data line
    if plotMaxAndMin == 1
        plot(maxs, ratios(:,2), 'k', 'linewidth', 0.5);
        plot(mins, ratios(:,2), 'k', 'linewidth', 0.5);
    end
    if plotP60and40 == 1
        plot(P60s, ratios(:,2), 'k:', 'linewidth',1);
        plot(P40s, ratios(:,2), 'k:', 'linewidth',1);
    end
    if plotP70and30 == 1
        plot(P70s, ratios(:,2), 'k:', 'linewidth',1);
        plot(P30s, ratios(:,2), 'k:', 'linewidth',1);
    end
    if plotP80and20 == 1
        plot(P80s, ratios(:,2), 'k:', 'linewidth',1);
        plot(P20s, ratios(:,2), 'k:', 'linewidth',1);
    end
    if plotP90and10 == 1
        plot(P90s, ratios(:,2), 'k:', 'linewidth',2);
        plot(P10s, ratios(:,2), 'k:', 'linewidth',2);
    end

    %figure formatting
    xlabel('Clustering Index');
    ylabel('Strat. Level (m)');
    title(figTitle);
%         xlim([0 5])
%         ylim([0 120])
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 3 3];

    %save figure
    if printFig == 1
    print(figTitle,'-dtiff','-r1200')
    end

end

%Plot the grid and the centroid points
if plotGrid == 1

        clear scatter ;
        figure; 
        scatter(windowPoints(:,1), windowPoints(:,2), 30, 'k.');
        title(gridTitle);
        for i=1:length(cells)
            rectangle('Position', [cells(i).x, cells(i).y, ...
                cells(i).width, cells(i).height]);
        end
        rectangle('Position', [winStart, base, movingWindowLength,...
            (max(windowPoints(:,2))-min(windowPoints(:,2)))]);   
                %use this for the full vert window

        xlabel('Distance (m)');
        ylabel('Strat. Level (m)');
        if printGrid == 1
            fig = gcf;
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 3.8 2];
            print(gridTitle,'-dtiff','-r1200')
        end
end
hold off;

% Keep master list of the simulations
allP90s = vertcat(allP90s, P90s);
allP10s = vertcat(allP10s, P10s);

end %end lateral window
exportResults = horzcat(ratiosChecks, maxminBreech) ;

    
%%
saveName = ['cc5-quadWithStats-',num2str(movingWindowLength),...
    'wind-',num2str(cw*100),...
    'cw-',num2str(ch*100),'ch-',num2str(movingWindowIncrement),...
    'latinc-',num2str(increment),'vertinc'];
    finalOutput = horzcat(ratiosChecks(:,1), ratiosChecks(:,3), ...
        ratiosChecks(:,2));
if saveVars == 1
    save(saveName);
end

% end %end chs
%%
if averageWidths == 1
    aveWidthRatios = horzcat(aveWidthRatios, ratiosChecks(:,4));
end
end %end cws
end %end movingWindowLengths
%%
if averageWidths == 1
    saveName = ['cc5flat-widthave-',num2str(movingWindowLength),...
        'wind-',num2str(ch*100),...
    'ch-',num2str(movingWindowIncrement),'latinc-',...
    num2str(increment),'vertinc'];
    clear mean;
	meanRatios = mean(aveWidthRatios, 2); %2 indicates the dimension 
        % to calculate across, in this case returns a column vector 
        % with the mean of each row. 
    %finalOutput format is; [lateral window, vertical window, breeches, 
        % width averaged cluster indices];
    finalOutput = horzcat(ratiosChecks(:,1), ratiosChecks(:,3), ...
        ratiosChecks(:,2), maxminBreech);
    save(saveName);
end

    finalOutput = horzcat(ratiosChecks(:,1), ratiosChecks(:,3), ...
        ratiosChecks(:,2), maxminBreech);
    save(saveName);

    

