% CellColorMap.m
% Developed by Alex Koch, alexrkoch@gmail.com

% This creates cells that can be colored by a map value, but don't need
% to occur in a regular grid. 

% Custom functions called by this script;
%   colorCell.m

clear;
% Decide what to plot. 1 on, 0 off.
plotKfn = 0;
plotQuad = 0;
plotNTG = 1;
plotWidth = 0;
plotThick = 0;
printFig = 1; % to save the figure or not.
singleWindow = 0;

figTitle = 'cc5-quadstd-1610wind-50x804cells';


%% results
results = []; % Paste results for the property to map
    % format: [lateral window, vertical window, property value]
%%
		
% Bulk shift of Y values. I used this because there was blank space
% from 0 to 20 m on the Y axis. Comment out if not necessary
for e = 1:length(results) 
resultsYfix = results(e,2)-20; 
results(e,2) = resultsYfix;
end


% Cell thickness needs to match whatever increment your Y-values
% changed by in your analysis.
if plotKfn == 1
cellThickness = (10); 
    % For K-function, use 2x whatever the K-scale increment is. I have
    % consistently used 5, so I use 10 here.
end

if plotNTG == 1 || plotQuad == 1 || plotWidth == 1 || plotThick == 1
    cellThickness = 15; % vertical window movement increment
end
% 
if singleWindow == 1 % single window would plot one lateral window, and
        % all of its vertical results as just one column of cells.
    limX = 2;
else
    limX = max(results(:,1))+1; % the X axis is just going to be the 
        % number of lateral windows there are. adding 1 to the end
        % because thats where the final cell will end.
end
limY = max(results(:,2))+cellThickness; %limY and minY are used to 
    % constrain the figure axes later.
minY = min(results(:,2))-cellThickness;

% cells is a variable to place the cell objects in
cells = colorCell();

% Construct the cells
if  singleWindow == 1
    for a = 1:length(results)
        cells(a) = colorCell(results(a,1), results(a,2), 3, cellThickness, results(a,3));
    end
else
    for a = 1:length(results)
    cells(a) = colorCell(results(a,1), results(a,2), 1, cellThickness, results(a,3));
    end
end

%% For Cluster
if plotKfn == 1 || plotQuad == 1
% Plot and color cells
figure;
for b=1:length(cells)
    if cells(b).z > 0 % if clustered make it black
        rectangle('Position', [cells(b).x, cells(b).y, cells(b).width, ...
            cells(b).height], 'FaceColor', [0 0 0], 'LineStyle', 'none');
    elseif cells(b).z < 0 % if uniform make it gray
        rectangle('Position', [cells(b).x, cells(b).y, cells(b).width, ...
            cells(b).height], 'FaceColor', [0.7 0.7 0.7], 'LineStyle', ...
            'none');
    end
end
xlim([1 limX]);
% xticklabels('auto');
% xticks(1:1:(limX-1));
ylim([minY limY]);

if printFig == 1
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 5 2.1]; % this controls the dimensions of
        % figure file. It took a lot of trial and error for me to 
        % find what worked for certain figures to make them all uniform
    daspect([0.18 1 1]); % this fixes the aspect ratio of the figure.
        % combined with fig.PaperPosition above, you'll just have to 
        % experiment until the figure comes out how you want it. Then 
        % make sure to document the settings you used for each type of
        % figure.

end
end

%% For Geo Stats
if plotNTG == 1 || plotWidth == 1 || plotThick == 1 
    figure;
for b=1:length(cells)
    c = 1-((cells(b).z-(min(results(:,3))))/ ( 833.8 - ... 
        min(results(:,3)))); %normalizes colorscale to 1 = blk, 0 = 
        % white. replace 833.8 with this for regular range: 
        % (max(results(:,3))
    rectangle('Position', [cells(b).x, cells(b).y, cells(b).width, ...
            cells(b).height], 'FaceColor', [c c c], 'LineStyle', 'none');
end
if printFig == 1
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 5 2.5]; % See comments above.
    daspect([0.132 1 1 ]); 
        
    xlim([0 limX])
    if singleWindow == 1
        xticks([0 2]);
        xticklabels({' ', ' '});
    end

    ylim([0 limY])
    set(gca,'fontsize',10);
    print(figTitle,'-dpng','-r600')
%     print(figTitle,'-dtiff','-r1200')
end
end



