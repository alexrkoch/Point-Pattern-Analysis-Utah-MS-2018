% KfunctionLatWindVE.m    

% Developed by Alex Koch, contact at alexrkoch@gmail.com
% Modified from Benhallam, W., 2015, "Spatial analysis of channel-belt
% stacking patterns: Metrics to discriminate between local and regional
% controls on deposition in the fluvial John Henry Member of the 
% Straight Cliffs Formation, southern Utah" Appendix C.

% This script performs Ripley's K-function over a lateral moving window 
% but with no vertical moving window. It can introduce varying  
% anisotropy to the analysis by changing vertical exaggeration of the 
% dataset. 

% Custom functions called by this script;
%   computeKfunctionForPattern.m
%   computeProportionOptimized.m
%   plotResults.m

clear ;

%% Load pointsets
load('CC_Centroids.mat'); %This is a variable saved with all the 
    % pointsets of the project. It is simply a 2-column dataset 
    % containing the x and y coordinates of the channel belt centroids. 
points = CC5_Centroids_Flattened; % assign the data to the variable 
    % points. Note that 'CC5_Centroids_Flattened' is one of two
    % datasets that I have saved within 'CC_Centroids.mat'

%% Input window lengths to run, and variation to vertical exaggeration.
movingWindowLengths = [700, 900, 1100, 1500, 1700, 1900, 2100, 2300,...
    2500, 2700, 2900]; % You can iterate through multiple window length
        % values here, or just pick one. must be evenly divisible by 
        % cell width (cw)
VEDivisor = [1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]; 
    % The VEDivisor determines the magnitude of vertical exaggeration
    % of the pointset. 0 is no VE. 1 is 100% VE based on mean belt 
    % width to thickness ratios (see Flood and Hampson 2015).
limitFactor = 0.5; %this value determines the allowance for the largest 
    % scale of observation. 0.5 means that the largest dimension of the
    % K-function search radius can be no larger than 0.5 of the 
    % shortest dimension of the dataset.
for ix = 1:length(movingWindowLengths) %These two for loops take it
    % through each moving window length, and each VE Divisor.
for ixx = 1:length(VEDivisor)
    
%% Vertically exaggerate the Y values (or don't)
veExag = 1; % this is an on/off switch for using Vertical Exaggeration. 
    % 1 is on, 0 is off.
if VEDivisor == 0 %if the VEDaivisor is 0, there isn't any V.E.
    veExag = 0;
end
if veExag == 1 % if using V.E., this loop runs. 
    fullVE = 32.2 ; % this number is the ratio of mean channel belt 
        %width to mean channel belt thickness (sensu flood and hampson)
    vertExag = fullVE * VEDivisor(ixx) ;
    veYs = points(:,2) .* vertExag ;
    veXs = points(:,1) ; % x values do not change
    vePoints = horzcat(veXs, veYs); % create new vertically exaggerated
        % set of points.
else
    vePoints = points ;
    vertExag = 0 ;
end
%% Code for lateral moving window
movingWindowIncrement = 50; % this is the amount that the moving window
    % jumps laterally after each run. 
movingWindowLength = movingWindowLengths(ix); %input Length of window 
    % must be divisible by cell width (cw)
xminimum = floor(min(vePoints(:, 1))); % find the range of ponits.
xmaximum = ceil(max(vePoints(:, 1)));
R = ceil(((xmaximum - xminimum) - movingWindowLength)/ ... 
    movingWindowIncrement); %determine total # of lateral windows
patternChecks = [];
saveVars = 1 ; % on/off switch to save variables/results to a .mat file
plotFig = 0 ; % on/off switch to plot a figure.

%%
for iii= 0:1:R  %Iterate through lateral windows
    
disp(['Lateral Window ',num2str(iii),' of ', num2str(R)])
%% Pull the window subset from points.
winStart = xminimum + (movingWindowIncrement*iii);  %Gets minimum X 
    % value of window
winEnd = (xminimum + movingWindowLength) + (movingWindowIncrement*iii);  
    % Gets maximum X value of window
pointPattern = [];  %this will be filled with only the points in the 
    % window
for i=1:size(vePoints,1) %determines if points are in window
    if vePoints(i,1)>winStart && vePoints(i,1)<winEnd 
       newRow = [vePoints(i,1) vePoints(i,2)];
       pointPattern = vertcat(pointPattern, newRow);    %adds points to
            % windowPoints.
    end
end

%% --------Find points and range----------

Ymin = min(pointPattern(:,2));
Ymax = max(pointPattern(:,2));
Xmin = min(pointPattern(:,1));
Xmax = max(pointPattern(:,1));

% These automatic figure titles below can help if you're making many 
% figures. But can be hinderances if not. Use with caution.
figTitle = ['CC1S-noVE-1-K-fn-', num2str(movingWindowLength), ...
    '-win-', num2str(winStart)]; % makes an automatic figure title.
totalFigTitle = ['CC1S-noVE-1-totalKfnFig-', ...
    num2str(movingWindowLength), '-win-', ...
    num2str(movingWindowIncrement),'inc']; 
%--------Required parameters-----------
N=size(pointPattern, 1); % Here we're setting up parameters to make a 
    % simulated point set with the same number of points and dimensions
    % as the real dataset.
width = Xmax - Xmin;
limity = (Ymax - Ymin); %scale limit, set to section thickness, this
    % will be modified by our limitFactor from earlier.
A=width*(limity); % area
if limity > width
    limit = width*limitFactor ;
else limit = limity*limitFactor ;
end
increment = 5; % increment by which we will increase K-function scale. 
allResults = [];
MonteCarloSimulations = 99; % number of simulations. 99 is common
    % because considered with the 1 real dataset, there are 100 total.
mySum = 0;
scaleIncrements = [];
for i=0:increment:limit
    scaleIncrements = vertcat(scaleIncrements, i); %appends 'i' onto 
        % last row of h
        % basically making a list of all scales to be investigated.
end

%% ------------generateSimulations--------------------
for i=1:MonteCarloSimulations
    randPoints = rand(N,2);     %makes an Nx2 martix of random numbers 
        % between 0 and 1
    randPoints(:,1) = Xmin + width.*randPoints(:,1); %multiplies width 
        % by each of the x Rands (from 0 to 1)
    randPoints(:,2) = Ymin + (Ymax - Ymin).*randPoints(:,2); %same as 
        % above but for y's.
    simulationResults = [];
    for h=0:increment:limit % this loop goes through the actual
            % K-function computation, on the random dataset.
        mySum = 0; Ko = 0; K = 0;
        for j=1:N
            for k=1:N
                if j~=k
                    distance = sqrt((randPoints(j,1) - ...
                        randPoints(k,1))^2 + ((randPoints(j,2) - ...
                        randPoints(k,2))^2)); %Eucledean distance 
                            % between sj and sk
                    realDistance = sqrt((randPoints(j,1) - ...
                        randPoints(k,1))^2 + (randPoints(j,2) - ...
                        randPoints(k,2))^2);
                    if distance <= h % this loop makes the edge
                            % correction using the script 
                            % 'computeProportionOptimized.m'
                        proportion = computeProportionOptimized...
                            (randPoints(j,1), randPoints(j,2), ...
                            distance, Xmin, Ymin, Xmax, Ymax);
                             % proportion = (2*pi - 0)/(2*pi); %use if 
                                % don't want to make correction.
                        mySum = mySum + (1/proportion);
                    end
                end
            end
        end
        %Compute K with corresponding h, append to simulationResults
        Ko = mySum * (A/(N*N));
        K = sqrt(Ko*sqrt(1)/pi) - h;
        simulationResults = [simulationResults; K];
    end
    allResults = horzcat(allResults, simulationResults);
end
simulations = allResults;

%------------generatePattern----------------
pattern = computeKfunctionForPattern(N, A, limit, increment, ...
    pointPattern, Xmin, Ymin, Xmax, Ymax); 
    % custom function that now computes the K-function for the real
    % data set.

%% -------Find where the pattern is greater than the simulations-------
% When the data K-function is greater than the maximum value of any
% simulation, or lower than the minimum value of any simulation, then 
% we know that result is statistically significant.
maxs = [];
mins = [];
for i=1:size(simulations, 1) %goes through each 'window' and takes the 
        % max and min out of any simulation.
    maxs = [maxs; max(simulations(i,:))];
    mins = [mins; min(simulations(i,:))];
end

for t=1:size(pattern, 1)
    patternCheck = [];
    if pattern(t,1) > maxs(t,1)
        patternCheck = [(iii), (pattern(t,1)-maxs(t,1)), ...
            scaleIncrements(t,1), pattern(t,1)];
        patternChecks = vertcat(patternChecks, patternCheck);
    elseif pattern(t,1) < mins(t,1)
        patternCheck = [(iii), (pattern(t,1)-mins(t,1)), ...
            scaleIncrements(t,1),pattern(t,1)];
        patternChecks = vertcat(patternChecks, patternCheck);
    else
        patternCheck = [(iii), 0, scaleIncrements(t,1), pattern(t,1)];
        patternChecks = vertcat(patternChecks, patternCheck);
    end
end

%% ------------plotResults--------------------
% calls on 'plotResults.m' which is a custom function to plot the  
% K-function.
if plotFig == 1 
hold off;
figure;
if veExag == 1
    veIncrement = increment/vertExag ;
    veLimit = limit/vertExag ;
else veIncrement = increment ;
    veLimit = limit;
end
plotResults(veLimit, veIncrement, simulations, pattern);
xlabel('Scale of Observation (m)'); %if using V.E. change to "Scale of
    % Observation Normalized to 1 Mean Channel Width"
ylabel('K-Function Value');
title(figTitle);
%     xlim([0 200]);
%     ylim([-20 30]);
%     fig = gcf;
%     fig.PaperUnits = 'inches';
%     fig.PaperPosition = [0 0 4 6];
%     print(figTitle,'-dpng','-r0');  %print saves this out to file
end %end plotFig
end %end lateral window
finalOutput = horzcat(patternChecks(:,1), (patternChecks(:,3).*2), ...
    patternChecks(:,2), patternChecks(:,4));
    % This makes a variable with all of the results.
if veExag == 1
    VEforName = fullVE*VEDivisor(ixx)*100 ;
else
    VEforName = 0 ;
end
saveName = ['cc5-kfn-',num2str(movingWindowLength),'wind-',...
    num2str(VEforName),'VE-',num2str(movingWindowIncrement),'latinc'];
    % this will be the name of the .mat file if you choose to save 
    % variables.
if saveVars == 1 % 1 to save the variables to file, 0 to not.
    save(saveName);
end
end %end VEDivisor
end %end windowLengths

