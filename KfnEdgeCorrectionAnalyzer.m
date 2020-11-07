% KfnEdgeCorrectionAnalyzer.m
% Developed by Alex Koch, contact at alexrkoch@gmail.com
% to analyze percentage of centroid points that require edge-correction 
% in the K-function, by K-radius scale and by moving window scale.
% Basically, this will pretend to run the K-function, but won't perform
% the actual calculation. Instead it will just count the centroids that
% would need correction.

clear;
runAnalysis = 1; % on/off switch for running the edge correction 
    % proportion analysis.
runAverages = 1; % on/off switch for creating summary figure. 

if runAnalysis == 1
load('CC_Centroids.mat'); % input centroid points
points = CC1_Centroids_Flattened; % set points to variable
movingWindowLengths = [600, 700, 800, 900, 1000, 1100, 1200, 1300];
    % choose moving window lengths to investigate. must be divisible by 
    % cell width (cw)
VEDivisor = 0;%[1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1];
    % The VEDivisor determines the magnitude of vertical exaggeration
    % of the pointset. 0 is no VE. 1 is 100% VE based on mean belt 
    % width to thickness ratios (see Flood and Hampson 2015).
limitFactor = 0.5; %this value determines the allowance for the largest 
    % scale of observation. 0.5 means that the largest dimension of the
    % K-function search radius can be no larger than 0.5 of the 
    % shortest dimension of the dataset.
tolerance = 0.5; % tolerance is the percentage of k-radius required to 
    % be out of bounds before it registers as "edge effect". So at 0.5
    % the code will only return an "edge effect" if the given point is 
    % closer than 1/2 of the k-distance from the edge of the dataset.
for ix = 1:length(movingWindowLengths)
for ixx = 1:length(VEDivisor)
finalOutput = [];
finalOutputAves = [];
    
%% Vertically exaggerate the Y values (or don't)
veExag = 0; % this is an on/off switch for using Vertical Exaggeration. 
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
movingWindowIncrement = 50 ; %this is the amount that the moving window
    % jumps laterally after each run. 
movingWindowLength = movingWindowLengths(ix) ; %input Length of window 
    % must be divisible by cell width (cw)
xminimum = floor(min(points(:, 1))) ; 
xmaximum = ceil(max(points(:, 1))) ;
R = ceil(((xmaximum - xminimum) - movingWindowLength)/ ...
    movingWindowIncrement); %determine total # of lateral windows
saveName = ['CC1-KfnEdgeEffect50T-',num2str(movingWindowLength), ...
    'wind-',num2str(vertExag*100),'VE.mat'];
for iii = 1:R
    
disp(['Lateral Window ',num2str(iii),' of ', num2str(R)])
%% Pull the window subset from points.
winStart = xminimum + (movingWindowIncrement*iii);  %Gets minimum X 
    % value of window
winEnd = (xminimum + movingWindowLength) + ...
    (movingWindowIncrement*iii);  %Gets maximum X value of window
windowPoints = [];  %this will be filled with only the points in 
    % the window
for i=1:size(vePoints,1) %determines if points are in window
    if vePoints(i,1)>winStart && vePoints(i,1)<winEnd 
       newRow = [vePoints(i,1) vePoints(i,2)];
       windowPoints = vertcat(windowPoints, newRow); %adds points 
          % to windowPoints.
    end
end
Ymin = min(windowPoints(:,2));
Ymax = max(windowPoints(:,2));
Xmin = min(windowPoints(:,1));
Xmax = max(windowPoints(:,1));
N=size(windowPoints, 1); % Here we're setting up parameters to make a 
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
outProportions = []; % will store the proportion of points meeting the
    % edge correction conditions.
for k = 0:increment:limit
    outPoints = 0;
    % the following loop calculates the distance of the points from all
    % data boundaries.
    for a = 1: length(windowPoints)
        leftDist = sqrt(((windowPoints(a,1) - winStart)^2) + ...
            ((windowPoints(a,2) - (windowPoints(a,2))^2)));
        rightDist = sqrt(((windowPoints(a,1) - winEnd)^2) + ...
            ((windowPoints(a,2) - (windowPoints(a,2))^2)));
        topDist = sqrt(((windowPoints(a,1) - windowPoints(a,1))^2)...
            + ((windowPoints(a,2) - Ymax)^2));
        baseDist = sqrt(((windowPoints(a,1) - windowPoints(a,1))^2)...
            + ((windowPoints(a,2) - Ymin)^2));
        if leftDist < (k*tolerance) || rightDist < (k*tolerance) ||...
                topDist < (k*tolerance) || baseDist < (k*tolerance)
            outPoints = outPoints + 1; % if any are greater then the 
                % distance * tolerance, they get counted with the
                % 'outPoints' variable.
        end
    end % end searching all points
    outProportion = outPoints / length(windowPoints);
    outProportions = vertcat(outProportions, outProportion); 
    newRow = [iii, k, outProportion];
    finalOutput = vertcat(finalOutput, newRow); % this will hold the 
        % results separated out by scale, as opposed to the average
        % method below (finalOutputAves).
end
aveOutProp = mean(outProportions); 
newRow = [iii, aveOutProp];
finalOutputAves = vertcat(finalOutputAves, newRow); % this will hold 
    % the average proportion of points needing edge correction from 
    % all scales for each moving window (iii). 


end %end lateral window
save(saveName, 'finalOutput', 'finalOutputAves');
end %end VEdivisor
end %end window length
end %end if runAnalysis

if runAverages == 1
%% Compile KfnEdgeCorrection Results into Master sheet to be plotted.
% Use this part of the script to make a figure which gives trendlines
% for various levels of V.E., which plot K-function scale on the X-axis
% and proportion of edge corrections on the Y-axis (see figure 11 of
% this thesis).

clear;
output = [];
movingWindowLengths = 1300; %[600, 800, 1000, 1200, 1400, 1600, 1800, 
    % 2000, 2200, 2400, 2600, 2800, 3000];
VEDivisor = [1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0];
for ix = 1:length(movingWindowLengths)
    figTitle = (['CC1-KfnEdgeEffect50T-', ... 
        num2str(movingWindowLengths(ix)),'wind-figure']);
    figure; hold on;
for ixx = 1:length(VEDivisor)
vertExag = 32.2 * VEDivisor(ixx) ; % again this number is the ratio of 
        % mean channel belt width to mean channel belt thickness (sensu 
        % Flood and Hampson 2015)
load(['CC1-KfnEdgeEffect50T-',num2str(movingWindowLengths(ix)),... 
    'wind-',num2str(vertExag*100),'VE.mat']);

finish = max(finalOutput(:,2));
scaleAverages = [];
for a = 5:5:finish
   sum = 0;
   counter = 0;
   for b = 1:length(finalOutput)
       if finalOutput(b,2) == a
           sum = sum + finalOutput(b,3);
           counter = counter + 1;
       end
   end
   average = sum / counter;
   newRow = [a, average];
   scaleAverages = vertcat(scaleAverages, newRow);
end
plot(scaleAverages(:,1),scaleAverages(:,2), 'LineWidth', 2)
% xlim([5, 300]) 
ylim([0 , 1])

aveProportion = mean(finalOutputAves(:,2));
newRow = [movingWindowLengths(ix), vertExag, aveProportion];
output = vertcat(output, newRow);

end %end VEDivisor
hold off;
legend('32.2x VE', '28.98x VE','25.76x VE', '22.54x VE','19.32x VE',...
    '16.1x VE','12.88x VE','9.66x VE','6.44x VE','3.22x VE',...
    '0x VE', 'Location','SouthEast');
% xlabel('K-function Width Scale (m)');
% ylabel('Percent of Centroids Needing Edge Correction');
% title(figTitle)
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 3 4];
print(figTitle,'-dtiff','-r1200')

end %end MWLength
end %end if runAverages
