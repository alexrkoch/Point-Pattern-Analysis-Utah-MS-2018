% MovingWindowCorrelation.m
% Developed by Alex Koch, contact at alexrkoch@gmail.com
% Repository at https://github.com/alexrkoch


% Use this code to statistically compare the results of the Quadrat
% analysis to any geologic statistic calculated over the same set of
% moving windows (both vertically and laterally). It calculates a
% Pearson's R^2 correlation coefficient between two properties for each
% moving window.

clear;

%% Control Panel

% what to run
% abbreviations; q = quadrat; n = ntg; w = width; t = thickness; 
%   wt = width-to-thickness ratio; ai = amalgamation index.
qVn = 1 ;
qVw = 0 ;
qVt = 0 ;
qVwt = 0 ;
qVai = 0 ;
nVw = 0 ;
nVt = 0 ;
nVwt = 0 ;
nVai = 0;
wVt = 0 ;
wVwt = 0 ;
wVai = 0 ;
tVwt = 0 ;
tVai = 0 ;
wtVai = 0 ;
SDwidthToThickVqtoNtg = 0;

% Naming info
OC = 'CC1';
version = '1000qVn-100thick'; 
xlimit = 13; 
ymin = -1;

% set histogram axis limits
hstXmin = -0.5;
hstXmax = 1;

saveVars = 1;
saveName = ([OC,'-qVn-',version,'-vars']);

%what to plot and print
plotLine = 0 ;
printLine = 0 ;
plotScatter = 1 ;
printScatter = 1 ;
plotHist = 1 ;
printHist = 1 ;
plotWindowNtgVsQuad = 0 ;
printWindowNtgVsQuad = 0;
plotMeanSD = 0;
printMeanSD = 0;

%% Input quadrat results
quadResults = [ ];

%% Input ntg results
ntgResults = [ ];

%% Input width results
wResults = [ ];

%% Input thickness results
tResults = [ ];

%% Input width:thickness results
wtResults = [ ];

%% Input amalgamation index results
aiResults = [ ];

%% Put quadrat results into a structure
quadS = struct('LatWind', [], 'Strat', [], 'Quad', []);

for a = min(quadResults(:,1)):1:max(quadResults(:,1)) %for each lateral 
        % window
    aa = a+1;
    quad = [];
    strat = [];
    for b = 1:length(quadResults) %for each vertical window result
        if quadResults(b,1) == a %is it within the lateral window in 
                % question?
            quadS(aa).LatWind = a;
            quad = vertcat(quad, quadResults(b,3));
            strat = vertcat(strat, quadResults(b,2));
        end
    end
    quadS(aa).Quad = quad;
    quadS(aa).Strat = strat;
end

%% Put ntg results into a structure
ntgS = struct('LatWind', [], 'Strat', [], 'Ntg', []);

for a = min(ntgResults(:,1)):1:max(ntgResults(:,1)) %for each lateral 
        % window
    aa = a+1;
    ntg = [];
    strat = [];
    for b = 1:length(ntgResults) %for each vertical window result
        if ntgResults(b,1) == a %is it within the lateral window in 
                % question?
            ntgS(aa).LatWind = a;
            ntg = vertcat(ntg, ntgResults(b,3));
            strat = vertcat(strat, ntgResults(b,2));
        end
    end
    ntgS(aa).Ntg = ntg;
    ntgS(aa).Strat = strat;
end

%% Put width results into a structure
wS = struct('LatWind', [], 'Strat', [], 'W', []);

for a = min(wResults(:,1)):1:max(wResults(:,1))%for each lateral window
    aa = a+1;
    w = [];
    strat = [];
    for b = 1:length(wResults) %for each vertical window result
        if wResults(b,1) == a %is it within the lateral window in 
                % question?
            wS(aa).LatWind = a;
            w = vertcat(w, wResults(b,3));
            strat = vertcat(strat, wResults(b,2));
        end
    end
    wS(aa).W = w;
    wS(aa).Strat = strat;
end

%% Put thickness results into a structure
tS = struct('LatWind', [], 'Strat', [], 'T', []);

for a = min(tResults(:,1)):1:max(tResults(:,1))%for each lateral window
    aa = a+1;
    t = [];
    strat = [];
    for b = 1:length(tResults) %for each vertical window result
        if tResults(b,1) == a %is it within the lateral window in 
                % question?
            tS(aa).LatWind = a;
            t = vertcat(t, tResults(b,3));
            strat = vertcat(strat, tResults(b,2));
        end
    end
    tS(aa).T = t;
    tS(aa).Strat = strat;
end

%% Put w:t results into a structure
wtS = struct('LatWind', [], 'Strat', [], 'WT', []);

for a = min(wtResults(:,1)):1:max(wtResults(:,1))%for each lateral window
    aa = a+1;
    wt = [];
    strat = [];
    for b = 1:length(wtResults) %for each vertical window result
        if wtResults(b,1) == a %is it within the lateral window in 
                % question?
            wtS(aa).LatWind = a;
            wt = vertcat(wt, wtResults(b,3));
            strat = vertcat(strat, wtResults(b,2));
        end
    end
    wtS(aa).WT = wt;
    wtS(aa).Strat = strat;
end

%% Put amalg index results into a structure
aiS = struct('LatWind', [], 'Strat', [], 'AI', []);

for a = min(aiResults(:,1)):1:max(aiResults(:,1)) %for each lateral 
        % window
    aa = a+1;
    ai = [];
    strat = [];
    for b = 1:length(aiResults) %for each vertical window result
        if aiResults(b,1) == a %is it within the lateral window in 
                % question?
            aiS(aa).LatWind = a;
            ai = vertcat(ai, aiResults(b,3));
            strat = vertcat(strat, aiResults(b,2));
        end
    end
    aiS(aa).AI = ai;
    aiS(aa).Strat = strat;
end

%% qVn

% Index for the scatter plots
latWindIndex = [];
for e = 1:1:length(quadS)
    latWindIndex = vertcat(latWindIndex, e);
end

%Compute and plot correlation coefficients
if qVn == 1
        qVnrVals = [];
        for c = min(ntgResults(:,1)):1:max(ntgResults(:,1)) %for each 
                % lateral window
            cc = c+1; %index plus 1 because lateral windows begin at 0
            corrCoefs = corrcoef(ntgS(cc).Ntg, quadS(cc).Quad);
            qVnrVals = vertcat(qVnrVals, (corrCoefs(1,2)));
        end


% Plots

    if plotScatter == 1
        figure;
        scatter(latWindIndex, qVnrVals, 75, 'k.');
%         title([OC,' Cluster Index to NTG Correlations Vs Lateral ... 
                %Window'],'FontSize', 8);
%         xlabel('Moving Window Position');
%         ylabel('Correlation Coefficient');
        xlim([1 xlimit]);
        ylim([ymin 1]);
        if printScatter == 1
            fig = gcf;
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 2 2];
            print([OC,'-qVn-',version,'-sctr'],'-dtiff','-r1200')
        end
    end
    if plotLine == 1
        figure;
        hold on;
        plot(qVnrVals)
        title('CC5F-quadratVsNtg-1-ln');
        xlabel('Moving Window Position');
        ylabel('Correlation Coefficient');
        plot([0,50],[0,0])
        if printLine == 1
            print('CC5F-qVn-1-ln','-dpng','-r0')
        end
    end

    if plotHist == 1
        figure;
        histogram(qVnrVals, 'FaceColor', 'k');
%         title([OC,' Cluster Index to NTG']);
%         xlabel('r^2 value');
%         ylabel('Frequency');
        xlim([-1 1])
        if printHist == 1
            fig = gcf;
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 2 2];
            print([OC,'-qVn-',version,'-hst'],'-dtiff','-r1200');
        end
    end

    if plotWindowNtgVsQuad == 1 
        figure;
        hold on;
        scatter((quadS(5).Quad), (ntgS(5).Ntg), 50, 'k.') %Insert 
            % lateral window index in the argument
        set(groot, 'DefaultLineLineWidth', 2) ;
        lsline
        title('Cluster Index Vs NTG Window 4')
        xlabel('Quadrat Clustering Index');
        ylabel('NTG');
        if printWindowNtgVsQuad == 1
            fig = gcf;
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 3 3];
            print('Cluster Index Vs NTG w4','-dtiff','-r1200');
        end
    end

    meanqVnrVal = mean(qVnrVals);
    SDqVnrVal = std(qVnrVals);
end

%% qVw 
% Index for the scatter plots
latWindIndex = [];
for e = 1:1:length(quadS)
    latWindIndex = vertcat(latWindIndex, e);
end

%Compute and plot correlation coefficients
if qVw == 1 
    if qVw == 1
        qVwrVals = [];
        for c = min(ntgResults(:,1)):1:max(ntgResults(:,1)) %for each 
                % lateral window
            cc = c+1;
            corrCoefs = corrcoef(quadS(cc).Quad, wS(cc).W);
            qVwrVals = vertcat(qVwrVals, (corrCoefs(1,2)));
        end
    end

    % qVw Plots

    if plotScatter == 1
        figure;
        scatter(latWindIndex, qVwrVals, 75, 'k.');
%         title('CC5 Cluster Index to Width Correlations Vs Lateral... 
                % Window');
%         xlabel('Moving Window Position');
%         ylabel('Correlation Coefficient');
        xlim([0 xlimit]);
        ylim([ymin 1]);
        if printScatter == 1
            fig = gcf;
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 2 2];
            print([OC,'-qVw-',version,'-sctr'],'-dtiff','-r1200')
        end
    end
    if plotLine == 1
        figure;
        plot(qVwrVals)
        title('CC5F-quadratVsWidth-1-ln');
        xlabel('Moving Window Position');
        ylabel('Correlation Coefficient');
        if printLine == 1
            print('CC5F-qVw-1-ln','-dpng','-r0')
        end
    end

    if plotHist == 1
        figure;
        histogram(qVwrVals, 'FaceColor', 'k');
%         title('CC5 Cluster Index to Width');
%         xlabel('Correlation Coefficients');
%         ylabel('Frequency');
        if printHist == 1
            fig = gcf;
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 2 2];
            print([OC,'-qVw-',version,'-hst'],'-dtiff','-r1200');
        end
    end

    if plotWindowNtgVsQuad == 1 
        figure;
        scatter((quadS(1).Quad), (wS(1).W)) %Insert lateral window 
            % index in the argument
        title('CC5F-quadratVsWidth-1');
        xlabel('Quadrat Clustering Index');
        ylabel('Width');
    end

    meanqVwrVal = mean(qVwrVals);
    SDqVwrVal = std(qVwrVals);
end

%% qVt 
% Index for the scatter plots
latWindIndex = [];
for e = 1:1:length(quadS)
    latWindIndex = vertcat(latWindIndex, e);
end

%Compute and plot correlation coefficients
if qVt == 1 
        qVtrVals = [];
        for c = min(ntgResults(:,1)):1:max(ntgResults(:,1)) %for each 
                % lateral window
            cc = c+1;
            corrCoefs = corrcoef(tS(cc).T, quadS(cc).Quad);
            qVtrVals = vertcat(qVtrVals, (corrCoefs(1,2)));
        end

    % qVt Plots

    if plotScatter == 1
        figure;
        scatter(latWindIndex, qVtrVals, 75, 'k.');
%         title('CC5 Cluster Index to Thickness Correlations Vs Lateral 
                % Window');
%         xlabel('Moving Window Position');
%         ylabel('Correlation Coefficient');
        xlim([0 xlimit]);
        ylim([ymin 1]);
        if printScatter == 1
            fig = gcf;
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 2 2];
            print([OC,'-qVt-',version,'-sctr'],'-dtiff','-r1200')
        end
    end
    if plotLine == 1
        figure;
        plot(qVtrVals)
        title('CC5F-quadratVsThickness-1-ln');
        xlabel('Moving Window Position');
        ylabel('Correlation Coefficient');
        if printLine == 1
            print('CC5F-qVt-1-ln','-dpng','-r0')
        end
    end

    if plotHist == 1
        figure;
        histogram(qVtrVals, 'FaceColor', 'k');
%         title('CC5 Cluster Index to Thickness');
%         xlabel('Correlation Coefficients');
%         ylabel('Frequency');
        if printHist == 1
            fig = gcf;
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 2 2];
            print([OC,'-qVt',version,'-hst'],'-dtiff','-r1200');
        end
    end

    if plotWindowNtgVsQuad == 1 
        figure;
        scatter((quadS(1).Quad), (tS(1).T)) %Insert lateral window 
            % index in the argument
        title('CC5F-quadratVsThickness-1');
        xlabel('Quadrat Clustering Index');
        ylabel('Thickness');
    end

    meanqVtrVal = mean(qVtrVals);
    SDqVtrVal = std(qVtrVals);
end

%% qVai 
% Index for the scatter plots
latWindIndex = [];
for e = 1:1:length(quadS)
    latWindIndex = vertcat(latWindIndex, e);
end
%Compute and plot correlation coefficients
if qVai == 1 
    qVairVals = [];
    for c = min(aiResults(:,1)):1:max(aiResults(:,1)) %for each lateral 
            % window
        cc = c+1;
        corrCoefs = corrcoef(aiS(cc).AI, quadS(cc).Quad);
        qVairVals = vertcat(qVairVals, (corrCoefs(1,2)));
    end
    
    % qVai Plots

    if plotScatter == 1
        figure;
        scatter(latWindIndex, qVairVals, 75, 'k.');
%         title('CC5 Cluster Index to Amalgamation Correlations Vs... 
                % Lateral Window');
%         xlabel('Moving Window Position');
%         ylabel('Correlation Coefficient');
        xlim([0 xlimit]);
        ylim([ymin 1]);
        if printScatter == 1
            fig = gcf;
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 2 2];
            print([OC,'-qVai-',version,'-sctr'],'-dtiff','-r1200')
        end
    end
    if plotLine == 1
        figure;
        plot(qVairVals)
        title('CC5F-quadratVsAmalgIndex-1-ln');
        xlabel('Moving Window Position');
        ylabel('Correlation Coefficient');
        if printLine == 1
            print('CC5F-qVai-1-ln','-dpng','-r0')
        end
    end

    if plotHist == 1
        figure;
        histogram(qVairVals, 'FaceColor', 'k');
%         title('CC5 Cluster Index to Amalgamation');
%         xlabel('Correlation Coefficients');
%         ylabel('Frequency');
        if printHist == 1
            fig = gcf;
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 2 2];
            print([OC,'-qVai-',version,'-hst'],'-dtiff','-r1200');
        end
    end

    if plotWindowNtgVsQuad == 1 
        figure;
        scatter((quadS(1).Quad), (aiS(1).AI))
        title('CC5F-quadratVsAmalgIndex-1');
        xlabel('Quadrat Clustering Index');
        ylabel('AmalgamationIndex');
    end

    meanqVairVal = mean(qVairVals);
    SDqVairVal = std(qVairVals);
end

%% nVw 
% Index for the scatter plots
latWindIndex = [];
for e = 1:1:length(ntgS)
    latWindIndex = vertcat(latWindIndex, e);
end

%Compute and plot correlation coefficients
if nVw == 1 
        nVwrVals = [];
        for c = min(ntgResults(:,1)):1:max(ntgResults(:,1)) %for each 
                % lateral window
            cc = c+1;
            corrCoefs = corrcoef(ntgS(cc).Ntg, wS(cc).W);
            nVwrVals = vertcat(nVwrVals, (corrCoefs(1,2)));
        end

    % nVw Plots

    if plotScatter == 1
        figure;
        scatter(latWindIndex, nVwrVals, 75, 'k.');
%         title('CC5F-ntgVsWidth-1-sctr');
%         xlabel('Moving Window Position');
%         ylabel('Correlation Coefficient');
        xlim([0 xlimit]);
        ylim([ymin 1]);
        if printScatter == 1
            fig = gcf;
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 3 2];
            print([OC,'-nVw-',version,'-sctr'],'-dtiff','-r1200')
        end
    end
    if plotLine == 1
        figure;
        plot(nVwrVals)
        title('CC5F-ntgVsWidth-1-ln');
        xlabel('Moving Window Position');
        ylabel('Correlation Coefficient');
        if printLine == 1
            print('CC5F-nVw-1-ln','-dpng','-r0')
        end
    end

    if plotHist == 1
        figure;
        histogram(nVwrVals, 'FaceColor', 'k');
%         title('CC5F-ntgVsWidth-1');
%         xlabel('Correlation Coefficients');
%         ylabel('Frequency');
        xlim([hstXmin hstXmax]);
        if printHist == 1
            fig = gcf;
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 2 2];
            print([OC,'-nVw-',version,'-hst'],'-dtiff','-r1200')
        end
    end

    if plotWindowNtgVsQuad == 1 
        figure
        scatter((ntgS(1).Ntg), (wS(1).W))
        title('CC5F-ntgVsWidth-1');
        xlabel('NTG');
        ylabel('Width');
    end

    meannVwrVal = mean(nVwrVals);
    SDnVwrVal = std(nVwrVals);
end

%% nVt 
% Index for the scatter plots
latWindIndex = [];
for e = 1:1:length(ntgS)
    latWindIndex = vertcat(latWindIndex, e);
end
%Compute and plot correlation coefficients
if nVt == 1 ;
        nVtrVals = [];
        for c = min(ntgResults(:,1)):1:max(ntgResults(:,1)) %for each 
                % lateral window
            cc = c+1;
            corrCoefs = corrcoef(ntgS(cc).Ntg, tS(cc).T);
            nVtrVals = vertcat(nVtrVals, (corrCoefs(1,2)));
        end

    % nVt Plots

    if plotScatter == 1;
        figure;
        scatter(latWindIndex, nVtrVals, 75, 'k.');
%         title('CC5F-ntgVsThickness-1-sctr');
%         xlabel('Moving Window Position');
%         ylabel('Correlation Coefficient');
        xlim([0 xlimit]);
        ylim([ymin 1]);
        if printScatter == 1
            fig = gcf;
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 3 2];
            print([OC,'-nVt-',version,'-sctr'],'-dtiff','-r1200')
        end
    end
    if plotLine == 1;
        figure;
        plot(nVtrVals)
        title('CC5F-ntgVsThickness-1-ln');
        xlabel('Moving Window Position');
        ylabel('Correlation Coefficient');
        if printLine == 1
            print('CC5F-nVt-1-ln','-dpng','-r0')
        end
    end

    if plotHist == 1;
        figure;
        histogram(nVtrVals, 'FaceColor', 'k');
%         title('CC5F-ntgVsThickness-1');
%         xlabel('Correlation Coefficients');
%         ylabel('Frequency');
        xlim([hstXmin hstXmax]);
        if printHist == 1;
            fig = gcf;
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 2 2];
            print([OC,'-nVt-',version,'-hst'],'-dtiff','-r1200');
        end
    end

    if plotWindowNtgVsQuad == 1 
        figure;
        scatter((ntgS(1).Ntg), (tS(1).T))
        title('CC5F-ntgVsThickness-1');
        xlabel('NTG');
        ylabel('Thickness');
    end

    meannVtrVal = mean(nVtrVals);
    SDnVtrVal = std(nVtrVals);
end

%% nVai 
% Index for the scatter plots
latWindIndex = [];
for e = 1:1:length(ntgS)
    latWindIndex = vertcat(latWindIndex, e);
end
%Compute and plot correlation coefficients
if nVai == 1 
        nVairVals = [];
        for c = min(ntgResults(:,1)):1:max(ntgResults(:,1)) %for each 
                % lateral window
            cc = c+1;
            corrCoefs = corrcoef(ntgS(cc).Ntg, aiS(cc).AI);
            nVairVals = vertcat(nVairVals, (corrCoefs(1,2)));
        end
        
    % nVai Plots

    if plotScatter == 1
        figure;
        scatter(latWindIndex, nVairVals, 75, 'k.');
%         title('CC5F-ntgVsAmalgIndex-1-sctr');
%         xlabel('Moving Window Position');
%         ylabel('Correlation Coefficient');
        xlim([0 xlimit]);
        ylim([ymin 1]);
        if printScatter == 1
            fig = gcf;
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 3 2];
            print([OC,'-nVai-',version,'-sctr'],'-dtiff','-r1200')
        end
    end
    if plotLine == 1
        figure;
        plot(nVairVals)
        title('CC5F-ntgVsAmalgIndex-1-ln');
        xlabel('Moving Window Position');
        ylabel('Correlation Coefficient');
        if printLine == 1
            print('CC5F-nVai-1-ln','-dpng','-r0')
        end
    end

    if plotHist == 1
        figure;
        histogram(nVairVals, 'FaceColor', 'k');
%         title('CC5F-ntgVsAmalgIndex-1');
%         xlabel('Correlation Coefficients');
%         ylabel('Frequency');
        xlim([hstXmin hstXmax]);
        if printHist == 1
            fig = gcf;
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 2 2];
            print([OC,'-nVai-',version,'-hst'],'-dtiff','-r1200');
        end
    end

    if plotWindowNtgVsQuad == 1 
        figure;
        scatter((ntgS(1).Ntg), (aiS(1).AI))
        title('CC5F-ntgVsAmalgIndex-1');
        xlabel('NTG');
        ylabel('AmalgamationIndex');
    end

    meannVairVal = mean(nVairVals);
    SDnVairVal = std(nVairVals);
end

%% tVai 
latWindIndex = [];
for e = 1:1:length(tS)
    latWindIndex = vertcat(latWindIndex, e);
end
%Compute and plot correlation coefficients
if tVai == 1 
        tVairVals = [];
        for c = min(ntgResults(:,1)):1:max(ntgResults(:,1)) %for each 
                % lateral window
            cc = c+1;
            corrCoefs = corrcoef(tS(cc).T, aiS(cc).AI);
            tVairVals = vertcat(tVairVals, (corrCoefs(1,2)));
        end
        
    % tVai Plots

    if plotScatter == 1
        figure;
        scatter(latWindIndex, tVairVals, 75, 'k.');
%         title('CC5F-thicknessVsAmalgIndex-1-sctr');
%         xlabel('Moving Window Position');
%         ylabel('Correlation Coefficient');
        xlim([0 xlimit]);
        ylim([ymin 1]);
        if printScatter == 1
            fig = gcf;
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 3 2];
            print([OC,'-tVai-',version,'-sctr'],'-dtiff','-r1200');
        end
    end
    if plotLine == 1
        figure;
        plot(tVairVals)
        title('CC5F-thicknessVsAmalgIndex-1-ln');
        xlabel('Moving Window Position');
        ylabel('Correlation Coefficient');
        if printLine == 1
            print('CC5F-tVai-1-ln','-dpng','-r0')
        end
    end

    if plotHist == 1
        figure;
        histogram(tVairVals, 'FaceColor', 'k');
%         title('CC5F-thicknessVsAmalgIndex-1');
%         xlabel('Correlation Coefficients');
%         ylabel('Frequency');
        xlim([hstXmin hstXmax]);
        if printHist == 1
            fig = gcf;
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 2 2];
            print([OC,'-tVai-',version,'-hst'],'-dtiff','-r1200');
        end
    end

    if plotWindowNtgVsQuad == 1 
        figure;
        scatter((tS(1).T), (aiS(1).AI))
        title('CC5F-thicknessVsAmalgIndex-1');
        xlabel('Thickness');
        ylabel('AmalgamationIndex');
    end

    meantVairVal = mean(tVairVals);
    SDtVairVal = std(tVairVals);
end

%% wVt 
% Index for the scatter plots
latWindIndex = [];
for e = 1:1:length(wS)
    latWindIndex = vertcat(latWindIndex, e);
end
%Compute and plot correlation coefficients
if wVt == 1 
    wVtrVals = [];
    for c = min(ntgResults(:,1)):1:max(ntgResults(:,1)) %for each 
            % lateral window
        cc = c+1;
        corrCoefs = corrcoef(wS(cc).W, tS(cc).T);
        wVtrVals = vertcat(wVtrVals, (corrCoefs(1,2)));
    end

    % wVt Plots

    if plotScatter == 1
        figure;
        scatter(latWindIndex, wVtrVals, 75, 'k.');
%         title('CC5F-widthVsThickness-1-sctr');
%         xlabel('Moving Window Position');
%         ylabel('Correlation Coefficient');
        xlim([0 xlimit]);
        ylim([ymin 1]);
        if printScatter == 1
            fig = gcf;
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 3 2];
            print([OC,'-wVt-',version,'-sctr'],'-dtiff','-r1200')
        end
    end
    if plotLine == 1
        figure;
        plot(wVtrVals)
        title('CC5F-widthVsThickness-1-ln');
        xlabel('Moving Window Position');
        ylabel('Correlation Coefficient');
        if printLine == 1
            print('CC5F-wVt-1-ln','-dpng','-r0')
        end
    end

    if plotHist == 1
        figure;
        histogram(wVtrVals, 'FaceColor', 'k');
%         title('CC5F-widthVsThickness-1');
%         xlabel('Correlation Coefficients');
%         ylabel('Frequency');
        xlim([hstXmin hstXmax]);
        if printHist == 1
            fig = gcf;
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 2 2];
            print([OC,'-wVt-',version,'-hst'],'-dtiff','-r1200');
        end
    end

    if plotWindowNtgVsQuad == 1 
        scatter((wS(1).W), (tS(1).T))
        title('CC5F-widthVsThickness-1');
        xlabel('Width');
        ylabel('Thickness');
    end

    meanwVtrVal = mean(wVtrVals);
    SDwVtrVal = std(wVtrVals);
end

%% wVai 
latWindIndex = [];
for e = 1:1:length(wS)
    latWindIndex = vertcat(latWindIndex, e);
end
%Compute and plot correlation coefficients
if wVai == 1 
        wVairVals = [];
        for c = min(ntgResults(:,1)):1:max(ntgResults(:,1)) %for each 
                % lateral window
            cc = c+1;
            corrCoefs = corrcoef(wS(cc).W, aiS(cc).AI);
            wVairVals = vertcat(wVairVals, (corrCoefs(1,2)));
        end
        
    % wVai Plots

    if plotScatter == 1
        figure;
        scatter(latWindIndex, wVairVals, 75, 'k.');
%         title('CC5F-widthVsAmalgIndex-1-sctr');
%         xlabel('Moving Window Position');
%         ylabel('Correlation Coefficient');
        xlim([0 xlimit]);
        ylim([ymin 1]);
        if printScatter == 1
            fig = gcf;
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 3 2];
            print([OC,'-wVai-',version,'-sctr'],'-dtiff','-r1200');
        end
    end
    if plotLine == 1;
        figure;
        plot(wVairVals)
        title('CC5F-widthVsAmalgIndex-1-ln');
        xlabel('Moving Window Position');
        ylabel('Correlation Coefficient');
        if printLine == 1
            print([OC,'-wVai-',version,'-hst'],'-dtiff','-r1200')
        end
    end

    if plotHist == 1;
        figure;
        histogram(wVairVals, 'FaceColor', 'k');
%         title('CC5F-widthVsAmalgIndex-1');
%         xlabel('Correlation Coefficients');
%         ylabel('Frequency');
        xlim([hstXmin hstXmax]);
        if printHist == 1;
            fig = gcf;
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 2 2];
            print([OC,'-wVai-',version,'-hst'],'-dtiff','-r1200');
        end
    end

    if plotWindowNtgVsQuad == 1 
        figure;
        scatter((wS(1).W), (aiS(1).AI))
        title('CC5F-widthVsAmalgIndex-1');
        xlabel('Width');
        ylabel('Amalgamation Index');
    end

    meanwVairVal = mean(wVairVals);
    SDwVairVal = std(wVairVals);
end



%% tVwt 
latWindIndex = [];
for e = 1:1:length(tS)
    latWindIndex = vertcat(latWindIndex, e);
end
%Compute and plot correlation coefficients
if tVwt == 1 ;
        tVwtrVals = [];
        for c = min(ntgResults(:,1)):1:max(ntgResults(:,1)) %for each 
                % lateral window
            cc = c+1;
            corrCoefs = corrcoef(tS(cc).T, wtS(cc).WT);
            tVwtrVals = vertcat(tVwtrVals, (corrCoefs(1,2)));
        end

    % tVwt Plots

    if plotScatter == 1;
        figure;
        scatter(latWindIndex, tVwtrVals, 75, 'k.');
        title('CC5F-ThicknessVsWidthToThickness-1-sctr');
        xlabel('Moving Window Position');
        ylabel('Correlation Coefficient');
        if printScatter == 1
            print('CC5F-tVwt-1-sctr','-dpng','-r0')
        end
    end
    if plotLine == 1;
        figure;
        plot(tVwtrVals)
        title('CC5F-ThicknessVsWidthToThickness-1-ln');
        xlabel('Moving Window Position');
        ylabel('Correlation Coefficient');
        if printLine == 1
            print('CC5F-tVwt-1-ln','-dpng','-r0')
        end
    end

    if plotHist == 1;
        figure;
        histogram(tVwtrVals, 'FaceColor', 'k');
        title('CC5F-ThicknessVsWidthToThickness-1');
        xlabel('Correlation Coefficients');
        ylabel('Frequency');
        if printHist == 1;
            print('CC5F-tVwt-1-hst','-dpng','-r0');
        end
    end

    if plotWindowNtgVsQuad == 1 
        figure;
        scatter((tS(1).T), (wtS(1).WT))
        title('CC5F-ThicknessVsWidthToThickness-1');
        xlabel('Thickness');
        ylabel('WidthtoThickness');
    end

    meantVwtrVal = mean(tVwtrVals);
    SDtVwtrVal = std(tVwtrVals);
end

%% wVwt
latWindIndex = [];
for e = 1:1:length(wS)
    latWindIndex = vertcat(latWindIndex, e);
end
%Compute and plot correlation coefficients
if wVwt == 1
        wVwtrVals = [];
        for c = min(ntgResults(:,1)):1:max(ntgResults(:,1)) %for each 
                % lateral window
            cc = c+1;
            corrCoefs = corrcoef(wS(cc).W, wtS(cc).WT);
            wVwtrVals = vertcat(wVwtrVals, (corrCoefs(1,2)));
        end

    % wVwt Plots

    if plotScatter == 1
        figure;
        scatter(latWindIndex, wVwtrVals, 75, 'k.');
%         title('CC5F-widthVsWidthToThickness-1-sctr');
%         xlabel('Moving Window Position');
%         ylabel('Correlation Coefficient');
        if printScatter == 1
            fig = gcf;
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 2 2];
            print([OC,'-wVwt-',version,'-sctr'],'-dtiff','-r1200')
        end
    end
    if plotLine == 1
        figure;
        plot(wVwtrVals)
        title('CC5F-widthVsWidthToThickness-1-ln');
        xlabel('Moving Window Position');
        ylabel('Correlation Coefficient');
        if printLine == 1
            print('CC5F-wVwt-1-ln','-dpng','-r0')
        end
    end

    if plotHist == 1
        figure;
        histogram(wVwtrVals, 'FaceColor', 'k');
%         title('CC5F-widthVsWidthToThickness-1');
%         xlabel('Correlation Coefficients');
%         ylabel('Frequency');
        if printHist == 1
            fig = gcf;
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 2 2];
            print([OC,'-wVwt-',version,'-hst'],'-dtiff','-r1200');
        end
    end

    if plotWindowNtgVsQuad == 1 
        figure;
        scatter((wS(1).W), (wtS(1).WT))
%         title('CC5F-widthVsWidthToThickness-1');
%         xlabel('Width');
%         ylabel('WidthToThickness');
    end

    meanwVwtrVal = mean(wVwtrVals);
    SDwVwtrVal = std(wVwtrVals);
end

%% nVwt 
% Index for the scatter plots
latWindIndex = [];
for e = 1:1:length(ntgS)
    latWindIndex = vertcat(latWindIndex, e);
end
%Compute and plot correlation coefficients
if nVwt == 1 
        nVwtrVals = [];
        for c = min(ntgResults(:,1)):1:max(ntgResults(:,1)) %for each 
                % lateral window
            cc = c+1;
            corrCoefs = corrcoef(ntgS(cc).Ntg, wtS(cc).WT);
            nVwtrVals = vertcat(nVwtrVals, (corrCoefs(1,2)));
        end

    % nVwt Plots

    if plotScatter == 1
        figure;
        scatter(latWindIndex, nVwtrVals, 75, 'k.');
%         title('CC5F-ntgVsWidthToThickness-1-sctr');
%         xlabel('Moving Window Position');
%         ylabel('Correlation Coefficient');
        if printScatter == 1
            print([OC,'-nVwt-',version,'-sctr'],'-dtiff','-r1200')
        end
    end
    if plotLine == 1
        figure;
        plot(nVwtrVals)
        title('CC5F-ntgVsWidthToThickness-1-ln');
        xlabel('Moving Window Position');
        ylabel('Correlation Coefficient');
        if printLine == 1
            print('CC5F-nVwt-1-ln','-dpng','-r0')
        end
    end

    if plotHist == 1
        figure;
        histogram(nVwtrVals, 'FaceColor', 'k');
%         title('CC5F-ntgVsWidthToThickness-1');
%         xlabel('Correlation Coefficients');
%         ylabel('Frequency');
        if printHist == 1
            print([OC,'-nVwt-',version,'-hst'],'-dtiff','-r1200');
        end
    end

    if plotWindowNtgVsQuad == 1 
        figure;
        scatter((ntgS(1).Ntg), (wtS(1).WT))
        title('CC5F-ntgVsWidthToThickness-1');
        xlabel('NTG');
        ylabel('WidthToThickness');
    end

    meannVwtrVal = mean(nVwtrVals);
    SDnVwtrVal = std(nVwtrVals);
end


%% wtVai 
latWindIndex = [];
for e = 1:1:length(wtS)
    latWindIndex = vertcat(latWindIndex, e);
end
%Compute and plot correlation coefficients
if wtVai == 1 ;
        wtVairVals = [];
        for c = min(ntgResults(:,1)):1:max(ntgResults(:,1)) 
            cc = c+1;
            corrCoefs = corrcoef(wtS(cc).WT, aiS(cc).AI);
            wtVairVals = vertcat(wtVairVals, (corrCoefs(1,2)));
        end

    % wtVai Plots
    if plotScatter == 1;
        figure;
        scatter(latWindIndex, wtVairVals, 75, 'k.');
        title('CC5F-widthToThicknessVAmalgIndex-1-sctr');
        xlabel('Moving Window Position');
        ylabel('Correlation Coefficient');
        if printScatter == 1
            print('CC5F-wtVai-1-sctr','-dpng','-r0')
        end
    end
    if plotLine == 1;
        figure;
        plot(wtVairVals)
        title('CC5F-widthToThicknessVAmalgIndex-1-ln');
        xlabel('Moving Window Position');
        ylabel('Correlation Coefficient');
        if printLine == 1
            print('CC5F-wtVai-1-ln','-dpng','-r0')
        end
    end

    if plotHist == 1;
        figure;
        histogram(wtVairVals, 'FaceColor', 'k');
        title('CC5F-widthToThicknessVAmalgIndex-1');
        xlabel('Correlation Coefficients');
        ylabel('Frequency');
        if printHist == 1;
            print('CC5F-wtVai-1-hst','-dpng','-r0');
        end
    end

    if plotWindowNtgVsQuad == 1 
        scatter((quadS(1).Quad), (ntgS(1).Ntg))
        title('CC5F-widthToThicknessVAmalgIndex-1');
        xlabel('Witdh to Thickness Ratio');
        ylabel('Amalgamation Index');
    end

    meanwtVairVal = mean(wtVairVals);
    SDwtVairVal = std(wtVairVals);
end

%% qVwt 
% Index for the scatter plots
latWindIndex = [];
for e = 1:1:length(quadS)
    latWindIndex = vertcat(latWindIndex, e);
end
%Compute and plot correlation coefficients
if qVwt == 1 
        qVwtrVals = [];
        for c = min(ntgResults(:,1)):1:max(ntgResults(:,1)) 
            cc = c+1;
            corrCoefs = corrcoef(wtS(cc).WT, quadS(cc).Quad);
            qVwtrVals = vertcat(qVwtrVals, (corrCoefs(1,2)));
        end

    % qVwt Plots
    if plotScatter == 1
        figure;
        scatter(latWindIndex, qVwtrVals, 75, 'k.');
%         title('CC5');
%         xlabel('Moving Window Position');
%         ylabel('Correlation Coefficient');
        if printScatter == 1
            fig = gcf;
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 2 2];
            print([OC,'-qVwt-',version,'-sctr'],'-dtiff','-r1200')
        end
    end
    if plotLine == 1
        figure;
        plot(qVwtrVals)
        title('CC5F-quadratVsWidthToThick-1-ln');
        xlabel('Moving Window Position');
        ylabel('Correlation Coefficient');
        if printLine == 1
            print('CC5F-qVwt-1-ln','-dpng','-r0')
        end
    end

    if plotHist == 1
        figure;
        histogram(qVwtrVals, 'FaceColor', 'k');
%         title('CC5 Cluster Index to W:T Ratio');
%         xlabel('Correlation Coefficients');
%         ylabel('Frequency');
        if printHist == 1
            fig = gcf;
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 2 2];
            print([OC,'-qVwt-',version,'-hst'],'-dtiff','-r1200');
        end
    end

    if plotWindowNtgVsQuad == 1 
        figure;
        scatter((quadS(1).Quad), (wtS(1).WT))
        title('CC5F-quadratVsNtg-1');
        xlabel('Quadrat Clustering Index');
        ylabel('WidthToThickness');
    end

    meanqVwtrVal = mean(qVwtrVals);
    SDqVwtrVal = std(qVwtrVals);
end


%% what percentage of windows are negatively correlated?
% belowZero = 0;
% for d = 1:length(qVnrVals)
%     if qVnrVals(d) <= 0
%         belowZero = belowZero + 1;
%     end
% end
% percentNegs = belowZero/length(qVnrVals);
% %%

%% To make a correlation with standard deviation error bars plot 
    % (q,n,w,t)
if plotMeanSD == 1 
    figure;
    x = [1,2,3,4,5,6];
    y = [meanqVnrVal, meanqVwrVal, meanqVtrVal, meannVwrVal, ...
        meannVtrVal, meanwVtrVal];
    sd = [SDqVnrVal/2, SDqVwrVal/2, SDqVtrVal/2, SDnVwrVal/2, ...
        SDnVtrVal/2, SDwVtrVal/2];
    f = errorbar(x,y,sd,'bp');
    xlim([0 7]);
    ylim([-0.5 1]);
    xticks([0 1 2 3 4 5 6])
    f.Marker = '.';
    f.MarkerSize = 15;
    f.Color = 'k';
end

if printMeanSD == 1
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 3 2];
    print([OC,'-MeanSD-',version],'-dtiff','-r1200');
end

if qVn == 1
    minqVnrVals = min(qVnrVals);
end
if qVw == 1
    minqVwrVals = min(qVwrVals);
end
if qVt == 1
    minqVtrVals = min(qVtrVals);
end
if nVw == 1
    minnVwrVals = min(nVwrVals);
end
if nVt == 1
    minnVtrVals = min(nVtrVals);
end
if wVt == 1
    minwVtrVals = min(wVtrVals);
end



if saveVars == 1
    save(saveName)
end











