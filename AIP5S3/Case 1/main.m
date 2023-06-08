clear all;
close all;

pathSubfunctions = '/Users/thomasvandendorpe/Dropbox/Thesis/Code/AIP5_gitkraken_S2b/Case 1/Subfunctions';
addpath(genpath(pathSubfunctions));
addpath(genpath(pwd));
saveFigures = 1;
saveStats = 1;

%% Inputs and outputs
filename = 'Case 3 - dummy data';
filenameInput = 'Case 3 - dummy data';
sheet = 'Parameters';

filenameOutput = 'Case 3 - dummy data - statistics.xlsx';
sheetOutput = 'All';


[~,parameters,~] = xlsread(filename,sheet,'A:A');
% [param_values,~,~] = xlsread(filename,sheet,'B:B');
[multi_param_values,textv,~] = xlsread(filename,sheet,'B:Z');
textv = textv(1,:);

param_values = multi_param_values(:,1);
settings = createSetting_Case3_parallelmachines(parameters, param_values);
settings.general.caseName = '';
settings.saveStats = saveStats;
settings.general.subName = textv{1,1};    

% Create the blocksections based on the infra data
blocksections = generateBlockSections_case3(settings);
stops = retrieveStationInformation_case3(parameters, param_values, settings);

% Generate running and clearing times
runningtimes = generateRunningTimes_case3(blocksections, settings);

% Generate timetable

sheet = 'TimetableComplete';
[rawdata.num,rawdata.text,~] = xlsread(filename,sheet);
rawdata.num(:,2) = [];
rawdata.text(1,:) = [];
rawdata.text(:,[1,3:end]) = [];
allowedtrack = zeros(size(rawdata.num,1),size(settings.tracks,1));
for tr = 1:size(settings.tracks,1)
    try
        allowedtrack(:,tr) = rawdata.num(:,4+tr);
        allowedtrack(find(isnan(allowedtrack(:,tr))),tr) = 1;
    catch
        allowedtrack(:,tr) = ones(size(rawdata.num,1),1);
    end
end
rawdata.allowedtrack = allowedtrack;
[base_tt, hour_tt, complete_tt, traininfo] = generateGivenTimetableComplete_case3(settings,blocksections,runningtimes,rawdata);


% No need to slim the timetable if we don't use that option.



%% Settings for the plots
plotSettings.include_blocks = 1;
HH = 15;
MM = 30;
SS = 0;
plotSettings.firstTime = HH * 3600 + MM * 60 + SS;
plotSettings.tracks = [1 2];
plotSettings.grouptracks = 0;

% HOUR
timetable = hour_tt;
plotSettings.type = 'hour';
%[line, blocks] = plotTT_case3(timetable, blocksections, traininfo, settings, plotSettings);
%     if saveFigures
%         figname = [settings.general.caseName ' ' settings.general.subName ' original - 1h.fig'];
%         savefig(gcf,figname);
%     end


% COMPLETE
timetable = complete_tt;
plotSettings.type = 'complete';
plotSettings.grouptracks = 0;
%[line, blocks] = plotTT_case3(timetable, blocksections, traininfo, settings, plotSettings);
%     if saveFigures
%         figname = [settings.general.caseName ' ' settings.general.subName ' original - complete.fig'];
%         savefig(gcf,figname);
%     enden

%% DISRUPTION MODELLING
timetable = updateRunningTimesAndTT_case3(timetable, runningtimes, settings, traininfo, blocksections);


% Build and solve the model
[new_timetable, solu, traininfo, statistics, measures] = modelCase3_parallelmachines_full(timetable, blocksections, traininfo, settings);

% Plot the resulting and adjusted timetable
plotSettings.noDots = 0;
plotSettings.grouptracks = 1;
plotSettings.origTrack = 1;
if ~plotSettings.grouptracks
    alltracks = 1:size(settings.tracks,1);
    for tr = 1:size(settings.tracks,1)
        plotSettings.tracks = alltracks(tr);
        [line_new, blocks_new] = plotTT_new_case3(new_timetable, blocksections, traininfo, settings, plotSettings);
%         if saveFigures
%             figname = [settings.general.caseName ' ' settings.general.subName ' adjusted (track ' int2str(tr) ').fig'];
%             savefig(gcf,figname);
%         end
    end
else
    [line_new, blocks_new] = plotTT_new_case3(new_timetable, blocksections, traininfo, settings, plotSettings);
    plotSettings.grouptracks = 0;
    [line_new, blocks_new] = plotTT_new_case3(new_timetable, blocksections, traininfo, settings, plotSettings);
end
% 
% % Save the statistics and results!
% if saveStats
%     statsname = [settings.general.caseName ' ' settings.general.subName ' adjusted - stats.mat'];
%     try
%         save(statsname,'settings','statistics','traininfo','measures','new_timetable');
%     catch
%         disp('Statistics not generated, impossible to save them');
%     end
% end


% Write the statistics to Excel!
if saveStats
    writeStatsToExcel(filenameOutput, sheetOutput, settings, statistics);
end

    

