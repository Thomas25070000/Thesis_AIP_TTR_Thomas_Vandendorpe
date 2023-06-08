.clear all;
close all;

pathSubfunctions = '/Users/thomasvandendorpe/Dropbox/Thesis/Code/AIP5_gitkraken_S2a/Case 1/Subfunctions';
addpath(genpath(pathSubfunctions));
addpath(genpath(pwd));
saveFigures = 1;
saveStats = 1;

%% Inputs and outputs
filename = 'Case 2 - dummy data';
filenameInput = 'Case 2 - dummy data';
sheet = 'Parameters';

filenameOutput = 'Case 2 - dummy data - statistics.xlsx';
sheetOutput = 'All';

[~,parameters,~] = xlsread(filenameInput,sheet,'A:A');
[multi_param_values,textv,~] = xlsread(filenameInput,sheet,'B:Z');
textv = textv(1,:);
    
close all
% Read the values
param_values = multi_param_values(:,1);
settings = createSetting_Case1_singlemachine(parameters, param_values);
% Give additional settings
settings.general.caseName = 'Dummy';
settings.saveStats = saveStats;
settings.general.subName = textv{1,1};

%% Now an initial timetable can be generated.
% Create the blocksections based on the infra data
blocksections = generateBlockSections(settings);
% Close blocksections
blocksections = closeBlockSections(blocksections,settings);
 
% Also generate the clearing times!
[runningtimes, clearingtimes] = generateRunningAndClearingTimes(blocksections, settings);

% Create the timetable
sheet = 'TimetableComplete';
[rawdata.num,rawdata.text,~] = xlsread(filename,sheet);
data.direction = rawdata.num(:,1);
data.entryHH = rawdata.num(:,3);
data.entryMM = rawdata.num(:,4);
rawdata.text(1,:) = [];
data.type = rawdata.text(:,2);
[base_tt, hour_tt, complete_tt] = generateGivenTimetableComplete_v2(settings,blocksections,runningtimes,clearingtimes,data);

timetable = complete_tt;
allblocks = blocksections;
complete_runningtimes = runningtimes;

% If there are still blocks A/B and C, remove these from the timetable!
[timetable, blocksections, runningtimes] = slimTimetable(timetable,blocksections,runningtimes);

HH = 15;
MM = 30;
SS = 0;
firstTime = HH * 3600 + MM * 60 + SS;
[arrtime, arrtimeHHMMSS, deptime, deptimeHHMMSS] = retrieveArrivalAndDepartureTimes(timetable, firstTime);
settings.firstTime = firstTime;

% % Extract a specific hour for the timetable
% hour = 1;
% % This one gives the most busy hour!
% hour_tt = getHourFromTimetable(complete_tt,settings);
% 
% % Plot and save figures
% include_blocks = 1;
% [line, blocks, YTick, YTickHHMMSS] = plotTT_0(complete_tt, allblocks, settings, 'complete', include_blocks, firstTime);
% [line, blocks] = plotTT_1(complete_tt, allblocks, settings, 'complete', include_blocks, firstTime, YTick, YTickHHMMSS);
% [line, blocks] = plotTT(hour_tt, allblocks, settings, 'hour', include_blocks, firstTime);

% if saveFigures
%     figname = [settings.general.caseName ' ' settings.general.subName ' original - 1h.fig'];
%     savefig(gcf,figname);
% end
% 
% regular.TT = timetable;

%[line, blocks] = plotTT_to_explain_setup_times(complete_tt, allblocks, settings, 'hour', 1, firstTime)


%% Apply blockage
% Due to the blockage, the running times will increase, which make the
% timetable (probably) infeasible.

%Update the running times
timetable = updateRunningTimes_v2(timetable, runningtimes, settings);
full_timetable = updateRunningTimes_v2(complete_tt,complete_runningtimes,settings);
regular.TT = timetable;
% 
% Schedule FIFO
[minHW_reg] = createHeadwayMatrixRegularTT(full_timetable, allblocks, settings);
settings.clearingtimes = clearingtimes;
regular.HW = minHW_reg;

[new_timetable, measures, statistics] = scheduleFIFO(full_timetable, allblocks, settings);

% Save the statistics and measures.
if saveStats
	statsname = [settings.general.caseName ' ' settings.general.subName ' adjusted - stats.mat'];
	try
		save(statsname,'settings','statistics','measures','new_timetable');
	catch
		disp('Statistics not generated, impossible to save them');
	end
end

% Write the statistics to Excel!
if saveStats
	writeStatsToExcel(filenameOutput, sheetOutput, settings, statistics);
end

% 
% % Schedule Machine
% [minHW, trains] = createHeadwayMatrixClosedSection(full_timetable, allblocks, settings);
% settings.clearingtimes = clearingtimes;        
% 
% [timetable, solu, statistics] = modelCase1_singleMachine_v2(full_timetable, allblocks, trains, minHW, settings);
% 
% 
% % Save the statistics and measures.
% if saveStats
% 	statsname = [settings.general.caseName ' ' settings.general.subName ' adjusted - stats.mat'];
% 	try
% 		save(statsname,'settings','statistics','measures','new_timetable');
% 	catch
% 		disp('Statistics not generated, impossible to save them');
% 	end
% end
% 
% % Write the statistics to Excel!
% if saveStats
% 	writeStatsToExcel(filenameOutput, sheetOutput, settings, statistics);
% end


    

