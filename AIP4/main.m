clear all;
close all;

addpath(genpath(pwd));
pathSubfunctions = '/Users/thomasvandendorpe/Dropbox/Thesis/Code/AIP4_S1/Subfunctions';
addpath(genpath(pathSubfunctions));
addpath(genpath(pwd));
saveFigures = 1;
saveStats = 1;

%% Inputs and outputs


maxCC = 4;



%     filename = 'Case 1 - input Tienen - morning';
filename = 'Case 1 - dummy data';
sheet = 'Parameters';

filenameOutput = 'Case 1 - dummy data - statistics.xlsx';
sheetOutput = 'All';

[~,parameters,~] = xlsread(filename,sheet,'A:A');
[multi_param_values,textv,~] = xlsread(filename,sheet,'B:Z');
textv = textv(1,:);

% Read the values
param_values = multi_param_values(:,1);
settings = createSetting_Case1_singlemachine(parameters, param_values);


% Give additional settings
settings.general.caseName = 'Tienen';
settings.saveStats = saveStats;
settings.general.subName = textv{1,1};


subName = settings.general.subName;

%% Now an initial timetable can be generated.
% Create the blocksections based on the infra data
blocksections = generateBlockSections(settings);
% Close blocksections
blocksections = closeBlockSections(blocksections,settings);

if (settings.trains.length.IC == 0) || (settings.trains.length.R == 0)
	% Generate running times for both train types and all possible
	% combinations! E.g. regular, disrupted, acc, decc, IC, L ...
	runningtimes = generateRunningTimes(blocksections, settings);
	Nblocks = size(blocksections,2);
	clearingtimes.L1.regular = repmat(settings.TT.blocktimes.afterR,1,Nblocks);
	clearingtimes.L0.regular = repmat(settings.TT.blocktimes.afterR,1,Nblocks);
	clearingtimes.L1.disrupted = repmat(settings.TT.blocktimes.afterR,1,Nblocks);
	clearingtimes.L0.disrupted = repmat(settings.TT.blocktimes.afterR,1,Nblocks);
	clearingtimes.IC1.regular = repmat(settings.TT.blocktimes.afterIC,1,Nblocks);
	clearingtimes.IC0.regular = repmat(settings.TT.blocktimes.afterIC,1,Nblocks);
	clearingtimes.IC1.disrupted = repmat(settings.TT.blocktimes.afterIC,1,Nblocks);
	clearingtimes.IC0.disrupted = repmat(settings.TT.blocktimes.afterIC,1,Nblocks);
else
	% Also generate the clearing times!
	[runningtimes, clearingtimes] = generateRunningAndClearingTimes(blocksections, settings);
end
% Create the timetable
if settings.TT.givenComplete
	sheet = 'TimetableComplete';
	[rawdata.num,rawdata.text,~] = xlsread(filename,sheet);
	data.direction = rawdata.num(:,1);
	data.entryHH = rawdata.num(:,3);
	data.entryMM = rawdata.num(:,4);
%         data.stops = rawdata.num(:,7);
%         data.stops(find(isnan(data.stops))) = 0;
	rawdata.text(1,:) = [];
	data.type = rawdata.text(:,2);
	[base_tt, hour_tt, complete_tt] = generateGivenTimetableComplete_v2(settings,blocksections,runningtimes,clearingtimes,data);
elseif settings.TT.givenHour
	sheet = 'TimetableHour';
	[rawdata.num,rawdata.text,~] = xlsread(filename,sheet);
	rawdata.num(:,2) = [];
	rawdata.text(1,:) = [];
	rawdata.text(:,[1,3]) = [];
	[base_tt, hour_tt, complete_tt] = generateGivenTimetableHour_v2(settings,blocksections,runningtimes,clearingtimes,rawdata);
else        
	[base_tt, hour_tt, complete_tt] = generateTimetable_v2(settings,blocksections,runningtimes,clearingtimes);
end

timetable = complete_tt;
allblocks = blocksections;

% If there are still blocks A and B, remove these from the timetable!
[timetable, blocksections, runningtimes] = slimTimetable(timetable,blocksections,runningtimes);


HH = 15;
MM = 54;
SS = 0;
firstTime = HH * 3600 + MM * 60 + SS;
[arrtime, arrtimeHHMMSS, deptime, deptimeHHMMSS] = retrieveArrivalAndDepartureTimes(timetable, firstTime);
settings.firstTime = firstTime;

% Extract a specific hour for the timetable
hour = 1;
hour_tt = getHourFromTimetable(complete_tt,settings,hour);
% This one gives the most busy hour!
hour_tt = getHourFromTimetable(complete_tt,settings);

% Plot and save figures
include_blocks = 1;
[line, blocks] = plotTT(hour_tt, allblocks, settings, 'hour', include_blocks, firstTime);
if saveFigures
	figname = [settings.general.caseName ' ' settings.general.subName ' original - 1h.fig'];
	savefig(gcf,figname);
end
[line, blocks] = plotTT(complete_tt, allblocks, settings, 'complete', include_blocks, firstTime);
if saveFigures
	figname = [settings.general.caseName ' ' settings.general.subName ' original - complete.fig'];
	savefig(gcf,figname);
end

%[line, blocks] = plotTT_to_explain_setup_times(complete_tt, allblocks, settings, 'hour', 1, firstTime)


regular.TT = timetable;

%% Apply blockage
% Due to the blockage, the running times will increase, which make the
% timetable (probably) infeasible.


% Update the running times
%     timetable = updateRunningTimes(timetable, blocksections, settings);
timetable = updateRunningTimes_v2(timetable, runningtimes, settings);

% Create headway matrix for closed section
[minHW, trains] = createHeadwayMatrixClosedSection(timetable, blocksections, settings);
[minHW_reg] = createHeadwayMatrixRegularTT(regular.TT, blocksections, settings);

settings.clearingtimes = clearingtimes;


regular.HW = minHW_reg;
%[new_timetable, measures, statistics] = scheduleFIFO(timetable, blocksections, settings)
[timetable, solu, statistics] = modelCase1_singleMachine_v2(timetable, blocksections, trains, minHW, settings);
%         figname = 'Tienen - evening - FIFO - low vD';ch
%         title('Tienen - evening - FIFO - low vD');

%         
% Save the adjusted figure
if saveFigures
	figname = [settings.general.caseName ' ' settings.general.subName ' adjusted.fig'];
	savefig(gcf,figname);
end

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
