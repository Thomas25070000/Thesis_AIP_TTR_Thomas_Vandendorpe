function settings = createSetting_Case3_parallelmachines(parameters, param_values)

which = 'timetable sheet';
row = find(strcmp(parameters,which));
if ~isempty(row)
    settings.general.TTsheet = param_values(row);
else
    settings.general.TTsheet = 0;
end



%% Get infra data
% Number blocks
which = 'number block sections A';
row = find(strcmp(parameters,which));
settings.infrastructure.blocks.A = param_values(row);
which = 'number block sections D';
row = find(strcmp(parameters,which));
settings.infrastructure.blocks.D = param_values(row);
which = 'number block sections B';
row = find(strcmp(parameters,which));
settings.infrastructure.blocks.B = param_values(row);
which = 'number block sections C';
row = find(strcmp(parameters,which));
settings.infrastructure.blocks.C = param_values(row);
% Length parts
which = 'length part A';
row = find(strcmp(parameters,which));
settings.infrastructure.length.A = param_values(row);
which = 'length S1';
row = find(strcmp(parameters,which));
settings.infrastructure.length.S1 = param_values(row);
which = 'length part D';
row = find(strcmp(parameters,which));
settings.infrastructure.length.D = param_values(row);
which = 'length S2';
row = find(strcmp(parameters,which));
settings.infrastructure.length.S2 = param_values(row);
which = 'length part B';
row = find(strcmp(parameters,which));
settings.infrastructure.length.B = param_values(row);
which = 'length part C';
row = find(strcmp(parameters,which));
settings.infrastructure.length.C = param_values(row);

% Signalling speed if yellow
which = 'speed yellow';
row = find(strcmp(parameters,which));
settings.infrastructure.vyellow = param_values(row);
% Max. speed on the blocks
which = 'max. speed blocks';
row = find(strcmp(parameters,which));
settings.infrastructure.vmax = param_values(row);



% Tracks
which = 'Track number';
row = find(strcmp(parameters,which));
% Length of row indicates number of tracks!
Ntracks = length(row);
tracks = [];
for tt = 1:Ntracks
    tracks(tt).nr = param_values(row(tt));
    tracks(tt).vreg = param_values(row(tt)+1);
    tracks(tt).vdis = param_values(row(tt)+2);
    if tracks(tt).vdis == 0
        tracks(tt).closed = 1;
    else
        tracks(tt).closed = 0;
    end
end
settings.tracks = struct2table(tracks);
settings.closedtracks = find([tracks.closed]);

% Trains and timetable
which = 'speed R';
row = find(strcmp(parameters,which));
settings.trains.speed.R = param_values(row);
which = 'speed IC';
row = find(strcmp(parameters,which));
settings.trains.speed.IC = param_values(row);
which = 'acc. R';
row = find(strcmp(parameters,which));
settings.trains.acc.R = param_values(row);
which = 'acc. IC';
row = find(strcmp(parameters,which));
settings.trains.acc.IC = param_values(row);
which = 'decc. R';
row = find(strcmp(parameters,which));
settings.trains.decc.R = param_values(row);
which = 'decc. IC';
row = find(strcmp(parameters,which));
settings.trains.decc.IC = param_values(row);
which = 'length R';
row = find(strcmp(parameters,which));
if isempty(row)
    settings.trains.length.R = 0;
else
    settings.trains.length.R = param_values(row);
end
which = 'length IC';
row = find(strcmp(parameters,which));
if isempty(row)
    settings.trains.length.IC = 0;
else
    settings.trains.length.IC = param_values(row);
end


which = 'TT given hour';
row = find(strcmp(parameters,which));
settings.TT.givenHour = param_values(row);
which = 'TT given complete';
row = find(strcmp(parameters,which));
settings.TT.givenComplete = param_values(row);
which = 'TT restrict tracks';
row = find(strcmp(parameters,which));
settings.TT.restrictTracks = param_values(row);

% Minimum headways
which = 'min headway same dir';
row = find(strcmp(parameters,which));
settings.TT.headway.same = param_values(row);
which = 'min headway other dir';
row = find(strcmp(parameters,which));
settings.TT.headway.other = param_values(row);
which = 'min headway end of corridor left';
row = find(strcmp(parameters,which));
if isempty(row)
    settings.TT.headway.end_left = 0;
else
    settings.TT.headway.end_left = param_values(row) * 60;
end
which = 'min headway end of corridor right';
row = find(strcmp(parameters,which));
if isempty(row)
    settings.TT.headway.end_right = 0;
else
    settings.TT.headway.end_right = param_values(row) * 60;
end
which = 'no first approach';
row = find(strcmp(parameters,which));
if isempty(row)
    settings.TT.noFirstApproach = 0;
else
    settings.TT.noFirstApproach = param_values(row);
end
which = 'before track switch';
row = find(strcmp(parameters,which));
if isempty(row)
    settings.TT.headway.beforetrackswitch = 0;
else
    settings.TT.headway.beforetrackswitch = param_values(row) * 60;
end
which = 'after track switch';
row = find(strcmp(parameters,which));
if isempty(row)
    settings.TT.headway.aftertrackswitch = 0;
else
    settings.TT.headway.aftertrackswitch = param_values(row) * 60;
end


% For the blocking times
which = {'sight and reaction time','setup time'};
row = find(ismember(parameters,which));
settings.TT.blocktimes.setup = sum(param_values(row));
which = {'release time','clearing time R'};
row = find(ismember(parameters,which));
settings.TT.blocktimes.afterR = sum(param_values(row));
which = {'release time'};
row = find(ismember(parameters,which));
settings.TT.blocktimes.releaseT = sum(param_values(row));
which = {'release time','clearing time IC'};
row = find(ismember(parameters,which));
settings.TT.blocktimes.afterIC = sum(param_values(row));
% After stops
which = 'reserve time after stop';
row = find(strcmp(parameters,which));
if isempty(row)
    settings.TT.blocktimes.afterstop = 0;
else
    settings.TT.blocktimes.afterstop = param_values(row);
end
% Buffer
which = 'buffer on running time';
row = find(strcmp(parameters,which));
settings.TT.buffer = param_values(row)/100;
% Constrain solution space (e.g. not 5 trains of same direction after each other)
which = 'max after each other';
row = find(strcmp(parameters,which));
settings.TT.orderlimit = param_values(row);
try
    which = 'order limit also before';
    row = find(strcmp(parameters,which));
    settings.TT.orderlimit_before = param_values(row);
end
if isempty(row)
    settings.TT.orderlimit_before = 0;
end
try
    which = 'order limit also after';
    row = find(strcmp(parameters,which));
    settings.TT.orderlimit_after = param_values(row);
end
if isempty(row)
    settings.TT.orderlimit_after = 0;
end



%% Disruption
which = 'time duration';
row = find(strcmp(parameters,which));
settings.disruption.duration = param_values(row);
which = 'signals present';
row = find(strcmp(parameters,which));
settings.disruption.signals = param_values(row);
which = 'allowed off-balance (per hour)';
row = find(strcmp(parameters,which));
settings.disruption.offbalance = param_values(row);
which = 'off-balance aggregated';
row = find(strcmp(parameters,which));
settings.disruption.aggregateOB = param_values(row);
which = 'off-balance rolling horizon';
row = find(strcmp(parameters,which));
settings.disruption.OBrollinghorizon = param_values(row);
which = 'include after-hours';
row = find(strcmp(parameters,which));
settings.disruption.nr_hours_after = param_values(row);

try
    which = 'max delay';
    row = find(strcmp(parameters,which));
    settings.disruption.maxDelay = param_values(row);
end
which = 'min LOS direction 0';
row = find(strcmp(parameters,which));
if ~isempty(row)
    settings.disruption.minLOS_dir0 = param_values(row);
else
    settings.disruption.minLOS_dir0 = -1;
end
which = 'min LOS direction 1';
row = find(strcmp(parameters,which));
if ~isempty(row)
    settings.disruption.minLOS_dir1 = param_values(row);
else
    settings.disruption.minLOS_dir1 = -1;
end

% OF weights
try
    which = 'exit delay';
    row = find(strcmp(parameters,which));
    settings.weights.exitDelay = param_values(row);
end
try
    which = 'cancel';
    row = find(strcmp(parameters,which));
    settings.weights.cancel = param_values(row);
end
try
    which = 'entrance delay';
    row = find(strcmp(parameters,which));
    settings.weights.entranceDelay = param_values(row);
end
try
    which = 'pax weight dir 0';
    row = find(strcmp(parameters,which));
    settings.weights.pax_dir0 = param_values(row);
end
if isempty(row)
	settings.weights.pax_dir0 = 1;
end
try
    which = 'pax weight dir 1';
    row = find(strcmp(parameters,which));
    settings.weights.pax_dir1 = param_values(row);
end
if isempty(row)
	settings.weights.pax_dir1 = 1;
end
which = 'penalty track change';
row = find(strcmp(parameters,which));
if isempty(row)
	settings.weights.trackchange = 0;
else
    settings.weights.trackchange = param_values(row);
end

%% Constraints
% MinLOS
try
    which = 'minLOS';
    row = find(strcmp(parameters,which));
    settings.constraints.minLOS = param_values(row);
end
if isempty(row) && (settings.disruption.minLOS_dir1 > 0 || settings.disruption.minLOS_dir0 > 0)
	settings.constraints.minLOS = 1;
elseif isempty(row)
	settings.constraints.minLOS = 0;
end
which = 'allow stopskip';
row = find(strcmp(parameters,which));
if isempty(row)
	settings.constraints.allowstopskip = 0;
else
	settings.constraints.allowstopskip = param_values(row);
end
% OB
try
    which = 'include OB';
    row = find(strcmp(parameters,which));
    settings.constraints.OB = param_values(row);
end
if isempty(row) && (settings.disruption.offbalance < 5)
	settings.constraints.OB = 1;
elseif isempty(row)
	settings.constraints.OB = 0;
end
which = 'deptime own track';
row = find(strcmp(parameters,which));
if isempty(row)
	settings.constraints.deptimeowntrack = 1;
else
	settings.constraints.deptimeowntrack = param_values(row);
end
% maxConsecutive
try
    which = 'maxConsecutive';
    row = find(strcmp(parameters,which));
    settings.constraints.maxConsecutive = param_values(row);
end
if isempty(row) && (settings.general.modelVersion < 9)
	settings.constraints.maxConsecutive = 0;
end
which = 'include after hours';
row = find(strcmp(parameters,which));
if ~isempty(row) && (param_values(row) == 0)
	settings.disruption.nr_hours_after = 0;
end
% fixEntranceOrder
try
    which = 'fix order at entrances';
    row = find(strcmp(parameters,which));
    settings.constraints.fixEntranceOrder = param_values(row);
end
if isempty(row)
	settings.constraints.fixEntranceOrder = 1;
end
try
    which = 'fix order threshold';
    row = find(strcmp(parameters,which));
    settings.constraints.fixEntranceOrderThreshold = param_values(row);
end
if isempty(row)
	settings.constraints.fixEntranceOrderThreshold = 3600;
end
try
    which = 'fix order same journey time';
    row = find(strcmp(parameters,which));
    settings.constraints.fixEntranceOrderSameTime = param_values(row);
end
if isempty(row)
	settings.constraints.fixEntranceOrderSameTime = 1;
end
% Handle THA
try
    which = 'max delay THA';
    row = find(strcmp(parameters,which));
    settings.constraints.maxDelayTHA = param_values(row);
end
if isempty(row)
	settings.constraints.maxDelayTHA = 3600;
end
try
    which = 'cancel THA allowed';
    row = find(strcmp(parameters,which));
    settings.constraints.cancelTHAallowed = param_values(row);
end
if isempty(row)
	settings.constraints.cancelTHAallowed = 1;
end
which = 'allow border crossing';
row = find(strcmp(parameters,which));
if ~isempty(row)
	settings.constraints.allowbordercross = param_values(row);
else
    settings.constraints.allowbordercross = 1;
end
which = 'no gamma';
row = find(strcmp(parameters,which));
if ~isempty(row)
	settings.constraints.nogamma = param_values(row);
else
    settings.constraints.nogamma = 0;
end

which = 'use min. buffer';
row = find(strcmp(parameters,which));
if ~isempty(row)
	settings.constraints.minbuffer = param_values(row);
else
    settings.constraints.minbuffer = 0;
end


%% Decisions
% WIC
which = 'WIC allowed';
row = find(strcmp(parameters,which));
if ~isempty(row)
    settings.disruption.WICallowed = param_values(row);
else
    settings.disruption.WICallowed = 0;
end
which = 'gamma at start';
row = find(strcmp(parameters,which));
if ~isempty(row)
    settings.disruption.gammAtStart = param_values(row);
else
    settings.disruption.gammAtStart = 0;
end
which = 'noCancel dir 0';
row = find(strcmp(parameters,which));
if ~isempty(row)
    settings.disruption.noCancel_dir0 = param_values(row);
else
    settings.disruption.noCancel_dir0 = 0;
end
which = 'noCancel dir 1';
row = find(strcmp(parameters,which));
if ~isempty(row)
    settings.disruption.noCancel_dir1 = param_values(row);
else
    settings.disruption.noCancel_dir1 = 0;
end

%% Deviation
which = 'deviation route available';
row = find(strcmp(parameters,which));
if ~isempty(row)
    settings.deviation.available = param_values(row);
else
    settings.deviation.available = 0;
end
% Direction 0
which = 'nr trains dir 0';
row = find(strcmp(parameters,which));
if ~isempty(row)
    settings.deviation.cap_dir0 = param_values(row);
else
    settings.deviation.cap_dir0 = 0;
end
which = 'delay dir 0';
row = find(strcmp(parameters,which));
if ~isempty(row)
    settings.deviation.delay_dir0 = param_values(row);
else
    settings.deviation.delay_dir0 = 0;
end
% Direction 1
which = 'nr trains dir 1';
row = find(strcmp(parameters,which));
if ~isempty(row)
    settings.deviation.cap_dir1 = param_values(row);
else
    settings.deviation.cap_dir1 = 0;
end
which = 'delay dir 1';
row = find(strcmp(parameters,which));
if ~isempty(row)
    settings.deviation.delay_dir1 = param_values(row);
else
    settings.deviation.delay_dir1 = 0;
end

which = 'time limit CPU (s)';
row = find(strcmp(parameters,which));
settings.solver.timelimit = param_values(row);
if ~isempty(row)
    settings.solver.timelimit = param_values(row);
else
    settings.solver.timelimit = 0;
end
which = 'compare with time limit';
row = find(strcmp(parameters,which));
settings.solver.comparetimelimit = param_values(row);
if ~isempty(row)
    settings.solver.comparetimelimit = param_values(row);
else
    settings.solver.comparetimelimit = 0;
end
which = 'MIP focus';
row = find(strcmp(parameters,which));
if ~isempty(row)
    settings.solver.mipfocus = param_values(row);
else
    settings.solver.mipfocus = 0;
end

which = 'adjust bigM';
row = find(strcmp(parameters,which));
if ~isempty(row)
    settings.solver.adjustBigM = param_values(row);
else
    settings.solver.adjustBigM = 0;
end



end