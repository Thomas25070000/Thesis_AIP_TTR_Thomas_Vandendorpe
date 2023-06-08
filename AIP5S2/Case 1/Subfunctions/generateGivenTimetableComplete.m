function [base_tt, hour_tt, complete_tt] = generateGivenTimetableComplete_v2(settings,blocksections,runningtimes,clearingtimes,data)

traintype = data.type;
direction = data.direction;
entrytime = data.entryHH * 60 + data.entryMM;   % in [min]
firsttimepoint = min(entrytime) * 60;           % in [s] 
maxHour = settings.disruption.duration + settings.disruption.nr_hours_after;

% For each train type, we have to generate the timetable.
timetable = [];
tt = 0; % counter
train_id_dir0 = 0;
train_id_dir1 = 0;

% Assuming full speed over all block sections
speedIC = settings.trains.speed.IC / 3.6; % convert to m/s
speedR = settings.trains.speed.R / 3.6; % convert to m/s

% Start with the trains in the direction towards right (direction 1)
offsetdirections = settings.TT.delta.opposite*60;
% -> this is the time that the IC in direction 1 leaves after the one in
% direction 0 (time = 0)!

% Start with IC trains in direction 1
rows = find(strcmp(traintype,'IC') & direction == 1);
runningIC = runningtimes.IC1.regular;
for rr = 1:length(rows)
    time = entrytime(rows(rr))*60;
    hour = floor((time - firsttimepoint)/3600);
    if hour < maxHour
        train_id_dir1 = train_id_dir1+1;
        for bb = 1:size(blocksections,1)
            tt = tt+1;
            train_id = hour * 1000 + train_id_dir1*10 + 1;
            timetable(tt).event_id = tt;
            timetable(tt).train_id = train_id;
            timetable(tt).train_type = 'IC';
            timetable(tt).direction = 1;
            timetable(tt).blocksection = bb;
            timetable(tt).arrival = time;
            time = time + runningIC(bb);
            timetable(tt).departure = time;
            timetable(tt).running = timetable(tt).departure - timetable(tt).arrival;
        end
    end
end

% R trains in direction 1
% Note that delta(dir1) may be negative!
rows = find(~strcmp(traintype,'IC') & direction == 1);
runningL = runningtimes.L1.regular;
for rr = 1:length(rows)
    time = entrytime(rows(rr))*60;
    hour = floor((time - firsttimepoint)/3600);
    if hour < maxHour
        train_id_dir1 = train_id_dir1+1;
        for bb = 1:size(blocksections,1)
            tt = tt+1;
            timetable(tt).event_id = tt;
            train_id = hour * 1000 + train_id_dir1*10 + 1;
            timetable(tt).train_id = train_id;
            timetable(tt).train_type = 'R';
            timetable(tt).direction = 1;
            timetable(tt).blocksection = bb;
            timetable(tt).arrival = time;
            time = time + runningL(bb);
            timetable(tt).departure = time;
            timetable(tt).running = timetable(tt).departure - timetable(tt).arrival;
        end
    end
end

%% Take the other direction!
% IC trains, assume these leave at the period exactly!
rows = find(strcmp(traintype,'IC') & direction == 0);
for rr = 1:length(rows)
    time =  entrytime(rows(rr))*60;
    hour = floor((time - firsttimepoint)/3600);
    runningIC = runningtimes.IC0.regular;
    if hour < maxHour
        train_id_dir0 = train_id_dir0+1;
        for bb = size(blocksections,1):-1:1
            tt = tt+1;
            timetable(tt).event_id = tt;
            train_id = hour * 1000 + train_id_dir0*10;
            timetable(tt).train_id = train_id;
            timetable(tt).train_type = 'IC';
            timetable(tt).direction = 0;
            timetable(tt).blocksection = bb;
            timetable(tt).arrival = time;
            time = time + runningIC(bb);
            timetable(tt).departure = time;
            timetable(tt).running = timetable(tt).departure - timetable(tt).arrival;
        end
    end
end

% R trains
rows = find(~strcmp(traintype,'IC') & direction == 0);
for rr = 1:length(rows)
    time = entrytime(rows(rr))*60;
    hour = floor((time - firsttimepoint)/3600);
    runningL = runningtimes.L0.regular;
    if hour < maxHour
        train_id_dir0 = train_id_dir0+1;
        for bb = size(blocksections,1):-1:1
            tt = tt+1;
            timetable(tt).event_id = tt;
            train_id = hour * 1000 + train_id_dir0*10;
            timetable(tt).train_id = train_id;
            timetable(tt).train_type = 'R';
            timetable(tt).direction = 0;
            timetable(tt).blocksection = bb;
            timetable(tt).arrival = time;
            time = time + runningL(bb);
            timetable(tt).departure = time;
            timetable(tt).running = timetable(tt).departure - timetable(tt).arrival;
        end
    end
end

% Make sure that the arrival times are all positive
add = -min([timetable.arrival]);
for ee = 1:size(timetable,2)
    timetable(ee).arrival = timetable(ee).arrival + add;
    timetable(ee).departure = timetable(ee).departure + add;
end


%% Calculate the blocking times
t_setup = settings.TT.blocktimes.setup;
t_release_R = settings.TT.blocktimes.afterR;
t_release_IC = settings.TT.blocktimes.afterIC;

% If there is no approach, we assume the same running time as in the first
% block section!
for ee = 1:size(timetable,2)
    % What is the approach time?
    try
        if timetable(ee).train_id ~= timetable(ee-1).train_id
            % This is the first action of the train.
            t_before = timetable(ee).departure - timetable(ee).arrival + t_setup;
        else
            % There has already been an event
            t_before = timetable(ee-1).departure - timetable(ee-1).arrival + t_setup;
        end
    catch
        % This is the first event of the timetable!
        t_before = timetable(ee).departure - timetable(ee).arrival + t_setup;
    end
    timetable(ee).start = timetable(ee).arrival - t_before;
    
    switch timetable(ee).train_type
        case 'IC'
            timetable(ee).finish = timetable(ee).departure + t_release_IC;
        case 'R'
            timetable(ee).finish = timetable(ee).departure + t_release_R;
    end
end


% This generated the basic timetable.
basic_timetable = struct2table(timetable);

base_tt = basic_timetable;
hour_tt = basic_timetable;
complete_tt = basic_timetable;











end