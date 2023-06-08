function [base_tt, hour_tt, complete_tt] = generateGivenTimetableHour(settings,blocksections,runningtimes,rawdata)

traintype = rawdata.text;
direction = rawdata.num(:,1);
entrytime = rawdata.num(:,2);

% For each train type, we have to generate the timetable.
timetable = [];
tt = 0; % counter
train_id_dir0 = 0;
train_id_dir1 = 0;

% Assuming full speed over all block sections
speedIC = settings.trains.speed.IC / 3.6; % convert to m/s
speedR = settings.trains.speed.R / 3.6; % convert to m/s

% Start with the trains in the direction towards right (direction 1)
% offsetdirections = settings.TT.delta.opposite*60;
% -> this is the time that the IC in direction 1 leaves after the one in
% direction 0 (time = 0)!

% Start with IC trains in direction 1
rows = find(strcmp(traintype,'IC') & direction == 1);
runningIC = runningtimes.IC1.regular;
for rr = 1:length(rows)
    time = entrytime(rows(rr))*60;
    train_id_dir1 = train_id_dir1+1;
    for bb = 1:size(blocksections,1)
        tt = tt+1;
        train_id = train_id_dir1*10+1;
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

% R trains in direction 1
% Note that delta(dir1) may be negative!
rows = find(~strcmp(traintype,'IC') & direction == 1);
runningL = runningtimes.L1.regular;
for rr = 1:length(rows)
    time = entrytime(rows(rr))*60;
    train_id_dir1 = train_id_dir1+1;
    for bb = 1:size(blocksections,1)
        tt = tt+1;
        timetable(tt).event_id = tt;
        train_id = train_id_dir1*10+1;
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

%% Take the other direction!
% IC trains, assume these leave at the period exactly!
rows = find(strcmp(traintype,'IC') & direction == 0);
for rr = 1:length(rows)
    time =  entrytime(rows(rr))*60;
    runningIC = runningtimes.IC0.regular;
    train_id_dir0 = train_id_dir0+1;
    for bb = size(blocksections,1):-1:1
        tt = tt+1;
        timetable(tt).event_id = tt;
        train_id = train_id_dir0*10;
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

% R trains
rows = find(~strcmp(traintype,'IC') & direction == 0);
for rr = 1:length(rows)
    time = entrytime(rows(rr))*60;
    runningL = runningtimes.L0.regular;
    train_id_dir0 = train_id_dir0+1;
    for bb = size(blocksections,1):-1:1
        tt = tt+1;
        timetable(tt).event_id = tt;
        train_id = train_id_dir0*10;
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

% Make sure that the arrival times are all positive
add = max(-min([timetable.arrival]),0);
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

%% Now duplicate this one to get the complete timetable.
base_tt = basic_timetable;
% Create a basic hour
BHP = base_tt;

% % Depends on the frequency per type
% freq_R_1 = settings.TT.frequency.R1;
% increment = 3600/freq_R_1;
% ff = 1;
% while ff < freq_R_1
%     
%     % Which rows?
%     rows = find(strcmp(base_tt.train_type,'R') & base_tt.direction == 1);
%     
%     new_part = base_tt(rows,:);
%     new_part.arrival = new_part.arrival + ff*increment;
%     new_part.departure = new_part.departure + ff*increment;
%     new_part.start = new_part.start + ff*increment;
%     new_part.finish = new_part.finish + ff*increment;
%     new_part.train_id = new_part.train_id+ff*100;
%     
%     
%     BHP = [BHP; new_part];
%     
%     ff = ff+1;
% end
% 
% % Depends on the frequency per type
% freq_IC_1 = settings.TT.frequency.IC1;
% increment = 3600/freq_IC_1;
% ff = 1;
% while ff < freq_IC_1
%     
%     % Which rows?
%     rows = find(strcmp(base_tt.train_type,'IC') & base_tt.direction == 1);
%     
%     new_part = base_tt(rows,:);
%     new_part.arrival = new_part.arrival + ff*increment;
%     new_part.departure = new_part.departure + ff*increment;
%     new_part.start = new_part.start + ff*increment;
%     new_part.finish = new_part.finish + ff*increment;
%     new_part.train_id = new_part.train_id+ff*100;
%     
%     
%     BHP = [BHP; new_part];
%     
%     ff = ff+1;
% end
% 
% % Depends on the frequency per type
% freq_R_0 = settings.TT.frequency.R0;
% increment = 3600/freq_R_0;
% ff = 1;
% while ff < freq_R_0
%     
%     % Which rows?
%     rows = find(strcmp(base_tt.train_type,'R') & base_tt.direction == 0);
%     
%     new_part = base_tt(rows,:);
%     new_part.arrival = new_part.arrival + ff*increment;
%     new_part.departure = new_part.departure + ff*increment;
%     new_part.start = new_part.start + ff*increment;
%     new_part.finish = new_part.finish + ff*increment;
%     new_part.train_id = new_part.train_id+ff*100;
%     
%     
%     BHP = [BHP; new_part];
%     
%     ff = ff+1;
% end
% 
% % Depends on the frequency per type
% freq_IC_0 = settings.TT.frequency.IC0;
% increment = 3600/freq_IC_0;
% ff = 1;
% while ff < freq_IC_0
%     
%     % Which rows?
%     rows = find(strcmp(base_tt.train_type,'IC') & base_tt.direction == 0);
%     
%     new_part = base_tt(rows,:);
%     new_part.arrival = new_part.arrival + ff*increment;
%     new_part.departure = new_part.departure + ff*increment;
%     new_part.start = new_part.start + ff*increment;
%     new_part.finish = new_part.finish + ff*increment;
%     new_part.train_id = new_part.train_id+ff*100;
%     
%     
%     BHP = [BHP; new_part];
%     
%     ff = ff+1;
% end
% 
% % Sort the trains by line id
% BHP = sortrows(BHP,'train_id','ascend');
% BHP.event_id = [1:size(BHP,1)]';

%% Finally, duplicate all this for the frequency
nr_hours = settings.disruption.duration + settings.disruption.nr_hours_after;
hh = 1;

complete_tt = BHP;
increment = 3600;
while hh < nr_hours
    new_part = BHP;
    new_part.arrival = new_part.arrival + hh*increment;
    new_part.departure = new_part.departure + hh*increment;
    new_part.start = new_part.start + hh*increment;
    new_part.finish = new_part.finish + hh*increment;
    new_part.train_id = new_part.train_id+hh*1000;
    
    hh = hh+1;
    complete_tt = [complete_tt; new_part];
end
complete_tt.event_id = [1:size(complete_tt,1)]';

%% Return these ones
base_tt = basic_timetable;
hour_tt = BHP;
complete_tt = complete_tt;






end