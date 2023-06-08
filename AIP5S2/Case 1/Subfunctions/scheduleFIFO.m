function [new_timetable, measures, statistics] = scheduleFIFO(timetable, blocksections, settings)

% Construct the full timetable, including NEW start and finish times.
newTT = [];
trains = unique(timetable.train_id);
row = 0;
t_before = settings.TT.blocktimes.setup;
t_after = settings.TT.blocktimes.afterIC;
disruption_dir = settings.disruption.direction;
crossing_threshold = 0;
for tt = 1:length(trains)
    train = trains(tt);
    TT_train = timetable(find(timetable.train_id == train),:);
    if floor(train/1000) >= settings.disruption.duration
        afterend = 1;
    else
        afterend = 0;
    end
    
    % First row starts at 0.
    for ee = 1:size(TT_train,1)
        row = row + 1;
        newTT(row).train_id = train;
        newTT(row).bs = TT_train.blocksection(ee);
        if ee == 1
            newTT(row).arrival = TT_train.arrival(1);
            [reg_app, ~] = returnFirstApproachTime(blocksections, settings, TT_train.train_type(1), mod(train,100));
        else
            newTT(row).arrival = newTT(row-1).departure;
            reg_app = newTT(row-1).running;
            app = newTT(row-1).running;
        end
        newTT(row).running = TT_train.running(ee);
        newTT(row).departure = newTT(row).arrival + newTT(row).running;
        newTT(row).start = newTT(row).arrival - reg_app - t_before;
        newTT(row).finish = newTT(row).departure + t_after;
    end
end
newTT = struct2table(newTT);

trains = unique(newTT.train_id);
Ntrains = length(trains);
release = zeros(Ntrains,1);
% Get the release times.
for tt = 1:Ntrains
    ev = newTT(find(newTT.train_id == trains(tt)),:);
    old_ev = timetable(find(timetable.train_id == trains(tt)),:);
    release(tt) = ev.arrival(1);
    departure(tt) = old_ev.departure(end);
    new_dep(tt) = ev.departure(end);
end
% Order trains according to release time.
[release,I] = sortrows(release);
trains = trains(I);
departure = departure(I);
direction = mod(trains,100);
new_dep = new_dep(I);

tic




% Now, schedule all trains by delaying them in the timetable!
% If they would depart from the corridor after more than w_cancel, this
% results in cancellation.
tt = 1;
% First one is schedule as regular!
prev_scheduled = trains(tt);
orig_delays(tt) = new_dep(tt) - departure(tt);
new_delay(tt) = 0;
tt = 2;
newTT.cancelled(1) = 0;
last_common_block = 3;
default_from_sander = 12 + 9;


% Now, schedule all trains by delaying them in the timetable!
% If they would depart from the corridor after more than w_cancel, this
% results in cancellation.
tt = 1;
% First one is schedule as regular!
prev_scheduled = trains(tt);
orig_delays(tt) = new_dep(tt) - departure(tt);
new_delay(tt) = 0;
tt = 2;
newTT.cancelled(1) = 0;
last_common_block = 3;
default_from_sander = 12 + 9;
default_from_sander = 0;


while tt <= Ntrains
    this = trains(tt);
    this_rows = find(newTT.train_id == this);
    for prev_index = 1:tt-1
        prev_scheduled = trains(prev_index);
        prev_rows = find(newTT.train_id == prev_scheduled);
        if newTT.cancelled(prev_rows(1)) == 0
            if mod(prev_scheduled, 100) == mod(this, 100)
                % Trains run in the same direction.
                % Minimum headway time?
                headway = max(max(newTT.finish(prev_rows) - newTT.start(this_rows)),0);
                newTT = delayTrain(newTT, this, headway, 1);
            end
            if mod(prev_scheduled, 100)+mod(this, 100) == 5        
                headway = max(max(newTT.finish(prev_rows(end-last_common_block:end)) - newTT.start(this_rows(end-last_common_block:end))),0);
                newTT = delayTrain(newTT, this, headway, 1);
            end
            if mod(prev_scheduled, 100)+mod(this, 100) == 25        
                if newTT.start(this_rows(last_common_block))-newTT.finish(prev_rows(last_common_block))<crossing_threshold
                    extra_delay = crossing_threshold - (newTT.start(this_rows(last_common_block))-newTT.finish(prev_rows(last_common_block)));
                    newTT = delayTrain(newTT,this,extra_delay,1);
                end
            end
            if abs(mod(prev_scheduled, 100)-mod(this, 100)) == 10 && (mod(prev_scheduled,100) == mod(disruption_dir,100) || mod(this,100) == mod(disruption_dir,100))
                headway = max(newTT.finish(prev_rows(end)) - newTT.start(this_rows(1))+default_from_sander,0);
                newTT = delayTrain(newTT, this, headway, 1);
            end
            if abs(mod(prev_scheduled, 100)-mod(this, 100)) == 10 && (mod(prev_scheduled,10) ~= mod(disruption_dir,10) || mod(this,100) ~= mod(disruption_dir,100))
                headway = 0;
                newTT = delayTrain(newTT, this, headway, 1);
            end    
            if abs(mod(prev_scheduled, 100)-mod(this, 100)) == 11
                if mod(this, 100) == 2  
                    if newTT.start(this_rows(end-last_common_block+1))-newTT.finish(prev_rows(last_common_block))<crossing_threshold
                        extra_delay = crossing_threshold - (newTT.start(this_rows(end-last_common_block+1))-newTT.finish(prev_rows(last_common_block)));
                        newTT = delayTrain(newTT,this,extra_delay,1);
                    end
                end
                if mod(this, 100) == 13  
                    if newTT.start(this_rows(last_common_block))-newTT.finish(prev_rows(end-last_common_block+1))<crossing_threshold
                        extra_delay = crossing_threshold - (newTT.start(this_rows(last_common_block))-newTT.finish(prev_rows(end-last_common_block+1)));
                        newTT = delayTrain(newTT,this,extra_delay,1);
                    end
                end
            end
            if abs(mod(prev_scheduled, 100)-mod(this, 100)) == 9 && mod(prev_scheduled,100)==03
                headway = max(newTT.finish(prev_rows(end)) - newTT.start(this_rows(1))+ default_from_sander,0) ;
                newTT = delayTrain(newTT, this, headway, 1);
            end
            if abs(mod(prev_scheduled, 100)-mod(this, 100)) == 9 && mod(prev_scheduled,100)==12
                headway = max(newTT.finish(prev_rows(last_common_block)) - newTT.start(this_rows(end-last_common_block+1)) + default_from_sander,0);
                newTT = delayTrain(newTT, this, headway, 1);
            end        
    %         if (abs(mod(prev_scheduled, 100)-mod(this, 100)) == 9 || abs(mod(prev_scheduled, 100)-mod(this, 100)) == 11) ...
    %                 && (mod(prev_scheduled, 100)==12 || mod(prev_scheduled, 100) ==13)
    %             headway = max(newTT.finish(prev_rows(last_common_block)) - newTT.start(this_rows(end-last_common_block+1)),0);
    %             newTT = delayTrain(newTT, this, headway, 1);
    %         end
            % Determine delay
            exit = newTT.departure(this_rows(end));
            delay = exit - departure(tt);
            orig_delays(tt) = delay;
            new_delay(tt) = exit - new_dep(tt);
            if delay >= settings.weights.cancel
                % If the resulting delay is too high, cancelled it.
                newTT.cancelled(this_rows) = 1;
                orig_delays(tt) = 0;
                new_delay(tt) = 0;
            end
%             % Otherwise, this is the previous train for the next iteration.
%             prev_scheduled = this;
%             orig_delays(tt) = delay;
%             new_delay(tt) = exit - new_dep(tt);
%         end
        end
    end

    tt = tt+1;
    
end

time = toc;
settings.FIFOtime = time;

% Collect everything in the full timetable.
timetable.cancelled(1) = 0;
timetable.adjusted_arrival(1) = 0;
timetable.adjusted_departure(1) = 0;
timetable.adjusted_start(1) = 0;
timetable.adjusted_finish(1) = 0;
timetable.deviated(1) = 0;
cancelled = zeros(Ntrains,1);
for tt = 1:Ntrains
    train = trains(tt);
    
    rowsTT = find(timetable.train_id == train);
    rows_newTT = find(newTT.train_id == train);
    
    if newTT.cancelled(rows_newTT)
        cancelled(tt) = 1;
    end
    timetable.cancelled(rowsTT) = newTT.cancelled(rows_newTT);
    timetable.adjusted_arrival(rowsTT) = newTT.arrival(rows_newTT);
    timetable.adjusted_departure(rowsTT) = newTT.departure(rows_newTT);
    timetable.adjusted_start(rowsTT) = newTT.start(rows_newTT);
    timetable.adjusted_finish(rowsTT) = newTT.finish(rows_newTT);

end

[measures, statistics] = deriveMeasuresFIFO(timetable, cancelled, settings, orig_delays, new_delay)

type = 'complete';
include_blocks = 1;
directions = [02 03 12 13];
%     directions = 0;
%[line_new, blocks_new] = plotTT_new_2(timetable, blocksections, settings, type, include_blocks, directions);
%[line_new, blocks_new] = plotTT_new_3(timetable, blocksections, settings, type, include_blocks, directions);
[line_new, blocks_new] = plotTT_full(timetable, blocksections, settings, type, include_blocks, directions);


new_timetable = timetable;






end