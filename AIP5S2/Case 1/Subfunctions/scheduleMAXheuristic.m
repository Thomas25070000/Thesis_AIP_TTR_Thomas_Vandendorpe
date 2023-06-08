function [new_timetable, measures, statistics] = scheduleMAXheuristic(timetable, blocksections, settings)

% Construct the full timetable, including NEW start and finish times.
newTT = [];
trains = unique(timetable.train_id);
row = 0;
t_before = settings.TT.blocktimes.setup;
t_after = settings.TT.blocktimes.afterIC;
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
            [reg_app, ~] = returnFirstApproachTime(blocksections, settings, TT_train.train_type(1), mod(train,10));
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
new_dep = new_dep(I);

dir0 = find(mod(trains,10) == 0);
dir1 = find(mod(trains,10) == 1);

trains0 = trains(dir0);
release0 = release(dir0);
departure0 = departure(dir0);
new_dep0 = new_dep(dir0);

trains1 = trains(dir1);
release1 = release(dir1);
departure1 = departure(dir1);
new_dep1 = new_dep(dir1);


maxSameDir = settings.solver.maxHeuristicNr;
waitingTime = settings.solver.maxHeuristicWaiting;

tic

counterDir0 = 0;
counterDir1 = 0;
% Now, schedule all trains by delaying them in the timetable!
% If they would depart from the corridor after more than w_cancel, this
% results in cancellation.
tt = 1;
% First one is scheduled as regular, take the earliest release!
currentDir = mod(trains(1),10);
if currentDir == 0
    prev_scheduled = trains0(1);
    orig_delays0(1) = new_dep0(1) - departure0(1);
    new_delay0(1) = 0;
    timing = new_dep0(1);
    counterDir0 = 1;
else
    prev_scheduled = trains1(1);
    orig_delays1(1) = new_dep1(1) - departure1(1);
    new_delay1(1) = 0;
    timing = new_dep1(1);
    counterDir1 = 1;
end
sameDir = 1;
prev_rows = find(newTT.train_id == prev_scheduled);



newTT.cancelled(1) = 0;

while counterDir0 < length(trains0) && counterDir1 < length(trains1)
    
    if counterDir0 == 17 || counterDir1 == 17
        disp('ok');
    end
    
    

    
    % If we go over the limit, switch directions!
    if sameDir >= maxSameDir
        sameDir = 0;
        currentDir = ~currentDir;
        justSwitched = 1;
    end
    
    flagSchedule = 0;   % Do we have to schedule a train, or just switch directions?
    if currentDir == 0
        % Currently serving direcion 0.
        % Take the next train in direction 0, if it is already waiting!
        bufSwitch = 0;
        if settings.constraints.minbuffer
            % Other direction!
            bufSwitch = bufSwitch + settings.TT.headway.other * 60;
        end
        
        if release0(counterDir0+1) <= timing + waitingTime + bufSwitch
            % Take this one!
            sameDir = sameDir + 1;
            counterDir0 = counterDir0+1;
            this = trains0(counterDir0);
            flagSchedule = 1;
        else 
            currentDir = 1;
            sameDir = 0;
        end
    else
        bufSwitch = 0;
        if settings.constraints.minbuffer
            % Other direction!
            bufSwitch = bufSwitch + settings.TT.headway.other * 60;
        end
        
        if release1(counterDir1+1) <= timing + waitingTime + bufSwitch
            % Take this one!
            sameDir = sameDir + 1;
            counterDir1 = counterDir1+1;
            this = trains1(counterDir1);
            flagSchedule = 1;
        end
    end
        
    if flagSchedule
        this_rows = find(newTT.train_id == this);
        if mod(prev_scheduled,10) == mod(this,10)
            % Trains run in the same direction.
            % Minimum headway time?
            headway = max(max(newTT.finish(prev_rows) - newTT.start(this_rows)),0);
            if settings.constraints.minbuffer
            % Same direction!
                headway = headway + settings.TT.headway.bufferafter * 60;
            end
            newTT = delayTrain(newTT, this, headway, 1);
        else
            headway = max(newTT.finish(prev_rows(end)) - newTT.start(this_rows(1)),0);
            if settings.constraints.minbuffer
                % Other direction!
                headway = headway + settings.TT.headway.other * 60;
            end
            newTT = delayTrain(newTT, this, headway, 1);
        end

        % Determine delay
        exit = newTT.departure(this_rows(end));
        if currentDir == 0
            delay = exit - departure0(counterDir0);
            if delay >= settings.weights.cancel
                % If the resulting delay is too high, cancelled it.
                newTT.cancelled(this_rows) = 1;
                orig_delays0(counterDir0) = 0;
                new_delay0(counterDir0) = 0;
                sameDir = sameDir - 1;
            else
                % Otherwise, this is the previous train for the next iteration.
                prev_scheduled = this;
                prev_rows = find(newTT.train_id == prev_scheduled);
                orig_delays0(counterDir0) = delay;
                new_delay0(counterDir0) = exit - new_dep0(counterDir0);
                timing = exit;
            end
        else
            delay = exit - departure1(counterDir1);
            if delay >= settings.weights.cancel
                % If the resulting delay is too high, cancelled it.
                newTT.cancelled(this_rows) = 1;
                orig_delays1(counterDir1) = 0;
                new_delay1(counterDir1) = 0;
                sameDir = sameDir - 1;
            else
                % Otherwise, this is the previous train for the next iteration.
                prev_scheduled = this;
                prev_rows = find(newTT.train_id == prev_scheduled);
                orig_delays1(counterDir1) = delay;
                new_delay1(counterDir1) = exit - new_dep1(counterDir1);
                timing = exit;
            end
        end
    else
        % We may have to progress the time, as we are switching directions!
        if sameDir >= maxSameDir
            % Definitely switch!
            sameDir = 0;
            if currentDir == 1
                timing = release0(counterDir0+1);
            else
                timing = release1(counterDir1+1);
            end
        else
            switch currentDir
                case 0
                    if release1(counterDir1+1) < release0(counterDir0+1)
                        % Switch direction!
                        currentDir = 1;
                        sameDir = 0;
                        timing = release1(counterDir1+1);
                    else
                        % Keep this direction!
                        timing = release0(counterDir0+1);
                    end
                case 1
                    if release1(counterDir1+1) < release0(counterDir0+1)
                        % Keep this direction!
                        timing = release1(counterDir1+1);
                    else
                        % Switch direction!
                        currentDir = 0;
                        sameDir = 0;
                        timing = release0(counterDir0+1);
                    end
            end
        end
        if currentDir == 1
            timing = release1(counterDir1+1);
        else
            timing = release0(counterDir0+1);
        end
    end
end

% Are there still unscheduled ones? Add them!
while counterDir0 < length(trains0)
    % Schedule the next train!
    counterDir0 = counterDir0 + 1;
    this = trains0(counterDir0);
    this_rows = find(newTT.train_id == this);
    if mod(prev_scheduled,10) == mod(this,10)
        % Trains run in the same direction.
        % Minimum headway time?
        headway = max(max(newTT.finish(prev_rows) - newTT.start(this_rows)),0);
        newTT = delayTrain(newTT, this, headway, 1);
    else
        headway = max(newTT.finish(prev_rows(end)) - newTT.start(this_rows(1)),0);
        newTT = delayTrain(newTT, this, headway, 1);
    end

    % Determine delay
    exit = newTT.departure(this_rows(end));
    delay = exit - departure0(counterDir0);
    if delay >= settings.weights.cancel
        % If the resulting delay is too high, cancelled it.
        newTT.cancelled(this_rows) = 1;
        orig_delays0(counterDir0) = 0;
        new_delay0(counterDir0) = 0;
    else
        % Otherwise, this is the previous train for the next iteration.
        prev_scheduled = this;
        prev_rows = find(newTT.train_id == prev_scheduled);
        orig_delays0(counterDir0) = delay;
        new_delay0(counterDir0) = exit - new_dep0(counterDir0);
        timing = exit;
    end
end
while counterDir1 < length(trains1)
    % Schedule the next train!
    counterDir1 = counterDir1 + 1;
    this = trains1(counterDir1);
    this_rows = find(newTT.train_id == this);
    if mod(prev_scheduled,10) == mod(this,10)
        % Trains run in the same direction.
        % Minimum headway time?
        headway = max(max(newTT.finish(prev_rows) - newTT.start(this_rows)),0);
        newTT = delayTrain(newTT, this, headway, 1);
    else
        headway = max(newTT.finish(prev_rows(end)) - newTT.start(this_rows(1)),0);
        newTT = delayTrain(newTT, this, headway, 1);
    end

    % Determine delay
    exit = newTT.departure(this_rows(end));
    delay = exit - departure1(counterDir1);
    if delay >= settings.weights.cancel
        % If the resulting delay is too high, cancelled it.
        newTT.cancelled(this_rows) = 1;
        orig_delays1(counterDir1) = 0;
        new_delay1(counterDir1) = 0;
    else
        % Otherwise, this is the previous train for the next iteration.
        prev_scheduled = this;
        prev_rows = find(newTT.train_id == prev_scheduled);
        orig_delays1(counterDir1) = delay;
        new_delay1(counterDir1) = exit - new_dep1(counterDir1);
        timing = exit;
    end
end
    
time = toc;

% Merge the data on delays, trains, ...
orig_delays = zeros(Ntrains,1);
orig_delays(find(ismember(trains,trains0))) = orig_delays0';
orig_delays(find(ismember(trains,trains1))) = orig_delays1';
new_delay = zeros(Ntrains,1);
new_delay(find(ismember(trains,trains0))) = new_delay0';
new_delay(find(ismember(trains,trains1))) = new_delay1';


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
        for ee = 1:length(rows_newTT)
            if newTT.stops(rows_newTT(ee))
                % Which penalty?
                bs = newTT.bs(rows_newTT(ee));
                ss = find([stops.bsID] == bs);
                timetable.stopskipped(rowsTT(ee)) = 0;
                if mod(newTT.train_id(rows_newTT(ee)),10) == 1
                    timetable.skpenalty(rowsTT(ee)) = settings.weights.cancel * ...
                        stops(ss).fractionpax_dir1 * stops(ss).weight_dir1;
                else
                    timetable.skpenalty(rowsTT(ee)) = settings.weights.cancel * ...
                        stops(ss).fractionpax_dir0 * stops(ss).weight_dir0;
                end
            end
        end

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
directions = 1;
%     directions = 0;
[line_new, blocks_new] = plotTT_new(timetable, blocksections, settings, type, include_blocks, directions);
    
new_timetable = timetable;






end