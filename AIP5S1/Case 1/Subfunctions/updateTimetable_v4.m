function new_timetable = updateTimetable_v4(timetable, blocksections, settings, trains, cancelled)

Ntrains = length(trains);

timetable.adjusted_start(1) = 0;
timetable.adjusted_finish(1) = 0;

% Approach time if they don't get stopped.
Nblocks = size(blocksections,1);
approachtime = zeros(Ntrains,Nblocks);
closedblocks = find(blocksections.closed);

for tt = 1:Ntrains
    rows_train = find(timetable.train_id == trains(tt));
    direction(tt) = timetable.direction(rows_train(1));
    traintype = timetable.train_type{rows_train(1)};
    deviated(tt) = timetable.deviated(rows_train(1));
    
    if direction(tt) == 12
        
        approachtime(tt,1) = returnFirstApproachTime(blocksections, settings, traintype, mod(direction(tt),100));
        
        for bb = 2:Nblocks
            try
                approachtime(tt,bb) = timetable.running(find(timetable.train_id == trains(tt)...
                                        & timetable.blocksection==bb-1));
            end
            if (direction(tt) == settings.disruption.direction) && ~settings.disruption.signals
                try
                    if (bb-1) == closedblocks(end)
                        approachtime(tt,bb) = sum(timetable.running(find(...
                                    timetable.train_id == trains(tt) & blocksections.closed(timetable.blocksection))));
                    end
                end
            end
        end
     end
     if direction(tt) == 13
        
        approachtime(tt,1) = returnFirstApproachTime(blocksections, settings, traintype, mod(direction(tt),100));
        
        for bb = 2:Nblocks
            try
                approachtime(tt,bb) = timetable.running(find(timetable.train_id == trains(tt)...
                                        & timetable.blocksection==bb-1));
            end
            if (direction(tt) == settings.disruption.direction) && ~settings.disruption.signals
                try
                    if (bb-1) == closedblocks(end)
                        approachtime(tt,bb) = sum(timetable.running(find(...
                                    timetable.train_id == trains(tt) & blocksections.closed(timetable.blocksection))));
                    end
                end
            end
        end
     end
    if direction(tt) == 02
         approachtime(tt,end) = returnFirstApproachTime(blocksections, settings, traintype, mod(direction(tt),100));
        
        
        for bb = (Nblocks-1):-1:1
            try
                approachtime(tt,bb) = timetable.running(find(timetable.train_id == trains(tt)...
                                        & timetable.blocksection==bb+1));
            end
            if (direction(tt) == settings.disruption.direction) && ~settings.disruption.signals
                try
                    if (bb+1) == closedblocks(1)
                        approachtime(tt,bb) = sum(timetable.running(find(...
                                    timetable.train_id == trains(tt) & blocksections.closed(timetable.blocksection))));
                    end
                end
            end
        end
    end
    if direction(tt) == 03
         approachtime(tt,end) = returnFirstApproachTime(blocksections, settings, traintype, mod(direction(tt),100));
        
        
        for bb = (Nblocks-1):-1:1
            try
                approachtime(tt,bb) = timetable.running(find(timetable.train_id == trains(tt)...
                                        & timetable.blocksection==bb+1));
            end
            if (direction(tt) == settings.disruption.direction) && ~settings.disruption.signals
                try
                    if (bb+1) == closedblocks(1)
                        approachtime(tt,bb) = sum(timetable.running(find(...
                                    timetable.train_id == trains(tt) & blocksections.closed(timetable.blocksection))));
                    end
                end
            end
        end
    end
end


% Update the other events that happen on the closed section!
for tt = 1:Ntrains
    train = trains(tt);

    if cancelled(tt) == 0 && deviated(tt) == 0
        rows_train = find(timetable.train_id == train);

        % Only update after it hits the first closed part.
        from_closed_on = 0;
        for rr = 1:length(rows_train)
            % What block are we on?
            bb = timetable.blocksection(rows_train(rr));

            if blocksections.closed(bb) || from_closed_on
                from_closed_on = 1;

                if timetable.adjusted_arrival(rows_train(rr)) == 0 && ~(timetable.arrival(rows_train(rr)) ==0)
                    % It has not yet been assigned!
                    timetable.adjusted_arrival(rows_train(rr)) = timetable.adjusted_departure(rows_train(rr-1));
                end
%                 if timetable.adjusted_departure(rows_train(rr)) == 0
                    % It has not yet been assigned!
                    timetable.adjusted_departure(rows_train(rr)) = ...
                        timetable.adjusted_arrival(rows_train(rr)) + timetable.running(rows_train(rr));
%                 end
            end
        end
        
        % Update start and finish times
        t_setup = settings.TT.blocktimes.setup;
        t_release_R = settings.TT.blocktimes.afterR;
        t_release_IC = settings.TT.blocktimes.afterIC;

        % If there is no approach, we assume the same running time as in the first
        % block section!
        for rr = 1:length(rows_train)
            ee = timetable.event_id(rows_train(rr));
            % What is the approach time?
            
%             if ee == 310
%                 disp('check');
%             end
            
            
%             try
%                 if timetable.train_id(ee) ~= timetable.train_id(ee-1)
%                     % This is the first action of the train.
% %                     t_before = timetable.adjusted_departure(ee) - timetable.adjusted_arrival(ee) + t_setup;
%                     t_before = timetable.running(ee) + t_setup;
%                 else
%                     % There has already been an event
%                     t_before = timetable.adjusted_departure(ee-1) - timetable.adjusted_arrival(ee-1) + t_setup;
%                 end
%             catch
%                 % This is the first event of the timetable!
%                 t_before = timetable.adjusted_departure(ee) - timetable.adjusted_arrival(ee) + t_setup;
%             end
            bs = timetable.blocksection(rows_train(rr));
            timetable.adjusted_start(ee) = timetable.adjusted_arrival(ee) ...
                            - t_setup - approachtime(tt,bs);
            
            try
                if timetable.train_id(ee) == timetable.train_id(ee+1)
                    switch timetable.train_type{ee}
                        case {'IC','THA'}
                            timetable.adjusted_finish(ee) = timetable.adjusted_arrival(ee+1) + t_release_IC;
                        case 'R'
                            timetable.adjusted_finish(ee) = timetable.adjusted_arrival(ee+1)+ t_release_R;
                    end
                else
                    switch timetable.train_type{ee}
                        case {'IC','THA'}
                            timetable.adjusted_finish(ee) = timetable.adjusted_departure(ee) + t_release_IC;
                        case 'R'
                            timetable.adjusted_finish(ee) = timetable.adjusted_departure(ee)+ t_release_R;
                    end
                end
            catch
                % If it fails, then we are at the end of all trains!
                switch timetable.train_type{ee}
                        case {'IC','THA'}
                            timetable.adjusted_finish(ee) = timetable.adjusted_departure(ee) + t_release_IC;
                        case 'R'
                            timetable.adjusted_finish(ee) = timetable.adjusted_departure(ee)+ t_release_R;
                end
            end
        end
    end
end

% In case we do not have signals in the closed direction, the starting and
% finish times have to be updated!
if settings.disruption.signals == 0
    % For the trains running in the closed direction, we have a problem.
    affected = timetable(find((timetable.direction == settings.disruption.direction) ...
                        & blocksections.closed(timetable.blocksection)),:);
    
    % Which train ids?
    tt_id = unique(affected.train_id);
    
    for tt = 1:length(tt_id)
        if cancelled(find(trains == tt_id(tt))) == 0 && deviated(find(trains == tt_id(tt)))
            rows_orig = find((timetable.train_id == tt_id(tt)) ...
                                & blocksections.closed(timetable.blocksection));
            start = timetable.adjusted_start(rows_orig(1));
            finish = timetable.adjusted_finish(rows_orig(end));
            timetable.adjusted_start(rows_orig) = start;
            timetable.adjusted_finish(rows_orig) = finish;
        end
    end

end



new_timetable = timetable;

end