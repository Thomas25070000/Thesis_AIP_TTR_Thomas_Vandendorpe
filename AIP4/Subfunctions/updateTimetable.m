function new_timetable = updateTimetable(timetable, blocksections, settings, trains, cancelled)

Ntrains = length(trains);

% Update the other events that happen on the closed section!
for tt = 1:Ntrains
    train = trains(tt);

    if cancelled(tt) == 0
        rows_train = find(timetable.train_id == train);

        % Only update after it hits the first closed part.
        from_closed_on = 0;
        for rr = 1:length(rows_train)
            % What block are we on?
            bb = timetable.blocksection(rows_train(rr));

    
                if blocksections.closed(bb) || from_closed_on
                    from_closed_on = 1;
    
                    if timetable.adjusted_arrival(rows_train(rr)) == 0
                        if timetable.arrival(rows_train(rr))==0
                            timetable.adjusted_arrival(rows_train(rr)) = 0;
                        else
                        % It has not yet been assigned!
                            timetable.adjusted_arrival(rows_train(rr)) = timetable.adjusted_departure(rows_train(rr-1));
                        end
                    end
                    if timetable.adjusted_departure(rows_train(rr)) == 0
                        % It has not yet been assigned!
                        timetable.adjusted_departure(rows_train(rr)) = ...
                            timetable.adjusted_arrival(rows_train(rr)) + timetable.running(rows_train(rr));
                    end
                else
                    timetable.adjusted_arrival(rows_train(rr)) = timetable.arrival(rows_train(rr));
                    timetable.adjusted_departure(rows_train(rr)) = ...
                            timetable.arrival(rows_train(rr)) + timetable.running(rows_train(rr));
                end

            % Update start and finish times
            t_setup = settings.TT.blocktimes.setup;
            t_release_R = settings.TT.blocktimes.afterR;
            t_release_IC = settings.TT.blocktimes.afterIC;

            % If there is no approach, we assume the same running time as in the first
            % block section!
            for ee = 1:size(timetable,1)
                % What is the approach time?
                try
                    if timetable.train_id(ee) ~= timetable.train_id(ee-1)
                        % This is the first action of the train.
                        t_before = timetable.adjusted_departure(ee) - timetable.adjusted_arrival(ee) + t_setup;
                    else
                        % There has already been an event
                        t_before = timetable.adjusted_departure(ee-1) - timetable.adjusted_arrival(ee-1) + t_setup;
                    end
                catch
                    % This is the first event of the timetable!
                    t_before = timetable.adjusted_departure(ee) - timetable.adjusted_arrival(ee) + t_setup;
                end
                timetable.adjusted_start(ee) = timetable.adjusted_arrival(ee) - t_before;

                switch timetable.train_type{ee}
                    case 'IC'
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
        if cancelled(find(trains == tt_id(tt))) == 0
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