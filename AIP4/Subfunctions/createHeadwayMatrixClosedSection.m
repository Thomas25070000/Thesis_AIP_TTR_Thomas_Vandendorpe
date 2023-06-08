function [minHW, trains] = createHeadwayMatrixClosedSection(timetable, blocksections, settings)

orig_tt = timetable;



trains = unique(timetable.train_id);

% Adjust the event times. Let all start at time 0!
timetable.arrival(1) = 0;
timetable.departure(1) = timetable.running(1);
for ee = 2:size(timetable,1)
    if timetable.train_id(ee) ~= timetable.train_id(ee-1)
        % We have a new train!
        timetable.arrival(ee) = 0;
    else
        timetable.arrival(ee) = timetable.departure(ee-1);
    end
    timetable.departure(ee) = timetable.arrival(ee) + timetable.running(ee);
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
            t_before = timetable.departure(ee) - timetable.arrival(ee) + t_setup;
        else
            % There has already been an event
            t_before = timetable.departure(ee-1) - timetable.arrival(ee-1) + t_setup;
        end
    catch
        % This is the first event of the timetable!
        t_before = timetable.departure(ee) - timetable.arrival(ee) + t_setup;
    end
    timetable.start(ee) = timetable.arrival(ee) - t_before;
    
    switch timetable.train_type{ee}
        case 'IC'
            timetable.finish(ee) = timetable.departure(ee) + t_release_IC;
        case 'R'
            timetable.finish(ee) = timetable.departure(ee)+ t_release_R;
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
        rows_orig = find((timetable.train_id == tt_id(tt)) ...
                            & blocksections.closed(timetable.blocksection));
        start = timetable.start(rows_orig(1));
        finish = timetable.finish(rows_orig(end));
        timetable.start(rows_orig) = start;
        timetable.finish(rows_orig) = finish;
    end

end


% Using these 'patterns', the min. headways in order to achieve separation
% can be calculated.

% Do this between each pair of train ids
trains = unique(timetable.train_id)';
direction = zeros(1,length(trains));
types = [];

% Remove the events outside the closed sections.
timetable(find(~blocksections.closed(timetable.blocksection)),:) = [];



for tt = 1:length(trains)
    ev = timetable(find(timetable.train_id == trains(tt)),:);
    direction(tt) = ev.direction(1);
    types{tt} = ev.train_type{1};
    
    % Put all first arrivals back to 0
    first_arrival = ev.arrival(1);
    timetable = delayTrain(timetable, trains(tt), -first_arrival, 1);
end

for ii = 1:length(trains)-1
    for jj = (ii+1):length(trains)
        % lambda equal to 1
        if direction(ii) == direction(jj)
            % Trains are running in the same direction. Regardless of the
            % train type (they are running at the same speed), we get the
            % same headway!).
            rows_ii = find(timetable.train_id == trains(ii));
            rows_jj = find(timetable.train_id == trains(jj));
            time_ii_jj = timetable.finish(rows_ii) - timetable.start(rows_jj);
            min_HW_ii_jj = max(time_ii_jj);
            time_jj_ii = timetable.finish(rows_jj) - timetable.start(rows_ii);
            min_HW_jj_ii = max(time_jj_ii);
            
            % If norms are to be used!
            % Timetable only has those running over the closed sections.
            norm_ii_jj = settings.TT.headway.same * 60;
            norm_jj_ii = settings.TT.headway.same * 60;
            
            minHW(ii,jj,2) = max(min_HW_ii_jj, norm_ii_jj);
            minHW(jj,ii,2) = max(min_HW_jj_ii, norm_jj_ii);
            
        else
            
            rows_ii = find(timetable.train_id == trains(ii));
            rows_jj = find(timetable.train_id == trains(jj));
            
            run_ii = sum(timetable.running(rows_ii));
            run_jj = sum(timetable.running(rows_jj));
            
            
            
        
            switch types{ii}
                case 'IC'
                    clear_ii = t_release_IC;
                case 'R'
                    clear_ii = t_release_R;
            end
            switch types{jj}
                case 'IC'
                    clear_jj = t_release_IC;
                case 'R'
                    clear_jj = t_release_R;
            end
        
            before_ii = timetable.arrival(rows_ii(1)) - timetable.start(rows_ii(1));
            before_jj = timetable.arrival(rows_jj(1)) - timetable.start(rows_jj(1));
            
            % If norms are to be used!
            % Timetable only has those running over the closed sections.            
            norm_ii_jj = 0;
            norm_jj_ii = 0;
%             norm_ii_jj = run_ii + settings.TT.headway.other * 60;
%             norm_jj_ii = run_jj + settings.TT.headway.other * 60;
            
            minHW(ii,jj,2) = max(run_ii + clear_ii + before_jj, norm_ii_jj);
            minHW(jj,ii,2) = max(run_jj + clear_jj + before_ii, norm_jj_ii);
        
        end
        % lambda equal to 0
        if direction(ii) == direction(jj)
            % Trains are running in the same direction. Regardless of the
%             % train type (they are running at the same speed), we get the
%             % same headway!).
%             rows_ii = find(timetable.train_id == trains(ii));
%             rows_jj = find(timetable.train_id == trains(jj));
%             time_ii_jj = timetable.finish(rows_ii) - timetable.start(rows_jj);
%             min_HW_ii_jj = max(time_ii_jj);
%             time_jj_ii = timetable.finish(rows_jj) - timetable.start(rows_ii);
%             min_HW_jj_ii = max(time_jj_ii);
%             
%             % If norms are to be used!
%             % Timetable only has those running over the closed sections.
%             norm_ii_jj = settings.TT.headway.same * 60;
%             norm_jj_ii = settings.TT.headway.same * 60;

            % not considere for now
            

            %minHW(ii,jj,1) = max(min_HW_ii_jj, norm_ii_jj);
            %minHW(jj,ii,1) = max(min_HW_jj_ii, norm_jj_ii);
            minHW(ii,jj,1) = -1000;
            minHW(jj,ii,1) = -1000;
            
        else
            
            rows_ii = find(timetable.train_id == trains(ii));
            rows_jj = find(timetable.train_id == trains(jj));
            
            run_ii_first = sum(timetable.running(rows_ii(1)));
            run_jj_first = sum(timetable.running(rows_jj(1)));
            run_ii_last = sum(timetable.running(rows_ii(end)));
            run_jj_last = sum(timetable.running(rows_jj(end)));
            run_ii_all = sum(timetable.running(rows_ii));
            run_jj_all = sum(timetable.running(rows_jj));
            
            
        
            switch types{ii}
                case 'IC'
                    clear_ii = t_release_IC;
                case 'R'
                    clear_ii = t_release_R;
            end
            switch types{jj}
                case 'IC'
                    clear_jj = t_release_IC;
                case 'R'
                    clear_jj = t_release_R;
            end
        
            approach_ii = timetable.arrival(rows_ii(end)) - timetable.start(rows_ii(end));
            approach_jj = timetable.arrival(rows_jj(end)) - timetable.start(rows_jj(end));
            
            % If norms are to be used!
            % Timetable only has those running over the closed sections.            
            norm_ii_jj = 0;
            norm_jj_ii = 0;
%             norm_ii_jj = run_ii + settings.TT.headway.other * 60;
%             norm_jj_ii = run_jj + settings.TT.headway.other * 60;
            
            minHW(ii,jj,1) = run_ii_all-clear_ii-approach_ii-run_ii_last-run_jj_first;
            minHW(jj,ii,1) = run_jj_all-clear_jj-approach_jj-run_jj_last-run_ii_first;
        
        end
        % lambda equal to 0
        if direction(ii) == direction(jj)
            % Trains are running in the same direction. Regardless of the
%             % train type (they are running at the same speed), we get the
%             % same headway!).
%             rows_ii = find(timetable.train_id == trains(ii));
%             rows_jj = find(timetable.train_id == trains(jj));
%             time_ii_jj = timetable.finish(rows_ii) - timetable.start(rows_jj);
%             min_HW_ii_jj = max(time_ii_jj);
%             time_jj_ii = timetable.finish(rows_jj) - timetable.start(rows_ii);
%             min_HW_jj_ii = max(time_jj_ii);
%             
%             % If norms are to be used!
%             % Timetable only has those running over the closed sections.
%             norm_ii_jj = settings.TT.headway.same * 60;
%             norm_jj_ii = settings.TT.headway.same * 60;

            % not considere for now
            

            %minHW(ii,jj,1) = max(min_HW_ii_jj, norm_ii_jj);
            %minHW(jj,ii,1) = max(min_HW_jj_ii, norm_jj_ii);
            minHW(ii,jj,3) = -1000;
            minHW(jj,ii,3) = -1000;
            
        else
            
            rows_ii = find(timetable.train_id == trains(ii));
            rows_jj = find(timetable.train_id == trains(jj));
            
            run_ii_first = sum(timetable.running(rows_ii(1)));
            run_jj_first = sum(timetable.running(rows_jj(1)));
            run_ii_last = sum(timetable.running(rows_ii(end)));
            run_jj_last = sum(timetable.running(rows_jj(end)));
            run_ii_all = sum(timetable.running(rows_ii));
            run_jj_all = sum(timetable.running(rows_jj));
            
            
        
            switch types{ii}
                case 'IC'
                    clear_ii = t_release_IC;
                case 'R'
                    clear_ii = t_release_R;
            end
            switch types{jj}
                case 'IC'
                    clear_jj = t_release_IC;
                case 'R'
                    clear_jj = t_release_R;
            end
        
            approach_ii = timetable.arrival(rows_ii(end)) - timetable.start(rows_ii(end));
            approach_jj = timetable.arrival(rows_jj(end)) - timetable.start(rows_jj(end));
            
            % If norms are to be used!
            % Timetable only has those running over the closed sections.            
            norm_ii_jj = 0;
            norm_jj_ii = 0;
%             norm_ii_jj = run_ii + settings.TT.headway.other * 60;
%             norm_jj_ii = run_jj + settings.TT.headway.other * 60;
            
            minHW(ii,jj,3) = run_ii_first+run_jj_last+approach_jj-run_jj_all+clear_ii;
            minHW(jj,ii,3) = run_jj_first+run_ii_last+approach_ii-run_ii_all+clear_jj;
        
        end
    end
end


end