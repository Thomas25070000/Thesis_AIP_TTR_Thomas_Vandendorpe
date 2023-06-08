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
%timetable(find(~blocksections.closed(timetable.blocksection)),:) = [];



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
        if abs(direction(ii)-direction(jj))<2
            % Trains are running in the same direction. Regardless of the
            % train type (they are running at the same speed), we get the
            % same headway!).
            rows_ii = find(timetable.train_id == trains(ii)& ismember(timetable.blocksection, blocksections.id(blocksections.closed==1)));
            rows_jj = find(timetable.train_id == trains(jj)& ismember(timetable.blocksection, blocksections.id(blocksections.closed==1)));
            time_ii_jj = timetable.finish(rows_ii) - timetable.start(rows_jj);
            min_HW_ii_jj = max(time_ii_jj);
            release_jj_ii = timetable.arrival(rows_jj(1))-timetable.arrival(rows_ii(1));
            time_jj_ii = timetable.finish(rows_jj) - timetable.start(rows_ii);
            min_HW_jj_ii = max(time_jj_ii);
            release_ii_jj = timetable.arrival(rows_ii(1))-timetable.arrival(rows_jj(1));

            
            % If norms are to be used!
            % Timetable only has those running over the closed sections.
            norm_ii_jj = settings.TT.headway.same * 60;
            norm_jj_ii = settings.TT.headway.same * 60;
            
            minHW(ii,jj) = max(min_HW_ii_jj, norm_ii_jj)+release_jj_ii;
            minHW(jj,ii) = max(min_HW_jj_ii, norm_jj_ii)+release_ii_jj;
            
        end
        default_from_sander = 12 + 9;
        %default_from_sander = 0;
        if abs(direction(ii)-direction(jj))==10
            if direction(ii)>10
                rows_ii = find(timetable.train_id == trains(ii));
                t_all_ii = sum(timetable.running(rows_ii));
                rows_jj_open = find(timetable.train_id == trains(jj)& ismember(timetable.blocksection, blocksections.id(blocksections.closed==0)));
                t_open_jj = sum(timetable.running(rows_jj_open));
                rows_jj_closed = find(timetable.train_id == trains(jj)& ismember(timetable.blocksection, blocksections.id(blocksections.closed==1)));
                t_closed_jj = sum(timetable.running(rows_jj_closed));
                first_jj = rows_jj(1);
                t_approach_jj = timetable.arrival(first_jj)-timetable.start(first_jj);
                first_ii = rows_ii(1);
                t_approach_ii = timetable.arrival(first_ii)-timetable.start(first_ii);
                minHW(ii,jj) = t_all_ii + t_open_jj + default_from_sander+t_approach_jj;
                minHW(jj,ii) = t_closed_jj + default_from_sander+t_approach_ii;
            end
            if direction(ii)<10
                rows_ii_open = find(timetable.train_id == trains(ii)& ismember(timetable.blocksection, blocksections.id(blocksections.closed==0)));
                t_open_ii = sum(timetable.running(rows_ii_open));
                rows_ii_closed = find(timetable.train_id == trains(ii)& ismember(timetable.blocksection, blocksections.id(blocksections.closed==1)));
                t_closed_ii = sum(timetable.running(rows_ii_closed));
                rows_jj = find(timetable.train_id == trains(jj));
                t_all_jj = sum(timetable.running(rows_jj));
                first_ii = rows_ii(1);
                t_approach_ii = timetable.arrival(first_ii)-timetable.start(first_ii);
                first_jj = rows_jj(1);
                t_approach_jj = timetable.arrival(first_jj)-timetable.start(first_jj);                
                minHW(ii,jj) = t_closed_ii + default_from_sander+t_approach_jj;
                minHW(jj,ii) = t_all_jj + t_open_ii + default_from_sander+t_approach_ii;
            end
        end
        if abs(direction(ii)-direction(jj))==9 || abs(direction(ii)-direction(jj))==11
            rows_ii = find(timetable.train_id == trains(ii)& ismember(timetable.blocksection, blocksections.id(blocksections.closed==1)));
            t_closed_ii = sum(timetable.running(rows_ii));
            rows_jj = find(timetable.train_id == trains(jj)& ismember(timetable.blocksection, blocksections.id(blocksections.closed==1)));
            t_closed_jj = sum(timetable.running(rows_jj));
            first_jj = rows_jj(1);
            t_approach_jj = timetable.arrival(first_jj)-timetable.start(first_jj);
            first_ii = rows_ii(1);
            t_approach_ii = timetable.arrival(first_ii)-timetable.start(first_ii);            
            minHW(ii,jj) = t_closed_ii + default_from_sander + t_approach_jj;
            minHW(jj,ii) = t_closed_jj + default_from_sander + t_approach_ii;
        end
    end
end
