function [minHWstack, minHWstack_reg, otherinfo] = createHeadwayStackDisruption(timetable, regular, blocksections, blocktype, bref, settings)

regTT = regular.TT;

trains = unique(timetable.train_id);
Ntrains = length(trains);
Nmachines = length(blocktype);
direction = zeros(Ntrains,1);

% Process times = total running time over blocks.
processtime = zeros(Ntrains,Nmachines);
regprocesstime = zeros(Ntrains,Nmachines);
% Approach time
approachtime = zeros(Ntrains,Nmachines);
regapproachtime = zeros(Ntrains,Nmachines);

% Arrival events on the machines
arr_ev = zeros(Ntrains,Nmachines);
% Release times = arrival at the start of the corridor.
releasetime = zeros(Ntrains,1);

%% Determine all the timings!
closed_blocks = find(blocksections.closed);
closedone = find(bref == 0);
for tt = 1:Ntrains
    id = trains(tt);
    allev = timetable(find(timetable.train_id == id),:);
    % Events on the closed section
    ev = allev(find(ismember(allev.blocksection,closed_blocks)),:);
    % Get the needed values
    processtime(tt,closedone) = sum(ev.running);
    % Regular values, i.e. not during disruption
    all_reg_ev = regTT(find(regTT.train_id == id),:);
    reg_ev = all_reg_ev(find(ismember(all_reg_ev.blocksection,closed_blocks)),:);
    regprocesstime(tt,closedone) = sum(reg_ev.running);
    % Process time for the other blocks
    for bb = 1:Nmachines
        block = bref(bb);
        if bref(bb) > 0
            processtime(tt,bb) = allev.running(find(allev.blocksection == block));
            regprocesstime(tt,bb) = all_reg_ev.running(find(all_reg_ev.blocksection == block));
            arr_ev(tt,bb) = allev.event_id(find(allev.blocksection == block));
        else
            arr_ev(tt,bb) = ev.event_id(1);
        end
    end
    
    release_ev(tt) = allev.event_id(1);        % First event of this train on the closed section.
    releasetime(tt) = allev.arrival(1);
    direction(tt) = allev.direction(1);
%     deptime(tt) = ev.departure(end);
    deptime(tt) = allev.arrival(1) + sum(allev.running);
    completion_ev(tt) = allev.event_id(end);   % Last event of this train on the closed section.

    switch ev.train_type{1}
        case 'IC'
            cleartime(tt) = settings.TT.blocktimes.afterIC;
        case 'R'
            cleartime(tt) = settings.TT.blocktimes.afterR;
    end
    
    if ev.arrival(1) >= settings.disruption.duration*3600
        outside_disruption_period(tt) = 1;
    end
end

otherinfo.arr_ev = arr_ev;
otherinfo.release_ev = release_ev;
otherinfo.direction = direction;
otherinfo.deptime = deptime;
otherinfo.completion_ev = completion_ev;
otherinfo.outside_disruption_period = outside_disruption_period;
otherinfo.releasetime = releasetime;

%% Fill the approach times
for tt = 1:Ntrains
    % Depending on direction!
    if direction(tt) == 1
        for mm = 2:Nmachines
            approachtime(tt,mm) = processtime(tt,mm-1);
            regapproachtime(tt,mm) = regprocesstime(tt,mm-1);
            % BUT: not if running on the closed section + signals are present!
            % flag indicates that we have to correct
            flag = ~(((direction(tt) == settings.disruption.direction) && ~settings.disruption.signals)...
                            || outside_disruption_period(tt));
            if flag && (bref(mm-1)==0)
                closedev = timetable(find(timetable.train_id == trains(tt) ...
                                        & blocksections.closed(timetable.blocksection)),:);
                approachtime(tt,mm) = closedev.running(end);
                closedev = regTT(find(regTT.train_id == trains(tt) ...
                                        & blocksections.closed(regTT.blocksection)),:);
                regapproachtime(tt,mm) = closedev.running(end);
            end
            if flag && (bref(mm)==0)
                closedev = timetable(find(timetable.train_id == trains(tt) ...
                                        & blocksections.closed(timetable.blocksection)),:);
                approachtime(tt,mm) = max(closedev.running);
            end
        end
    else
        for mm = (Nmachines-1):-1:1
            approachtime(tt,mm) = processtime(tt,mm+1);
            regapproachtime(tt,mm) = regprocesstime(tt,mm+1);
            % BUT: not if running on the closed section + signals are present!
            % flag indicates that we have to correct
            flag = ~((direction(tt) == settings.disruption.direction) && ~settings.disruption.signals) ...
                            || outside_disruption_period(tt);
            if flag && (bref(mm+1)==0)
                closedev = timetable(find(timetable.train_id == trains(tt) ...
                                        & blocksections.closed(timetable.blocksection)),:);
                approachtime(tt,mm) = closedev.running(end);
                closedev = regTT(find(regTT.train_id == trains(tt) ...
                                        & blocksections.closed(regTT.blocksection)),:);
                regapproachtime(tt,mm) = closedev.running(end);
            end
            if flag && (bref(mm)==0)
                closedev = timetable(find(timetable.train_id == trains(tt) ...
                                        & blocksections.closed(timetable.blocksection)),:);
                approachtime(tt,mm) = max(closedev.running);
            end
        end
    end
end

%% What headway do we need between two trains on each considered machine?
minHWstack = zeros(Ntrains,Ntrains,Nmachines);
minHWstack_reg = zeros(Ntrains,Ntrains,Nmachines);

% Start with the closed part.
% Headways needed for trains in both directions, depending on the presence
% of signals.
mm = find(bref == 0);

t_setup = settings.TT.blocktimes.setup;
for ii = 1:Ntrains
    for jj = 1:Ntrains
        if ii ~= jj
            % After event i and before event j:
            % running_i + clearing_i + setup + approach_j
            hw_ii_jj = processtime(ii,mm) + cleartime(ii) + t_setup ...
                            + approachtime(jj,mm);
            reg_hw_ii_jj = regprocesstime(ii,mm) + cleartime(ii) + t_setup ...
                            + regapproachtime(jj,mm);
            minHWstack(ii,jj,mm) = hw_ii_jj;
            minHWstack_reg(ii,jj,mm) = reg_hw_ii_jj;
        end
    end
end

% For the other machines, we only have to do it if they happen in the same
% direction.
for mm = 1:Nmachines
    if mm ~= find(bref == 0)
        for ii = 1:Ntrains
            for jj = 1:Ntrains
                if ii ~= jj && direction(ii) == direction(jj)
                    % After event i and before event j:
                    % running_i + clearing_i + setup + approach_j
                    hw_ii_jj = processtime(ii,mm) + cleartime(ii) + t_setup ...
                                    + approachtime(jj,mm);
                    reg_hw_ii_jj = regprocesstime(ii,mm) + cleartime(ii) + t_setup ...
                                    + regapproachtime(jj,mm);
                    minHWstack(ii,jj,mm) = hw_ii_jj;
                    minHWstack_reg(ii,jj,mm) = reg_hw_ii_jj;
                end
            end
        end
    end
end

otherinfo.processtime = processtime;
otherinfo.regprocesstime = regprocesstime;
otherinfo.approachtime = approachtime;
otherinfo.regapproachtime = regapproachtime;
otherinfo.cleartime = cleartime;



end