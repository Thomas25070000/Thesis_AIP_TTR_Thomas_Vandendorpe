function [timetable, solu] = modelCase1_singleMachine_v5(timetable, blocksections, trains, minHW, settings)
%% This version also accounts for the increase in blocking times next to the closed section.


% Path CPLEX
pathCPLEX = 'C:\Program Files\IBM\ILOG\CPLEX_Studio126\cplex\matlab\x64_win64';
addpath(genpath(pathCPLEX));
% Path YALMIP files
pathYALMIP = 'D:\OneDrive\Documenten\Thesis\Models\YALMIP';
addpath(genpath(pathYALMIP));

% Create different a machine for each block section. Keep the closed one as
% together.
blocktype = [];
bref = [];
bb = 0;
while ~blocksections.closed(bb+1)
    bb = bb + 1;
    blocktype = [blocktype blocksections.type(bb)];
    bref = [bref blocksections.id(bb)];
end
blocktype = [blocktype 'D'];
bref = [bref 0];
closedone = length(bref);
bb = find(blocksections.closed);
bb = bb(end);
while bb+1 <= size(blocksections,1)
    bb = bb + 1;
    blocktype = [blocktype blocksections.type(bb)];
    bref = [bref blocksections.id(bb)];
end

Nmachines = length(blocktype);

%% Process until we have a machine scheduling problem
Ntrains = length(trains);
direction = zeros(Ntrains,1);
% Process times = total running time over blocks.
processtime = zeros(Ntrains,Nmachines);
% Arrival events on the machines
arr_ev = zeros(Ntrains,Nmachines);
% Release times = arrival at the start of the corridor.
releasetime = zeros(Ntrains,1);
% Setup times are given by the minimum headways
setuptime = minHW;      % ONLY FOR THE CLOSED PART!
% The original departure from the corridor
deptime = zeros(Ntrains,1);
% What is the clearing time? This depends on the train type.
cleartime = zeros(Ntrains,1);

closed_blocks = find(blocksections.closed);
for tt = 1:Ntrains
    id = trains(tt);
    allev = timetable(find(timetable.train_id == id),:);
    % Events on the closed section
    ev = allev(find(ismember(allev.blocksection,closed_blocks)),:);
    % Get the needed values
    processtime(tt,closedone) = sum(ev.running);
    % Process time for the other blocks
    for bb = 1:Nmachines
        block = bref(bb);
        if bref(bb) > 0
            processtime(tt,bb) = allev.running(find(allev.blocksection == block));
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
end



%% Create the variables
Ntrains = length(trains);
% Starting time on a machine
t = sdpvar(Ntrains,Nmachines);
% Completion time on the corridor
C = sdpvar(Ntrains,1);
% Cancellation of full train
x = binvar(Ntrains,1);
% Order variables on the machines
q = binvar(Ntrains,Ntrains);
% Delay upon completion
d = sdpvar(Ntrains,1);
% Indicate on which machines the train has been stopped.
stop = binvar(Ntrains,Nmachines);
% Deviation from original release ('entrance') time
e = sdpvar(Ntrains,1);

%% Constraints
bigM = settings.disruption.duration * 3600 * 20;     % Just make it very big!
Constraints = [];

% Start time >= release time on the first section of the train
for tt = 1:Ntrains
    if direction(tt)
        % First block is the left one
        fb = 1;
    else
        % First block is at the end.
        fb = Nmachines;
    end
    minStartTime = [t(tt,fb) - releasetime(tt) * (1-x(tt)) >= 0] : ['start_' int2str(trains(tt))];
    Constraints = [Constraints, minStartTime];
end

% Time to start on the next blocksection > time on previous one + process
% time on previous one.
for tt = 1:Ntrains   
    if direction(tt)
        % Running from left to right.
        for bb = 2:Nmachines
            minStartTime = [t(tt,bb) - t(tt,bb-1) - processtime(tt,bb-1) * (1-x(tt)) >= 0]: ['intermediate_' int2str(trains(tt)) '_' int2str(bb)];
            % Note that if the time on the next block is larger, this means
            % it has been stopped before the signal!
            maxStartTime = [t(tt,bb) - t(tt,bb-1) - processtime(tt,bb-1) * (1-x(tt)) - bigM * stop(tt,bb-1) <= 0]: ['intermediate_stop_' int2str(trains(tt)) '_' int2str(bb)];
            
            % Do not allow stop on closed section!
            if bref(bb-1) == 0
                noStop = [stop(tt,bb-1) == 0]: ['no_stop_' int2str(trains(tt))];
                Constraints = [Constraints, noStop];
            end
            
            
            Constraints = [Constraints, minStartTime, maxStartTime];
        end
    else
        % From right to left, start on one-but-last.
        for bb = (Nmachines-1):-1:1
            minStartTime = [t(tt,bb) - t(tt,bb+1) - processtime(tt,bb+1) * (1-x(tt)) >= 0]: ['intermediate_' int2str(trains(tt)) '_' int2str(bb)];
            % Note that if the time on the next block is larger, this means
            % it has been stopped before the signal!
            maxStartTime = [t(tt,bb) - t(tt,bb+1) - processtime(tt,bb+1) * (1-x(tt)) - bigM * stop(tt,bb+1) <= 0]: ['intermediate_stop_' int2str(trains(tt)) '_' int2str(bb)];
          
            % Do not allow stop on closed section!
            if bref(bb+1) == 0
                noStop = [stop(tt,bb+1) == 0]: ['no_stop_' int2str(trains(tt))];
                Constraints = [Constraints, noStop];
            end
            
            
            Constraints = [Constraints, minStartTime, maxStartTime];
        end
    end
end

% Completion time > start time on one but last + process time
% What is the delay?
for tt = 1:Ntrains
    if direction(tt)
        % Last event takes place on the last machine.
        minCompletionTime = [C(tt) - t(tt,end) - processtime(tt,end) * (1-x(tt)) >= 0]: ['completion_' int2str(trains(tt))];
    else
        % Last event takes place on the first machine.
        minCompletionTime = [C(tt) - t(tt,1) - processtime(tt,1) * (1-x(tt)) >= 0]: ['completion_' int2str(trains(tt))];
    end
    delay = [d(tt) + deptime(tt) * (1-x(tt)) - C(tt) >= 0]: ['delay_' int2str(trains(tt))];
    if direction(tt)
        entrancedelay = [e(tt) + releasetime(tt) * (1-x(tt)) - t(tt,1) >= 0]: ['entrancedelay_' int2str(trains(tt))];
    else
        entrancedelay = [e(tt) + releasetime(tt) * (1-x(tt)) - t(tt,end) >= 0]: ['entrancedelay_' int2str(trains(tt))];
    end
%     Constraints = [Constraints, minCompletionTime, delay];
    Constraints = [Constraints, minCompletionTime, delay, entrancedelay];
end


% Sequencing variables: which comes after which?
% Note that some will have to be fixed to 0 or 1 for the same direction!
seqConstraints = [];
closedM = find(bref==0);
for ii = 1:Ntrains
    for jj = ii:Ntrains
        if jj>=ii
            % Decision on the order will be made on the closed part.
            label = ['seq_' int2str(ii) '_' int2str(jj) '_disrupted'];
            seq1 = [t(jj,closedM) + bigM * (1 - q(ii,jj)) - t(ii,closedM) ...
                    - setuptime(ii,jj) * (1 - x(ii)) >= 0]: label;
            seq2 = [t(ii,closedM) + bigM * q(ii,jj) - t(jj,closedM) ...
                    - setuptime(jj,ii) * (1 - x(jj)) >= 0]: label;
            seqConstraints = [seqConstraints, seq1, seq2];
       
            if direction(jj) == direction(ii)
                % We have to fix the q-value!
                label = ['fix q(' int2str(ii) ',' int2str(jj) ')'];
                if releasetime(ii) < releasetime(jj)
                    Q = [q(ii,jj) >= 1 - x(ii) - x(jj)] : label;
                else
                    Q = [q(ii,jj) <= x(ii) + x(jj)] : label;
                end
                seqConstraints = [seqConstraints, Q];

            end
        end
    end
end
Constraints = [Constraints, seqConstraints];

t_before = settings.TT.blocktimes.setup;
t_after_R = settings.TT.blocktimes.afterR;
t_after_IC = settings.TT.blocktimes.afterIC;
% We still need the minimum separation on the other machines of the
% corridor!
corridorSeq = [];
for ii = 1:Ntrains
    for jj = ii:Ntrains
        if jj>ii
            if direction(ii) == direction(jj)
                % Only if they are in the same direction!
                if direction(ii) == 1
                    % Both trains run from left to right!
                    if releasetime(jj) > releasetime(ii)
                        
                        for mm = 1:Nmachines-1
                            if bref(mm) > 0
                                try
                                    t_setup = cleartime(ii) + t_before + processtime(jj,mm-1);
                                    % Gives wrong result if signals are
                                    % present!
                                    if bref(mm-1) == 0 && (direction(ii) ~= settings.disruption.direction || settings.disruption.signals)
                                        % Previous blocksection is the
                                        % closed one, but blocking time is
                                        % acceptable.
                                        index = find(timetable.train_id == trains(jj) & timetable.blocksection == bref(mm-1));
                                        t_setup = cleartime(ii) + t_before + timetable.running(index);
                                    end
                                catch
                                    t_setup = cleartime(ii) + t_before + processtime(jj,1);
                                end
                               
                                seqConstraint = [t(jj,mm) - t(ii,mm+1) - t_setup*(1-x(ii)) >= 0];
                                corridorSeq = [corridorSeq seqConstraint];
                            end
                        end
                        % Also for the completing one!
                        seqConstraint = [t(jj,end) - C(ii) - t_setup*(1-x(ii)) >= 0];
                        corridorSeq = [corridorSeq seqConstraint];
                    else
  
                        for mm = 1:Nmachines-1
                            if bref(mm) > 0         % Don't do it for the disrupted section!
                                try
                                    t_setup = cleartime(jj) + t_before + processtime(ii,mm-1);
                                    % Gives wrong result if signals are
                                    % present!
                                    if bref(mm-1) == 0 && (direction(jj) ~= settings.disruption.direction || settings.disruption.signals)
                                        % Previous blocksection is the
                                        % closed one, but blocking time is
                                        % acceptable.
                                        index = find(timetable.train_id == trains(ii) & timetable.blocksection == bref(mm-1));
                                        t_setup = cleartime(ii) + t_before + timetable.running(index);
                                    end
                                catch
                                    t_setup = cleartime(jj) + t_before + processtime(ii,1);
                                end
                                
                                seqConstraint = [t(ii,mm) - t(jj,mm+1) - t_setup*(1-x(jj)) >= 0];
                                corridorSeq = [corridorSeq seqConstraint];
                            end
                        end
                        % Also for the completing one!
                        seqConstraint = [t(ii,end) - C(jj) - t_setup*(1-x(jj)) >= 0];
                        corridorSeq = [corridorSeq seqConstraint];
                    end
                else
                    % Trains go in the other direction!
                    if releasetime(jj) > releasetime(ii)
%                         t_setup = cleartime(ii) + t_before + processtime(jj,end);
                        for mm = Nmachines:-1:2
                            if bref(mm)>0
                                try
                                    t_setup = cleartime(ii) + t_before + processtime(jj,mm+1);
                                    % Gives wrong result if signals are
                                    % present!
                                    if bref(mm+1) == 0 && (direction(ii) ~= settings.disruption.direction || settings.disruption.signals)
                                        % Previous blocksection is the
                                        % closed one, but blocking time is
                                        % acceptable.
                                        index = find(timetable.train_id == trains(jj) & timetable.blocksection == bref(mm+1));
                                        t_setup = cleartime(ii) + t_before + timetable.running(index);
                                    end
                                catch
                                    t_setup = cleartime(ii) + t_before + processtime(jj,end);
                                end
                                
                                seqConstraint = [t(jj,mm) - t(ii,mm-1) - t_setup*(1-x(ii)) >= 0];
                                corridorSeq = [corridorSeq seqConstraint];
                            end
                        end
                        % Also for the completing one!
                        seqConstraint = [t(jj,1) - C(ii) - t_setup*(1-x(ii)) >= 0];
                        corridorSeq = [corridorSeq seqConstraint];
                    else
%                         t_setup = cleartime(jj) + t_before + processtime(ii,end);
                        for mm = Nmachines:-1:2
                            if bref(mm)>0
                                try
                                    t_setup = cleartime(jj) + t_before + processtime(ii,mm+1);
                                    % Gives wrong result if signals are
                                    % present!
                                    if bref(mm+1) == 0 && (direction(ii) ~= settings.disruption.direction || settings.disruption.signals)
                                        % Previous blocksection is the
                                        % closed one, but blocking time is
                                        % acceptable.
                                        index = find(timetable.train_id == trains(ii) & timetable.blocksection == bref(mm+1));
                                        t_setup = cleartime(jj) + t_before + timetable.running(index);
                                    end
                                catch
                                    t_setup = cleartime(jj) + t_before + processtime(ii,end);
                                end
                                
                                seqConstraint = [t(ii,mm) - t(jj,mm-1) - t_setup*(1-x(jj)) >= 0];
                                corridorSeq = [corridorSeq seqConstraint];
                            end
                        end
                        % Also for the completing one!
                        seqConstraint = [t(ii,1) - C(jj) - t_setup*(1-x(jj)) >= 0];
                        corridorSeq = [corridorSeq seqConstraint];
                    end
                end
            end
        end
    end
end
Constraints = [Constraints, corridorSeq];



% Off-balance constraint
% Which trains run in direction 1?
if settings.disruption.aggregateOB
    dir1 = find(mod(trains,10)==1);
    dir0 = setdiff([1:Ntrains],dir1);
    delta = settings.disruption.duration * settings.disruption.offbalance;
    LBbalance = [sum(x(dir1)) - sum(x(dir0)) >= -delta]: 'LB_balance';
    UBbalance = [sum(x(dir1)) - sum(x(dir0)) <= delta]: 'UB_balance';
    Constraints = [Constraints, LBbalance, UBbalance];
else
    % If we want this per hour
    delta = settings.disruption.offbalance;
    for hh = 1:settings.disruption.duration
        dir1 = find((mod(trains,10)==1) & (floor(trains/1000)==(hh-1)));
        dir0 = find((mod(trains,10)==0) & (floor(trains/1000)==(hh-1)));
        LBbalance = [sum(x(dir1)) - sum(x(dir0)) >= -delta]: ['LB_balance_hour_' int2str(hh)];
        UBbalance = [sum(x(dir1)) - sum(x(dir0)) <= delta]: ['UB_balance_hour_' int2str(hh)];
        Constraints = [Constraints, LBbalance, UBbalance];
    end
end


%% Build the objective function
Objective = 0;

switch settings.TT.paxdirection
    case -1
        w_dir0 = 1;
        w_dir1 = 1;
    case 0
        w_dir0 = 2;
        w_dir1 = 1;
    case 1
        w_dir0 = 1;
        w_dir1 = 2;
end

% Penalize the delays
for tt = 1:Ntrains
%     Objective = Objective + C(tt) - deptime(tt);
    if mod(trains(tt),10) == 1
        Objective = Objective + w_dir1 * d(tt);
        Objective = Objective + w_dir1 * e(tt)/1000;     % Very small penalty for an entrance delay.
    else
        Objective = Objective + w_dir0 * d(tt);
        Objective = Objective + w_dir0 * e(tt)/1000;
    end
end

% Penalize the cancellation of trains
w_cancel = 2400;
for tt = 1:Ntrains
    if mod(trains(tt),10) == 1
        Objective = Objective + w_cancel * w_dir1 * x(tt);
    else
        Objective = Objective + w_cancel * w_dir0 * x(tt);
    end
end

%% Build and solve the model
ops = sdpsettings('solver','cplex');
ops.savesolveroutput = 1;

[solu] = optimize(Constraints,Objective,ops);
status = solu.solveroutput.output.cplexstatus;
disp(int2str(status));

if status == 103
    disp('Infeasible!')
    pause();
else
    arrivals =  round(value(t));
    departures = round(value(C));
    cancelled = round(value(x));
    
    cancelled_ev = [];
    for xx = 1:length(cancelled)
        if cancelled(xx)
            rows = find(timetable.train_id == trains(xx));
            cancelled_ev = [cancelled_ev, rows'];
        end
    end
    timetable.cancelled(1) = 0;
    timetable.cancelled(cancelled_ev) = 1;
    
    timetable.adjusted_arrival(1) = 0;
    timetable.adjusted_departure(1) = 0;
    for ii = 1:size(arr_ev,1)
        for jj = 1:size(arr_ev,2)
            timetable.adjusted_arrival(arr_ev(ii,jj)) = arrivals(ii,jj);
            timetable.adjusted_departure(arr_ev(ii,jj)) = arrivals(ii,jj)...
                            + timetable.running(arr_ev(ii,jj));
        end
    end
    timetable.adjusted_departure(completion_ev) = departures';
    
    % Update all timings
    new_timetable = updateTimetable_v3(timetable, blocksections, settings, trains, cancelled);
    
    
    % Plot new timetable
%     type = 'hour';
    type = 'complete';
%     type = 'base';
    include_blocks = 1;
    directions = 1;
    directions = 0;
    [line_new, blocks_new] = plotTT_new(new_timetable, blocksections, settings, type, include_blocks, directions);
    
end

end
