function [timetable, solu, measures, statistics] = modelCase1_singleMachine_v11(timetable, regular, blocksections, trains, minHW, settings, runningtimes)
%% This version allows to change the orders of the incoming trains during 
% the disruption.

try
    warmstart = settings.useWarmstart;
catch
    warmstart = 0;
end

try
    maxD = settings.disruption.maxDelay;
catch
    maxD = 3600;
end
if isempty(maxD)
    maxD = 3600;
end
try
    w_cancel = settings.weights.cancel;
catch
    w_cancel = 3600;
end
if isempty(w_cancel)
    w_cancel = 3600;
end
try
    w_late_entry = settings.weights.entranceDelay;
catch
    w_late_entry = 1/1000;
end
if isempty(w_late_entry)
    w_late_entry = 1/1000;
end
try
    w_delay = settings.weights.exitDelay;
catch
    w_delay = 1;
end
if isempty(w_delay)
    w_delay = 1;
end
try
    WIC = settings.disruption.WICallowed;
catch
    WIC = 0;
end
if ~WIC
    w_late_entry = 0;
end

clearingtimes = settings.clearingtimes;


%% This version also accounts for the increase in blocking times next to the closed section.
max_E = 5000;
% Penalty for cancellation
%w_cancel = 2400;
%w_late_entry = 1/1000;

regTT = regular.TT;
regHW = regular.HW;

% Path CPLEX
pathCPLEX = 'C:\Program Files\IBM\ILOG\CPLEX_Studio126\cplex\matlab\x64_win64';
addpath(genpath(pathCPLEX));
addpath(genpath('C:\gurobi651\win64'));
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

% Latest point we can change, is right before the corridor! So for trains
% coming from direction 1, it is the first block of D. For trains from
% direction 0, the last one
S1 = find(strcmp(blocktype,'D'));
S1 = S1(1);
S2 = find(strcmp(blocktype,'D'));
S2 = S2(end);


%% Process until we have a machine scheduling problem
Ntrains = length(trains);
direction = zeros(Ntrains,1);
% Process times = total running time over blocks.
processtime = zeros(Ntrains,Nmachines);
regprocesstime = zeros(Ntrains,Nmachines);
% Arrival events on the machines
arr_ev = zeros(Ntrains,Nmachines);
% Release times = arrival at the start of the corridor.
releasetime = zeros(Ntrains,1);
% Setup times are given by the minimum headways
setuptime = minHW;      % ONLY FOR THE CLOSED PART!
% The original departure from the corridor
deptime = zeros(Ntrains,1);
orig_deptime = zeros(Ntrains,1);
% What is the clearing time? This depends on the train type.
cleartime = zeros(Ntrains,Nmachines);
% Approach time
approachtime = zeros(Ntrains,Nmachines);    % If both in the same direction.
approach_other = zeros(Ntrains,1);          % If in the other direction, i.e. on the closed section.
regapproachtime = zeros(Ntrains,Nmachines);
% Outside disruption period
outside_disruption_period = zeros(Ntrains,1);
% End of disruption
t_end = settings.disruption.duration * 3600;


closed_blocks = find(blocksections.closed);
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
			switch ev.train_type{1}
                case 'IC'
                    if direction
                        cleartime(tt,bb) = clearingtimes.IC1.disrupted(bref(bb));
                    else
                        cleartime(tt,bb) = clearingtimes.IC0.disrupted(bref(bb));
                    end
                case 'R'
                    if direction
                        cleartime(tt,bb) = clearingtimes.L1.disrupted(bref(bb));
                    else
                        cleartime(tt,bb) = clearingtimes.L0.disrupted(bref(bb));
                    end
            end
        else
            arr_ev(tt,bb) = ev.event_id(1);
			switch ev.train_type{1}
                case 'IC'
                    if direction
                        cleartime(tt,bb) = max(clearingtimes.IC1.disrupted);
                    else
                        cleartime(tt,bb) = max(clearingtimes.IC0.disrupted);
                    end
                case 'R'
                    if direction
                        cleartime(tt,bb) = max(clearingtimes.L1.disrupted);
                    else
                        cleartime(tt,bb) = max(clearingtimes.L0.disrupted);
                    end
            end
        end
    end
    
    release_ev(tt) = allev.event_id(1);        % First event of this train on the closed section.
    releasetime(tt) = allev.arrival(1);
    direction(tt) = allev.direction(1);
    traintype(tt) = allev.train_type(1);
%     deptime(tt) = ev.departure(end);
    deptime(tt) = allev.arrival(1) + sum(allev.running);
    orig_deptime(tt) = allev.departure(end);
    completion_ev(tt) = allev.event_id(end);   % Last event of this train on the closed section.

%    switch ev.train_type{1}
 %       case 'IC'
  %          cleartime(tt) = settings.TT.blocktimes.afterIC;
  %      case 'R'
  %          cleartime(tt) = settings.TT.blocktimes.afterR;
  %  end
    
    if ev.arrival(1) >= settings.disruption.duration*3600
        outside_disruption_period(tt) = 1;
    end
end

% Fill the approach times
for tt = 1:Ntrains
    % Depending on direction!
    if direction(tt) == 1
        for mm = 1:Nmachines
            if mm == 1
            	[approachtime(tt,1),regapproachtime(tt,1)] = returnFirstApproachTime(blocksections,settings, traintype(tt), direction(tt));
            else
                approachtime(tt,mm) = processtime(tt,mm-1);
				regapproachtime(tt,mm) = regprocesstime(tt,mm-1);
            end
            % BUT: not if running on the closed section + signals are present!
            % flag indicates that we have to correct
            flag = ~(((direction(tt) == settings.disruption.direction) && ~settings.disruption.signals)...
                            || outside_disruption_period(tt));
            if mm > 1 && flag && (bref(mm-1)==0)
                closedev = timetable(find(timetable.train_id == trains(tt) ...
                                        & blocksections.closed(timetable.blocksection)),:);
                approachtime(tt,mm) = closedev.running(end);

                closedev = regTT(find(regTT.train_id == trains(tt) ...
                                        & blocksections.closed(regTT.blocksection)),:);
                regapproachtime(tt,mm) = max(closedev.running);
            elseif mm > 1 && (bref(mm-1)==0)
                closedev = regTT(find(regTT.train_id == trains(tt) ...
                                        & blocksections.closed(regTT.blocksection)),:);
                regapproachtime(tt,mm) = max(closedev.running);
            end 
		end


  %          if flag && (bref(mm)==0)
%                 closedev = timetable(find(timetable.train_id == trains(tt) ...
%                                         & blocksections.closed(timetable.blocksection)),:);
%                 approachtime(tt,mm) = max(closedev.running);
		flag = ~((direction(tt) == settings.disruption.direction) && ~settings.disruption.signals);
        if flag
            mm = find(bref == 0);
            closedev = find(timetable.train_id == trains(tt) ...
                                        & blocksections.closed(timetable.blocksection));
            if mm == 1
                % It is the first machine, watch out!
                summed = timetable.running(closedev(2:end)) + timetable.running(closedev(2:end)-1);
                maxsummed = closedev(find(summed == max(summed)));
                firstmach = timetable.running(closedev(1)) + approachtime(tt,1);
                if firstmach < max(summed)
                    approachtime(tt,mm) = max(timetable.running(maxsummed(1) - 1), ...
                                                approachtime(tt,mm));
                end
                approach_other(tt) = approachtime(tt,1);
            else
                summed = timetable.running(closedev) + timetable.running(closedev-1);

                maxsummed = closedev(find(summed == max(summed)));
%                 approachtime(tt,mm) = max(timetable.running(maxsummed(1) - 1), ...
%                                             approachtime(tt,mm));
                approachtime(tt,mm) = timetable.running(maxsummed(1) - 1);

                approach_other(tt) = timetable.running(closedev(1)-1);

            end
        end
    else
        for mm = Nmachines:-1:1
            if mm == Nmachines
                [approachtime(tt,end),regapproachtime(tt,end)] = returnFirstApproachTime(blocksections,settings, traintype(tt), direction(tt));
            else
                 approachtime(tt,mm) = processtime(tt,mm+1);
				 regapproachtime(tt,mm) = regprocesstime(tt,mm+1);
            end
            % BUT: not if running on the closed section + signals are present!
            % flag indicates that we have to correct
            flag = ~((direction(tt) == settings.disruption.direction) && ~settings.disruption.signals) ...
                            || outside_disruption_period(tt);
            if mm < Nmachines && flag && (bref(mm+1)==0)
                closedev = timetable(find(timetable.train_id == trains(tt) ...
                                        & blocksections.closed(timetable.blocksection)),:);
                approachtime(tt,mm) = closedev.running(end);
                closedev = regTT(find(regTT.train_id == trains(tt) ...
                                        & blocksections.closed(regTT.blocksection)),:);
				regapproachtime(tt,mm) = max(closedev.running);

            
            elseif mm < Nmachines && (bref(mm+1)==0)
                closedev = regTT(find(regTT.train_id == trains(tt) ...
                                        & blocksections.closed(regTT.blocksection)),:);
                regapproachtime(tt,mm) = max(closedev.running);
            end 
		end
		flag = ~((direction(tt) == settings.disruption.direction) && ~settings.disruption.signals);
        if flag
            mm = find(bref == 0);

            closedev = find(timetable.train_id == trains(tt) ...
                                        & blocksections.closed(timetable.blocksection));
            if mm == Nmachines
                % It is the first machine, watch out!
                summed = timetable.running(closedev(2:end)) + timetable.running(closedev(2:end)-1);
                maxsummed = closedev(find(summed == max(summed)));
                lastmach = timetable.running(closedev(1)) + approachtime(tt,end);
                if lastmach < max(summed)
                    approachtime(tt,mm) = max(timetable.running(maxsummed(1) - 1),...
                                            approachtime(tt,mm));
                end
                approach_other(tt) = approachtime(tt,end);
            else



                summed = timetable.running(closedev) + timetable.running(closedev-1);


                maxsummed = closedev(find(summed == max(summed)));
%                 approachtime(tt,mm) = max(timetable.running(maxsummed(1) - 1),...
%                                             approachtime(tt,mm));
                approachtime(tt,mm) = timetable.running(maxsummed(1) - 1);





                approach_other(tt) = timetable.running(closedev(1)-1);

            end
        end
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
% Binaries to indicate whether the first arrival on the corridor happens
% after the end of the disruption or not.
gamma = binvar(Ntrains,1);

if warmstart
    assign(x,0);
end

%% Constraints
bigM = settings.disruption.duration * 3600 * 2;     % Just make it very big!
Constraints = [];

% Start time >= release time on the first section of the train + indicator
% whether it is after the disruption ending.
for tt = 1:Ntrains
    if direction(tt)
        % First block is the left one
        fb = 1;
        sw = S1;
    else
        % First block is at the end.
        fb = Nmachines;
        sw = S2;
    end
    minStartTime = [t(tt,fb) - releasetime(tt) * (1-x(tt)) >= 0] : ['start_' int2str(trains(tt))];
    % Refer to the disruption end at the time of entry on the switches!
    % This is the last moment on which the route can be changed.
    disruptionend1 = [t(tt,sw) - t_end * gamma(tt) >= 0]: ['end_of_disruption_' int2str(trains(tt))];
    disruptionend2 = [t(tt,sw) - bigM * gamma(tt) <= t_end]: ['end_of_disruption_' int2str(trains(tt))];
    Constraints = [Constraints, minStartTime, disruptionend1, disruptionend2];
    
    % Fix gamma?
    if releasetime(tt) + w_cancel < t_end
        % Never able to cross the disruption border, better to cancel it
        % than delay it that much.
         fixGamma = [gamma(tt) == 0]: ['fix_gamma_' int2str(trains(tt))];
         Constraints = [Constraints, fixGamma];
     elseif releasetime(tt) >= t_end
         fixGamma = [gamma(tt) == 1 - x(tt)]: ['fix_gamma_' int2str(trains(tt))];
         Constraints = [Constraints, fixGamma];
    end
    
	% If the trains is cancelled, gamma is always 0!
    fixGamma = [gamma(tt) + x(tt) <= 1];
    Constraints = [Constraints, fixGamma];
end

% Is WIC allowed?
if ~WIC
    noStop = [stop == 0];
    Constraints = [Constraints, noStop];

end

% Time to start on the next blocksection > time on previous one + process
% time on previous one.
for tt = 1:Ntrains   
    if direction(tt)
        % Running from left to right.
        for bb = 2:Nmachines
            delta_proc = regprocesstime(tt,bb-1) - processtime(tt,bb-1);
            
%             minStartTime = [t(tt,bb) - t(tt,bb-1) - processtime(tt,bb-1) * (1-x(tt)) >= 0]: ['intermediate_' int2str(trains(tt)) '_' int2str(bb)];
            minStartTime = [t(tt,bb) - t(tt,bb-1) - processtime(tt,bb-1) ...
                * (1-x(tt)) - delta_proc * gamma(tt) >= 0]: ['intermediate_' ...
                int2str(trains(tt)) '_' int2str(bb)];
            % Note that if the time on the next block is larger, this means
            % it has been stopped before the signal!
%             maxStartTime = [t(tt,bb) - t(tt,bb-1) - processtime(tt,bb-1) * (1-x(tt)) - bigM * stop(tt,bb-1) <= 0]: ['intermediate_stop_' int2str(trains(tt)) '_' int2str(bb)];
            maxStartTime = [t(tt,bb) - t(tt,bb-1) - processtime(tt,bb-1) ...
                * (1-x(tt)) - delta_proc * gamma(tt) - bigM * stop(tt,bb-1) <= 0]: ...
                ['intermediate_stop_' int2str(trains(tt)) '_' int2str(bb)];
            
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
            delta_proc = regprocesstime(tt,bb+1) - processtime(tt,bb+1);
            
            
%             minStartTime = [t(tt,bb) - t(tt,bb+1) - processtime(tt,bb+1) * (1-x(tt)) >= 0]: ['intermediate_' int2str(trains(tt)) '_' int2str(bb)];
            minStartTime = [t(tt,bb) - t(tt,bb+1) - processtime(tt,bb+1) ...
                * (1-x(tt)) - delta_proc * gamma(tt) >= 0]: ['intermediate_' ...
                int2str(trains(tt)) '_' int2str(bb)];
            
            % Note that if the time on the next block is larger, this means
            % it has been stopped before the signal!
%             maxStartTime = [t(tt,bb) - t(tt,bb+1) - processtime(tt,bb+1) * (1-x(tt)) - bigM * stop(tt,bb+1) <= 0]: ['intermediate_stop_' int2str(trains(tt)) '_' int2str(bb)];
            maxStartTime = [t(tt,bb) - t(tt,bb+1) - processtime(tt,bb+1) ...
                * (1-x(tt)) - delta_proc * gamma(tt) - bigM * stop(tt,bb+1) <= 0]: ...
                ['intermediate_stop_' int2str(trains(tt)) '_' int2str(bb)];
            
            
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
       delta_proc = processtime(tt,end) - regprocesstime(tt,end);
        minCompletionTime = [C(tt) - t(tt,end) - processtime(tt,end) * (1-x(tt))...
             + delta_proc * gamma(tt) >= 0]: ['completion_' int2str(trains(tt))];
    else
        delta_proc = processtime(tt,1) - regprocesstime(tt,1);
        % Last event takes place on the first machine.
        minCompletionTime = [C(tt) - t(tt,1) - processtime(tt,1) * (1-x(tt)) ...
            + delta_proc * gamma(tt) >= 0]: ['completion_' int2str(trains(tt))];
    end
%     delay = [d(tt) + orig_deptime(tt) * (1-x(tt)) - C(tt) >= 0]: ['delay_' int2str(trains(tt))];
    if releasetime(tt) >= t_end
        delay = [d(tt) + orig_deptime(tt) * (1-x(tt)) - C(tt) >= 0]: ['delay_' int2str(trains(tt))];
    else
        delay = [d(tt) + deptime(tt) * (1-x(tt)) - C(tt) >= 0]: ['delay_' int2str(trains(tt))]; 
    end
        
        %     if direction(tt)
%         entrancedelay = [e(tt) + releasetime(tt) * (1-x(tt)) - t(tt,1) >= 0]: ['entrancedelay_' int2str(trains(tt))];
%     else
%         entrancedelay = [e(tt) + releasetime(tt) * (1-x(tt)) - t(tt,end) >= 0]: ['entrancedelay_' int2str(trains(tt))];
%     end
	delayPositive = [d(tt) >= 0]: ['positiveDelay_' int2str(trains(tt))];
    maxDelay = [d(tt) <= maxD];% ['maxDelay_' int2str(trains(tt))];
    if direction(tt)
        entrancedelay = [e(tt) + releasetime(tt) * (1-x(tt)) - t(tt,1) >= 0]: ['entrancedelay_' int2str(trains(tt))];
    else
        entrancedelay = [e(tt) + releasetime(tt) * (1-x(tt)) - t(tt,end) >= 0]: ['entrancedelay_' int2str(trains(tt))];
    end
    minEtt = [e(tt) >= 0];
%     maxEtt = [e(tt) <= max_E];
%     Constraints = [Constraints, minCompletionTime, delay];
    Constraints = [Constraints, minCompletionTime, delay, maxDelay, delayPositive, entrancedelay, minEtt];
end


% Sequencing variables: which comes after which?
% Note that some will have to be fixed to 0 or 1 for the same direction!
seqConstraints = [];
closedM = find(bref==0);
t_before = settings.TT.blocktimes.setup;

for ii = 1:Ntrains
    for jj = ii:Ntrains
        if jj>ii
            
            if (direction(ii) == settings.disruption.direction) && ~settings.disruption.signals
                % We have to compare with the actual exit of ii!
                flag_exit_ii = 1;
            else
                flag_exit_ii = 0;
            end
            if (direction(jj) == settings.disruption.direction) && ~settings.disruption.signals
                % We have to compare with the actual exit of jj!
                flag_exit_jj = 1;
            else
                flag_exit_jj = 0;
            end
            
            
            if direction(ii) == 1
                first_machine_after_ii = find(bref==0) + 1;
            else
                first_machine_after_ii = find(bref==0) - 1;
            end
            if direction(jj) == 1
                first_machine_after_jj = find(bref==0) + 1;
            else
                first_machine_after_jj = find(bref==0) - 1;
            end
            
            if first_machine_after_ii > Nmachines
                % In the other direction, it will be 0 anyway!
                first_machine_after_ii = 0;
            end
            if first_machine_after_jj > Nmachines
                first_machine_after_jj = 0;
            end
            
            if direction(ii) ~= direction(jj)
                t_setup_ii_jj = cleartime(ii,closedM) + approach_other(jj) + t_before;
                t_setup_jj_ii = cleartime(jj,closedM) + approach_other(ii) + t_before;
            else
                t_setup_ii_jj = cleartime(ii,closedM) + approachtime(jj,closedM) + t_before;
                t_setup_jj_ii = cleartime(jj,closedM) + approachtime(ii,closedM) + t_before;
            end
            
            label = ['seq_' int2str(ii) '_' int2str(jj) '_disrupted'];
            
            % What is the reference event?
            if (flag_exit_ii || (direction(ii) ~= direction(jj))) && first_machine_after_ii==0
                ref_event_ii = C(ii);
            elseif (flag_exit_ii || (direction(ii) ~= direction(jj)))
                ref_event_ii = t(ii,first_machine_after_ii);
            else
                % Look at the time ii exits the first block section of the
                % closed part.
                closedev = timetable(find(timetable.train_id == trains(ii) ...
                                        & blocksections.closed(timetable.blocksection)),:);
                ref_event_ii = t(ii,closedM) + max(closedev.running);
            end
            
            closedev = regTT(find(regTT.train_id == trains(ii) ...
                                        & blocksections.closed(regTT.blocksection)),:);
            ref_event_ii_end = t(ii,closedM) + closedev.running(1);
            closedev = regTT(find(regTT.train_id == trains(jj) ...
                                        & blocksections.closed(regTT.blocksection)),:);
            ref_event_jj_end = t(jj,closedM) + closedev.running(1);
            
            
            % What is the reference event?
            if (flag_exit_jj || (direction(ii) ~= direction(jj))) && first_machine_after_jj==0
                ref_event_jj = C(jj);
            elseif (flag_exit_jj || (direction(ii) ~= direction(jj)))
                ref_event_jj = t(jj,first_machine_after_jj);
            else
                % Look at the time ii exits the first block section of the
                % closed part.
                closedev = timetable(find(timetable.train_id == trains(jj) ...
                                        & blocksections.closed(timetable.blocksection)),:);
                ref_event_jj = t(jj,closedM) + max(closedev.running);
            end
            
            % If trains fall outside of the disruption period, it does not
            % matter anymore if they go in different directions!
            flag_using_same_track = 0;
            if direction(ii) == direction(jj)
                flag_using_same_track = 1;
            elseif direction(ii) ~= direction(jj) && ~outside_disruption_period(ii) ...
                        && ~outside_disruption_period(jj)
                flag_using_same_track = 1;
            end
               

            if direction(ii) ~= direction(jj)
                if direction(ii) == settings.disruption.direction
                    seq1 = [t(jj,closedM) + bigM * (1 - q(ii,jj) + gamma(ii) + x(ii)) ...
                            - ref_event_ii - t_setup_ii_jj >= 0]: label;
                    seq2 = [t(ii,closedM) + bigM * (q(ii,jj) + gamma(ii) + x(jj)) ...
                            - ref_event_jj - t_setup_jj_ii >= 0]: label;    
                else
                    seq1 = [t(jj,closedM) + bigM * (1 - q(ii,jj) + gamma(jj) + x(ii)) ...
                            - ref_event_ii - t_setup_ii_jj  >= 0]: label;
                    seq2 = [t(ii,closedM) + bigM * (q(ii,jj) + gamma(jj) + x(jj)) ...
                            - ref_event_jj - t_setup_jj_ii  >= 0]: label;
                end
            else
                seq1 = [t(jj,closedM) + bigM * (1 - q(ii,jj) + gamma(ii) + gamma(jj) + x(ii)) ...
                            - ref_event_ii - t_setup_ii_jj >= 0]: label;
                seq2 = [t(ii,closedM) + bigM * (q(ii,jj) + gamma(jj) + gamma(ii) + x(jj)) ...
                            - ref_event_jj - t_setup_jj_ii >= 0]: label;  
            end
            
            seqConstraints = [seqConstraints, seq1, seq2];
%             end



       
%             if direction(jj) == direction(ii)
%                 % We have to fix the q-value!
%                 label = ['fix q(' int2str(ii) ',' int2str(jj) ')'];
%                 if releasetime(ii) < releasetime(jj)
%                     Q = [q(ii,jj) >= 1 - x(ii) - x(jj)] : label;
%                 else
%                     Q = [q(ii,jj) <= x(ii) + x(jj)] : label;
%                 end
%                 seqConstraints = [seqConstraints, Q];
% 
%             end
        end
    end
end
Constraints = [Constraints, seqConstraints];

t_before = settings.TT.blocktimes.setup;
t_after_R = settings.TT.blocktimes.afterR;
t_after_IC = settings.TT.blocktimes.afterIC;
left_end_HW = settings.TT.headway.end_left;
right_end_HW = settings.TT.headway.end_right;
% We still need the minimum separation on the other machines of the
% corridor!
corridorSeq = [];
for ii = 1:Ntrains
    for jj = ii:Ntrains
        if jj>ii
            if direction(ii) == direction(jj)
                % Only if they are in the same direction!
                if direction(ii) == 1
                    for mm = 1:Nmachines-1
                        if bref(mm) > 0
                            try
                                delta_approach_ii_jj = regapproachtime(jj,mm) - approachtime(jj,mm);
                                delta_approach_jj_ii = regapproachtime(ii,mm) - approachtime(ii,mm);
                                t_setup_ii_jj = cleartime(ii,mm) + approachtime(jj,mm) + t_before;
                                t_setup_jj_ii = cleartime(jj,mm) + approachtime(ii,mm) + t_before;
                            catch
                                delta_approach_ii_jj = regapproachtime(jj,mm) - approachtime(jj,1);
                                delta_approach_jj_ii = regapproachtime(ii,mm) - approachtime(ii,1);
                                t_setup_ii_jj = cleartime(ii,1) + approachtime(jj,11) + t_before;
                                t_setup_jj_ii = cleartime(jj,1) + approachtime(ii,1) + t_before;
                            end
                            
                            seqConstraint1 = [t(jj,mm) - t(ii,mm+1) - t_setup_ii_jj ...
                                    + bigM * (1-q(ii,jj)+x(ii)+x(jj)) - delta_approach_ii_jj * gamma(jj) >= 0];   % GAMMA
                            seqConstraint2 = [t(ii,mm) - t(jj,mm+1) - t_setup_jj_ii ...
                                    + bigM * (q(ii,jj)+x(ii)+x(jj)) - delta_approach_jj_ii * gamma(ii) >= 0];   % GAMMA
                            corridorSeq = [corridorSeq seqConstraint1 seqConstraint2];
                        end
                    end
                    % Also for the completing one!

                    delta_approach_ii_jj = regapproachtime(jj,end) - approachtime(jj,end);
                    delta_approach_jj_ii = regapproachtime(ii,end) - approachtime(ii,end);
                    t_setup_ii_jj = cleartime(ii,end) + approachtime(jj,end) + t_before;
                    t_setup_jj_ii = cleartime(jj,end) + approachtime(ii,end) + t_before;
                    % If needed, we give a minimum headway time
                    minHW_right_ii_jj = right_end_HW - processtime(jj,end);
                    minHW_right_jj_ii = right_end_HW - processtime(ii,end);
                    if minHW_right_ii_jj >= t_setup_ii_jj;
                        t_setup_ii_jj = minHW_right_ii_jj;
                        delta_approach_ii_jj = regprocesstime(jj,end) - processtime(jj,end);
                    end
                    if minHW_right_jj_ii >= t_setup_jj_ii; 
                        t_setup_jj_ii = minHW_right_jj_ii;
                        delta_approach_jj_ii = regprocesstime(ii,end) - processtime(ii,end);
                    end
                    
                    if closedM ~= Nmachines
                        seqConstraint1 = [t(jj,end) - C(ii) - t_setup_ii_jj ...
                                + bigM * (1-q(ii,jj)+x(ii)+x(jj)) - delta_approach_ii_jj * gamma(jj) >= 0];   % GAMMA
                        seqConstraint2 = [t(ii,end) - C(jj) - t_setup_jj_ii ...
                                + bigM * (q(ii,jj)+x(ii)+x(jj)) - delta_approach_jj_ii * gamma(ii) >= 0];   % GAMMA
                        corridorSeq = [corridorSeq seqConstraint1 seqConstraint2];
                    end

                else
                    % Trains go in the other direction!

%                         t_setup = cleartime(ii) + t_before + processtime(jj,end);
                    for mm = Nmachines:-1:2
                        if bref(mm)>0
                            try
                                delta_approach_ii_jj = regapproachtime(jj,mm) - approachtime(jj,mm);
                                delta_approach_jj_ii = regapproachtime(ii,mm) - approachtime(ii,mm);
                                t_setup_ii_jj = cleartime(ii,mm) + approachtime(jj,mm) + t_before;
                                t_setup_jj_ii = cleartime(jj,mm) + approachtime(ii,mm) + t_before;
                            catch
                                delta_approach_ii_jj = regapproachtime(jj,mm) - approachtime(jj,1);
                                delta_approach_jj_ii = regapproachtime(ii,mm) - approachtime(ii,1);
                                t_setup_ii_jj = cleartime(ii,1) + approachtime(jj,11) + t_before;
                                t_setup_jj_ii = cleartime(jj,1) + approachtime(ii,1) + t_before;
                            end

                            seqConstraint1 = [t(jj,mm) - t(ii,mm-1) - t_setup_ii_jj ...
                                    + bigM * (1-q(ii,jj)+x(ii)+x(jj)) - delta_approach_ii_jj * gamma(jj) >= 0];   % GAMMA
                            seqConstraint2 = [t(ii,mm) - t(jj,mm-1) - t_setup_jj_ii ...
                                    + bigM * (q(ii,jj)+x(ii)+x(jj)) - delta_approach_jj_ii * gamma(ii) >= 0];   % GAMMA
                            corridorSeq = [corridorSeq seqConstraint1 seqConstraint2];
                        end
                    end
                        
                    % Also for the completing one!
                    delta_approach_ii_jj = regapproachtime(jj,1) - approachtime(jj,1);
                    delta_approach_jj_ii = regapproachtime(ii,1) - approachtime(ii,1);
                    t_setup_ii_jj = cleartime(ii,1) + approachtime(jj,1) + t_before;
                    t_setup_jj_ii = cleartime(jj,1) + approachtime(ii,1) + t_before;
                    % If needed, we give a minimum headway time
                    minHW_left_ii_jj = left_end_HW - processtime(ii,1);
                    if minHW_left_ii_jj >= t_setup_ii_jj
                        t_setup_ii_jj = minHW_left_ii_jj;
                        delta_approach_ii_jj = regprocesstime(ii,1) - processtime(ii,1);
                    end
                    minHW_left_jj_ii = left_end_HW - processtime(jj,1);
                    if minHW_left_jj_ii >= t_setup_jj_ii
                        t_setup_jj_ii = minHW_left_jj_ii;
                        delta_approach_jj_ii = regprocesstime(jj,1) - processtime(jj,1);
                    end
                    
                    if closedM ~= 1
                        seqConstraint1 = [t(jj,1) - C(ii) - t_setup_ii_jj ...
                                + bigM * (1-q(ii,jj)+x(ii)+x(jj)) - delta_approach_ii_jj * gamma(jj) >= 0];   % GAMMA
                        seqConstraint2 = [t(ii,1) - C(jj) - t_setup_jj_ii ...
                                + bigM * (q(ii,jj)+x(ii)+x(jj)) - delta_approach_jj_ii * gamma(ii) >= 0];   % GAMMA
                        corridorSeq = [corridorSeq seqConstraint1 seqConstraint2];
                    end

                end
            end
        end
    end
end
Constraints = [Constraints, corridorSeq];

% But this means that we may loose track of the order on the corridor,
% which we need for the consecutivity constraints! 
% ==> fix q_ij if gamma_j crosses the border!
for ii = 1:Ntrains
    for jj = 1:Ntrains
        if jj > ii
            fixGamma = [q(ii,jj) - gamma(jj) + gamma(ii) >= 0];
            Constraints = [Constraints fixGamma];
            
            if settings.TT.orderlimit_after
                fixGamma = [q(ii,jj) - gamma(jj) + gamma(ii) <= 1];
                Constraints = [Constraints fixGamma];
            end
%         elseif jj < ii && settings.TT.orderlimit_after       % Otherwise it is not needed
%             % What to do with this one????
%             if releasetime(ii) < t_end && releasetime(jj) < t_end
%                 fixGamma = [q(ii,jj) - gamma(jj) + gamma(ii) <= 1];
%                 Constraints = [Constraints fixGamma];
%             end
        end
    end
end
% If a train gets cancelled, set q_ij to 1.
Qs = [];
for ii = 1:Ntrains
    for jj = ii+1:Ntrains
        cancelledQ1 = [q(ii,jj) >= x(ii)];
        cancelledQ2 = [1 - q(ii,jj) >= x(jj) - x(ii)];
        Qs = [Qs cancelledQ1, cancelledQ2];
    end
end
Constraints = [Constraints, Qs];

% FIX ONE SPECIFIC Q VALUE AND IT MAY NOT BE CANCELLED
noCancel1 = [x(1) == 0];
noCancel2 = [x(6) == 0];
fixQ = [q(1,6) == 0];
Constraints = [Constraints, fixQ, noCancel2, noCancel1];


% %% Cancel all
% cancelAll = [x == 1];
% Constraints = [Constraints cancelAll];
% %% Cancel all
% noCancel = [x == 0];
% Constraints = [Constraints noCancel];
% settings.disruption.offbalance = 30;

% Off-balance constraint OR minLOS constraints
minLOS_dir0 = settings.disruption.minLOS_dir0;
minLOS_dir1 = settings.disruption.minLOS_dir1;
if minLOS_dir0 < 0 && minLOS_dir1 < 0
	if settings.disruption.aggregateOB
		dir1 = find(mod(trains,10)==1);
		dir0 = setdiff([1:Ntrains],dir1);
		delta = settings.disruption.duration * settings.disruption.offbalance;
		LBbalance = [sum(x(dir1)) - sum(x(dir0)) >= -delta]: 'LB_balance';
		UBbalance = [sum(x(dir1)) - sum(x(dir0)) <= delta]: 'UB_balance';
		Constraints = [Constraints, LBbalance, UBbalance];
	else
		% If we want this per hour
		if settings.disruption.OBrollinghorizon
			balanceConstraints = generateOffBalanceConstraintsRollingHorizon(trains, releasetime, x, settings);
			Constraints = [Constraints balanceConstraints];
		else
			delta = settings.disruption.offbalance;
			for hh = 1:settings.disruption.duration
				dir1 = find((mod(trains,10)==1) & (floor(trains/1000)==(hh-1)));
				dir0 = find((mod(trains,10)==0) & (floor(trains/1000)==(hh-1)));
				LBbalance = [sum(x(dir1)) - sum(x(dir0)) >= -delta]: ['LB_balance_hour_' int2str(hh)];
				UBbalance = [sum(x(dir1)) - sum(x(dir0)) <= delta]: ['UB_balance_hour_' int2str(hh)];
				Constraints = [Constraints, LBbalance, UBbalance];
			end
		end
	end
else
	% Consider the level of service as the minimum number of trains that has to be
	% running in each direction, per hour.
	for hh = 1:settings.disruption.duration
		dir1 = find((mod(trains,10)==1) & (floor(trains/1000)==(hh-1)));
		dir0 = find((mod(trains,10)==0) & (floor(trains/1000)==(hh-1)));
		if minLOS_dir0 >= 0
			LBbalance = [sum(1 - x(dir0)) >= minLOS_dir0]: ['LB_dir0_hour_' int2str(hh)];
			Constraints = [Constraints, LBbalance];
		end
		if minLOS_dir1 >= 0
			LBbalance = [sum(1 - x(dir1)) >= minLOS_dir1]: ['LB_dir1_hour_' int2str(hh)];
			Constraints = [Constraints, LBbalance];
		end
	end
end

% Max. limit on the number of trains in the same direction right after each
% other.
% consecutiveConstraints = generateConsecutiveConstraints(timetable, trains, settings, q, x);
try
    if settings.TT.orderlimit_before && ~settings.TT.orderlimit_after
        consecutiveConstraints = generateConsecutiveConstraints_v3(timetable, trains, settings, q, x);
    elseif settings.TT.orderlimit_before && settings.TT.orderlimit_after
        consecutiveConstraints = generateConsecutiveConstraints_v4(timetable, trains, settings, q, x, gamma);
    else
        consecutiveConstraints = generateConsecutiveConstraints_v2(timetable, trains, settings, q, x);
    end
catch
    consecutiveConstraints = generateConsecutiveConstraints_v2(timetable, trains, settings, q, x);
end

Constraints = [Constraints, consecutiveConstraints];



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
        Objective = Objective + w_delay * w_dir1 * d(tt);
        Objective = Objective + w_dir1 * w_late_entry * e(tt);     % Very small penalty for an entrance delay.
    else
        Objective = Objective + w_dir0 * w_delay * d(tt);
        Objective = Objective + w_dir0 * w_late_entry * e(tt);
    end
end

% Penalize the cancellation of trains
for tt = 1:Ntrains
    if mod(trains(tt),10) == 1
        Objective = Objective + w_cancel * w_dir1 * x(tt);
    else
        Objective = Objective + w_cancel * w_dir0 * x(tt);
    end
end

%% Build and solve the model
ops = sdpsettings('solver','cplex','usex0',warmstart);
% ops = sdpsettings('solver','gurobi','usex0',warmstart);
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
    vgamma = round(value(gamma));
    
    % Delays?
    orig_delays = departures - orig_deptime;
    
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
            
            if vgamma(find(timetable.train_id(arr_ev(ii,jj)) == trains))
                timetable.adjusted_departure(arr_ev(ii,jj)) = arrivals(ii,jj)...
                                + regTT.running(arr_ev(ii,jj));
            else
                timetable.adjusted_departure(arr_ev(ii,jj)) = arrivals(ii,jj)...
                                + timetable.running(arr_ev(ii,jj));
            end
        end
    end
    timetable.adjusted_departure(completion_ev) = departures';
    
    % Update all timings
    new_timetable = updateTimetable_v5(timetable, regular, vgamma, blocksections, settings, trains, cancelled);
    
    % What type of measures do we take?
    t_closed = arrivals(:,closedM);
%     measures = deriveMeasures(timetable, t_closed, cancelled, releasetime, direction, settings);
    [measures, statistics] = deriveMeasures_v2(timetable, t_closed, cancelled, releasetime, direction, settings, orig_delays, solu)

    
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
