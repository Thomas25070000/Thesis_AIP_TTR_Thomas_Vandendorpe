% This function generates a MILP model based on the job shop scheduling problem
% for a single machine. It then solves the model and returns the resulting timetable,
% solution, measures, and statistics.

function [new_timetable, solu, measures, statistics] = modelCase1_singleMachine_full(timetable, regular, blocksections, trains, minHW, settings, runningtimes)

% Set default values for warmstart, maxD, maxD_THA, w_cancel, w_late_entry,and w_delay
warmstart = 0;
maxD = 3600;
maxD_THA = settings.constraints.maxDelayTHA;
if isempty(maxD)
    maxD_THA = 3600;
end
w_cancel = 3600;
if isempty(w_cancel)
    w_cancel = 3600;
end
w_late_entry = 1/1000;
if isempty(w_late_entry)
    w_late_entry = 1/1000;
end
w_delay = 1;
if isempty(w_delay)
    w_delay = 1;
end
WIC = 0;

% Overwrite default values if specified in settings
try
    warmstart = settings.useWarmstart;
catch
    % Use default value
end
try
    maxD = settings.disruption.maxDelay;
catch
    % Use default value
end
try
    w_cancel = settings.weights.cancel;
catch
    % Use default value
end
try
    w_late_entry = settings.weights.entranceDelay;
catch
    % Use default value
end
try
    w_delay = settings.weights.exitDelay;
catch
    % Use default value
end
try
    WIC = settings.disruption.WICallowed;
catch
    % Use default value
end

% Set w_late_entry to 0 if WIC is not allowed
if ~WIC
    w_late_entry = 0;
end

%% This version also accounts for the increase in blocking times next to the closed section.
max_E = 5000;
% Penalty for cancellation
%w_cancel = 2400;
%w_late_entry = 1/1000;

regTT = regular.TT;
regHW = regular.HW;
freeOrders_0_counter = 0;
freeOrders_1_counter = 0;
freeOrders_0 = [];
freeOrders_1 = [];

% Path CPLEX
pathCPLEX = '/Applications/CPLEX_Studio_Community2211/cplex/bin/x86-64_osx';
addpath(genpath(pathCPLEX));
% Path GUROBI
pathGUROBI = '/Library/gurobi1000/macos_universal2/matlab';
addpath(genpath(pathGUROBI));
% Path YALMIP files
pathYALMIP = '/Users/thomasvandendorpe/Dropbox/Thesis/Code/YALMIP-master';
addpath(genpath(pathYALMIP));
% pathYALMIP = 'D:\OneDrive\Documenten\Thesis\Models\YALMIP';
% addpath(genpath(pathYALMIP))



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
approachtime_end = zeros(Ntrains,1);
approach_other = zeros(Ntrains,1);          % If in the other direction, i.e. on the closed section.
regapproachtime = zeros(Ntrains,Nmachines);
regapproachtime_end = zeros(Ntrains,1);
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
        if bref(bb) > 0 && ismember(block, allev.blocksection)
            processtime(tt,bb) = allev.running(find(allev.blocksection == block));
            regprocesstime(tt,bb) = all_reg_ev.running(find(all_reg_ev.blocksection == block));
            arr_ev(tt,bb) = allev.event_id(find(allev.blocksection == block));
        else
           % processtime(tt,bb) = inf;
            arr_ev(tt,bb) = ev.event_id(1);
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
    
    if ev.arrival(1) >= settings.disruption.duration*3600
        outside_disruption_period(tt) = 1;
    end
end

% DETERMINE SETUPTIMES ...
if settings.disruption.signals
    [setuptimes, proctimes] = createSetupTimeMatrix(timetable, regTT, blocksections, settings);
else
    [setuptimes, proctimes] = createSetupTimeMatrixNoSignals(timetable, regTT, blocksections, settings);
end
    
% Fill the approach times
for tt = 1:Ntrains
    % Depending on direction!
    if direction(tt) == 12 | direction(tt) == 13
        for mm = 1:Nmachines
            if mm == 1
            	[approachtime(tt,1),regapproachtime(tt,1)] = returnFirstApproachTime(blocksections,settings, traintype(tt), mod(direction(tt),100));
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
                [approachtime(tt,end),regapproachtime(tt,end)] = returnFirstApproachTime(blocksections,settings, traintype(tt), mod(direction(tt),100));
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
% t = intvar(Ntrains,Nmachines);
% Completion time on the corridor
C = sdpvar(Ntrains,1);
% C = intvar(Ntrains,1);
% Cancellation of full train
x = binvar(Ntrains,1);
% Order variables on the machines
q = binvar(Ntrains,Ntrains);
% Delay upon completion
d = sdpvar(Ntrains,1);
% Deviation from original release ('entrance') time
e = sdpvar(Ntrains,1);
% Binaries to indicate whether the first arrival on the corridor happens
% after the end of the disruption or not.
gamma = binvar(Ntrains,1);
% Deviation variables
dev = binvar(Ntrains,1);

if warmstart
    assign(x,0);
end

%% Constraints
if settings.solver.adjustBigM
    bigM = setuptimes.maxSetup * 20;
else
    bigM = (settings.disruption.duration+1) * 3600;     % Just make it very big!
end

Constraints = [];
% Start time >= release time on the first section of the train + indicator
% whether it is after the disruption ending.
for tt = 1:Ntrains
	if  direction(tt) == 12 | direction(tt) == 13
		% First block is the left one
		fb = 1;
		if settings.disruption.gammAtStart
			sw = fb;
		else
			sw = S1;
		end
	else
		% First block is at the end.
		fb = Nmachines;
		if settings.disruption.gammAtStart
			sw = fb;
		else
			sw = S2;
		end
	end
	minStartTime = [t(tt,fb) - releasetime(tt) * (1-x(tt)-dev(tt)) >= 0] : ['start_' int2str(trains(tt))];
	% Refer to the disruption end at the time of entry on the switches!
	% This is the last moment on which the route can be changed.
	if settings.disruption.nr_hours_after > 0 && settings.constraints.allowbordercross ...
            && ~settings.constraints.nogamma
		disruptionend1 = [t(tt,sw) - t_end * gamma(tt) >= 0]: ['end_of_disruption_' int2str(trains(tt))];
		disruptionend2 = [t(tt,sw) - bigM * gamma(tt) <= t_end]: ['end_of_disruption_' int2str(trains(tt))];
		Constraints = [Constraints, minStartTime, disruptionend1, disruptionend2];
	else
		Constraints = [Constraints, minStartTime];
	end
	
	% Fix gamma?
	if settings.disruption.nr_hours_after > 0 && ~settings.constraints.nogamma
		if releasetime(tt) + w_cancel < t_end
			% Never able to cross the disruption border, better to cancel it
			% than delay it that much.
			 fixGamma = [gamma(tt) == 0]: ['fix_gamma_' int2str(trains(tt))];
			 Constraints = [Constraints, fixGamma];
		 elseif releasetime(tt) >= t_end
			 fixGamma = [gamma(tt) >= 1 - x(tt) - dev(tt)]: ['fix_gamma_' int2str(trains(tt))];
			 Constraints = [Constraints, fixGamma];
        end
        if settings.constraints.allowbordercross == 0
            if releasetime(tt) < t_end
                fixGamma = [gamma(tt) == 0];
                Constraints = [Constraints, fixGamma];
            else
                % If the trains is cancelled, gamma is always 0!
                fixGamma = [gamma(tt) + x(tt) + dev(tt) <= 1];
                Constraints = [Constraints, fixGamma];
            end
        else
            % If the trains is cancelled, gamma is always 0!
            fixGamma = [gamma(tt) + x(tt) + dev(tt) <= 1];
            Constraints = [Constraints, fixGamma];
        end
    else
        for tt = 1:Ntrains
            if releasetime(tt) < t_end
                fixGamma = [gamma(tt) == 0];
                Constraints = [Constraints, fixGamma];
            end
        end
                
		
	end
end




% Time to start on the next blocksection > time on previous one + process
% time on previous one.
for tt = 1:Ntrains   
    if direction(tt) == 12
        % Running from left to right.
        for bb = [2,3,4,5]

            delta_proc = regprocesstime(tt,bb-1) - processtime(tt,bb-1);
            minStartTime = [t(tt,bb) - t(tt,bb-1) - processtime(tt,bb-1) ...
                * (1-x(tt)-dev(tt)) - delta_proc * gamma(tt) >= 0]: ['intermediate_' ...
                int2str(trains(tt)) '_' int2str(bb)];

            Constraints = [Constraints, minStartTime];
        end
    end
    if direction(tt) == 13
        % Running from left to right.
        for bb = [6,7,8,9]
            if bb == 6
                delta_proc = regprocesstime(tt,1) - processtime(tt,1);
                minStartTime = [t(tt,bb) - t(tt,1) - processtime(tt,1) ...
                    * (1-x(tt)-dev(tt)) - delta_proc * gamma(tt) >= 0]: ['intermediate_' ...
                    int2str(trains(tt)) '_' int2str(bb)]; 
            else
                delta_proc = regprocesstime(tt,bb-1) - processtime(tt,bb-1);
                minStartTime = [t(tt,bb) - t(tt,bb-1) - processtime(tt,bb-1) ...
                    * (1-x(tt)-dev(tt)) - delta_proc * gamma(tt) >= 0]: ['intermediate_' ...
                    int2str(trains(tt)) '_' int2str(bb)];
            end
            Constraints = [Constraints, minStartTime];
        end
    end
    if direction(tt) == 02
        % From right to left, start on one-but-last.
        for bb = [4,3,2,1]

            delta_proc = regprocesstime(tt,bb+1) - processtime(tt,bb+1);        
            minStartTime = [t(tt,bb) - t(tt,bb+1) - processtime(tt,bb+1) ...
                * (1-x(tt)-dev(tt)) - delta_proc * gamma(tt) >= 0]: ['intermediate_' ...
                int2str(trains(tt)) '_' int2str(bb)];    

            Constraints = [Constraints, minStartTime];
        end
    end
    if direction(tt) == 03
        % From right to left, start on one-but-last.
        for bb = [8,7,6,1]
            if bb == 1
                delta_proc = regprocesstime(tt,6) - processtime(tt,6);
                minStartTime = [t(tt,bb) - t(tt,6) - processtime(tt,6) ...
                    * (1-x(tt)-dev(tt)) - delta_proc * gamma(tt) >= 0]: ['intermediate_' ...
                    int2str(trains(tt)) '_' int2str(bb)];  
            else
                delta_proc = regprocesstime(tt,bb+1) - processtime(tt,bb+1);
                minStartTime = [t(tt,bb) - t(tt,bb+1) - processtime(tt,bb+1) ...
                    * (1-x(tt)-dev(tt)) - delta_proc * gamma(tt) >= 0]: ['intermediate_' ...
                    int2str(trains(tt)) '_' int2str(bb)];     
            end
            Constraints = [Constraints, minStartTime];
        end
    end
end

% Completion time > start time on one but last + process time
% What is the delay?
for tt = 1:Ntrains
    if direction(tt) == 12
        % Last event takes place on the last machine.
        delta_proc = processtime(tt,7) - regprocesstime(tt,7);
        minCompletionTime = [C(tt) - t(tt,7) - processtime(tt,7) * (1-x(tt)-dev(tt))...
             + delta_proc * gamma(tt) >= 0]: ['completion_' int2str(trains(tt))];
    end
    if direction(tt) == 13
        % Last event takes place on the last machine.
        delta_proc = processtime(tt,end) - regprocesstime(tt,end);
        minCompletionTime = [C(tt) - t(tt,end) - processtime(tt,end) * (1-x(tt)-dev(tt))...
             + delta_proc * gamma(tt) >= 0]: ['completion_' int2str(trains(tt))];
    end
    if direction(tt) == 02 || direction(tt) == 03
        delta_proc = processtime(tt,1) - regprocesstime(tt,1);
        % Last event takes place on the first machine.
        minCompletionTime = [C(tt) - t(tt,1) - processtime(tt,1) * (1-x(tt)-dev(tt)) ...
            + delta_proc * gamma(tt) >= 0]: ['completion_' int2str(trains(tt))];
    end
    if releasetime(tt) >= t_end
        delay = [d(tt) + orig_deptime(tt) * (1-x(tt)-dev(tt)) - C(tt) >= 0]: ['delay_' int2str(trains(tt))];
    else
        delay = [d(tt) + deptime(tt) * (1-x(tt)-dev(tt)) - C(tt) >= 0]: ['delay_' int2str(trains(tt))]; 
    end
	delayPositive = [d(tt) >= 0]: ['positiveDelay_' int2str(trains(tt))];
    if strcmp(traintype(tt),'THA')
        maxDelay = [e(tt) <= maxD_THA]: ['maxDelay_' int2str(trains(tt))];
    else
        maxDelay = [d(tt) <= maxD]: ['maxDelay_' int2str(trains(tt))];
    end
    if direction(tt) == 12 | direction(tt) == 13
        entrancedelay = [e(tt) + releasetime(tt) * (1-x(tt)-dev(tt)) - t(tt,1) >= 0]: ['entrancedelay_' int2str(trains(tt))];
    else
        entrancedelay = [e(tt) + releasetime(tt) * (1-x(tt)-dev(tt)) - t(tt,end) >= 0]: ['entrancedelay_' int2str(trains(tt))];
    end
    minEtt = [e(tt) >= 0];
    Constraints = [Constraints, minCompletionTime, delay, maxDelay, delayPositive, entrancedelay, minEtt];
end


% Sequencing variables: which comes after which?
% Note that some will have to be fixed to 0 or 1 for the same direction!
seqConstraints = [];
closedM = find(bref==0);
t_before = settings.TT.blocktimes.setup;
cleartime = settings.TT.blocktimes.afterIC;
for ii = 1:Ntrains
    for jj = ii:Ntrains
        if jj>ii

            if direction(ii) == 12 | direction(ii) == 13
                if closedM == Nmachines
                    ref_event_ii = C(ii);
                else
                    ref_event_ii = t(ii,closedM+1);
                end
            else
                if closedM == 1
                    ref_event_ii = C(ii);
                else
                    ref_event_ii = t(ii,closedM-1);
                end
            end
            % Ref_event_jj
            if direction(jj) == 12 | direction(jj) == 13
                if closedM == Nmachines
                    ref_event_jj = C(jj);
                else
                    ref_event_jj = t(jj,closedM+1);
                end
            else
                if closedM == 1
                    ref_event_jj = C(jj);
                else
                    ref_event_jj = t(jj,closedM-1);
                end
            end            
            mm = closedM;

            t_setup_ii_jj = setuptimes.disrupted(ii,jj,closedM);
            t_setup_jj_ii = setuptimes.disrupted(jj,ii,closedM);
            delta_1gamma_ii_jj = setuptimes.disrupted(ii,jj,closedM)...
                                - setuptimes.onegamma(ii,jj,closedM);
            delta_1gamma_jj_ii = setuptimes.disrupted(jj,ii,closedM)...
                                - setuptimes.onegamma(jj,ii,closedM);
            delta_2gamma_ii_jj = setuptimes.disrupted(ii,jj,closedM)...
                                - setuptimes.bothgammas(ii,jj,closedM);
            delta_2gamma_jj_ii = setuptimes.disrupted(jj,ii,closedM)...
                                - setuptimes.bothgammas(jj,ii,closedM);

            
            % Do not add them in case both are after the hour!
            if releasetime(ii) < t_end && releasetime(jj) < t_end
          
                seqConstraint1 = [t(jj,mm) - ref_event_ii - t_setup_ii_jj ...
                                + bigM * (1-q(ii,jj)+x(ii)+x(jj)+dev(ii)+dev(jj)) ...
                                + delta_1gamma_ii_jj * (gamma(jj) - gamma(ii)) ...
                                + delta_2gamma_ii_jj * gamma(ii) >= 0];   % GAMMA
                seqConstraint2 = [t(ii,mm) - ref_event_jj - t_setup_jj_ii ...
                                + bigM * (q(ii,jj)+x(ii)+x(jj)+dev(ii)+dev(jj))...
                                + delta_1gamma_jj_ii * (gamma(ii) - gamma(jj)) ...
                                + delta_2gamma_jj_ii * gamma(jj) >= 0];   % GAMMA
                seqConstraints = [seqConstraints, seqConstraint1, seqConstraint2];
            end
            % Fixing the orders between trains!
            % CRITERIA
            flag_fixOrder = 0;
            threshold = settings.constraints.fixEntranceOrderThreshold;
            if abs(direction(ii)-direction(jj))<2
                entrancediff = releasetime(jj) - releasetime(ii);
                tjourney_ii = deptime(ii) - releasetime(ii);
                tjourney_jj = deptime(jj) - releasetime(jj);
                if settings.constraints.fixEntranceOrder
                    flag_fixOrder = 1;
                elseif abs(entrancediff) >= threshold
                    % Time difference is very large!
                    flag_fixOrder = 1;
                elseif tjourney_ii == tjourney_jj && settings.constraints.fixEntranceOrderSameTime
                    flag_fixOrder = 1;
                else
                    flag_fixOrder = 0;
                end
            else
                % Different directions!
                flag_fixOrder = 0;
            end
            
            if flag_fixOrder
                % We have to fix the q-value!
                label = ['fix q(' int2str(ii) ',' int2str(jj) ')'];
                if releasetime(ii) < releasetime(jj)
                    Q = [q(ii,jj) >= 1 - x(ii) - x(jj) - dev(ii) - dev(jj)] : label;
                else
                    Q = [q(ii,jj) <= x(ii) + x(jj) + dev(ii) + dev(jj)] : label;
                end
                seqConstraints = [seqConstraints, Q];
            elseif abs(direction(ii)-direction(jj))<2

                if direction(jj) == 12 | direction(jj) == 13
                    freeOrders_1_counter = freeOrders_1_counter + 1;
                    freeOrders_1 = [freeOrders_1; ii jj];
                else
                    freeOrders_0_counter = freeOrders_0_counter + 1;
                    freeOrders_0 = [freeOrders_0; ii jj];
                end
            end
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
            if abs(direction(ii)-direction(jj))<2
                % Only if they are in the same direction!
                if direction(ii) == 12 | direction(ii) == 13
                     for mm = 1:Nmachines
                        if mm == Nmachines
                            ref_event_ii = C(ii);
                            ref_event_jj = C(jj);
                            delta_1gamma_ii_jj = setuptimes.disrupted(ii,jj,mm)...
                                                - setuptimes.onegamma(ii,jj,mm);
                            delta_1gamma_jj_ii = setuptimes.disrupted(jj,ii,mm)...
                                                - setuptimes.onegamma(jj,ii,mm);
                            delta_2gamma_ii_jj = setuptimes.disrupted(ii,jj,mm)...
                                                - setuptimes.bothgammas(ii,jj,mm);
                            delta_2gamma_jj_ii = setuptimes.disrupted(jj,ii,mm)...
                                                - setuptimes.bothgammas(jj,ii,mm);
                        else
                            ref_event_ii = t(ii,mm+1);
                            ref_event_jj = t(jj,mm+1);
                            delta_1gamma_ii_jj = setuptimes.disrupted(ii,jj,mm)...
                                                - setuptimes.onegamma(ii,jj,mm);
                            delta_1gamma_jj_ii = setuptimes.disrupted(jj,ii,mm)...
                                                - setuptimes.onegamma(jj,ii,mm);
                            delta_2gamma_ii_jj = setuptimes.disrupted(ii,jj,mm)...
                                                - setuptimes.bothgammas(ii,jj,mm);
                            delta_2gamma_jj_ii = setuptimes.disrupted(jj,ii,mm)...
                                                - setuptimes.bothgammas(jj,ii,mm);
                        end
                
                        t_setup_ii_jj = setuptimes.disrupted(ii,jj,mm);
                        t_setup_jj_ii = setuptimes.disrupted(jj,ii,mm);
                       
					    seqConstraint1 = [t(jj,mm) - ref_event_ii - t_setup_ii_jj ...
								+ bigM * (1-q(ii,jj)+x(ii)+x(jj)+dev(ii)+dev(jj)) ...
                                + delta_1gamma_ii_jj * (gamma(jj) - gamma(ii)) ...
                                + delta_2gamma_ii_jj * gamma(ii) >= 0];   % GAMMA
						seqConstraint2 = [t(ii,mm) - ref_event_jj - t_setup_jj_ii ...
								+ bigM * (q(ii,jj)+x(ii)+x(jj)+dev(ii)+dev(jj))...
                                + delta_1gamma_jj_ii * (gamma(ii) - gamma(jj)) ...
                                + delta_2gamma_jj_ii * gamma(jj) >= 0];   % GAMMA
						corridorSeq = [corridorSeq seqConstraint1 seqConstraint2];
                    end
                else
                    % Trains go in the other direction!
					for mm = Nmachines:-1:1
                        if mm == 1
                            ref_event_ii = C(ii);
                            ref_event_jj = C(jj);
                            delta_2gamma_ii_jj = setuptimes.disrupted(ii,jj,mm)...
                                                - setuptimes.bothgammas(ii,jj,mm);
                            delta_2gamma_jj_ii = setuptimes.disrupted(jj,ii,mm)...
                                                - setuptimes.bothgammas(jj,ii,mm);
                        else
                            ref_event_ii = t(ii,mm-1);
                            ref_event_jj = t(jj,mm-1);
                            delta_2gamma_ii_jj = setuptimes.disrupted(ii,jj,mm)...
                                                - setuptimes.bothgammas(ii,jj,mm);
                            delta_2gamma_jj_ii = setuptimes.disrupted(jj,ii,mm)...
                                                - setuptimes.bothgammas(jj,ii,mm);
                        end
                        t_setup_ii_jj = setuptimes.disrupted(ii,jj,mm);
                        t_setup_jj_ii = setuptimes.disrupted(jj,ii,mm);
                        delta_1gamma_ii_jj = setuptimes.disrupted(ii,jj,mm)...
                                                - setuptimes.onegamma(ii,jj,mm);
                        delta_1gamma_jj_ii = setuptimes.disrupted(jj,ii,mm)...
                                                - setuptimes.onegamma(jj,ii,mm);  


                        seqConstraint1 = [t(jj,mm) - ref_event_ii - t_setup_ii_jj ...
								+ bigM * (1-q(ii,jj)+x(ii)+x(jj)+dev(ii)+dev(jj)) ...
                                + delta_1gamma_ii_jj * (gamma(jj) - gamma(ii)) ...
                                + delta_2gamma_ii_jj * gamma(ii) >= 0];   % GAMMA
						seqConstraint2 = [t(ii,mm) - ref_event_jj - t_setup_jj_ii ...
								+ bigM * (q(ii,jj)+x(ii)+x(jj)+dev(ii)+dev(jj))...
                                + delta_1gamma_jj_ii * (gamma(ii) - gamma(jj)) ...
                                + delta_2gamma_jj_ii * gamma(jj) >= 0];   % GAMMA
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
if  settings.constraints.maxConsecutive
    Qs = [];
    for ii = 1:Ntrains
        for jj = ii+1:Ntrains
            cancelledQ1 = [q(ii,jj) >= x(ii) + dev(ii)];
            cancelledQ2 = [1 - q(ii,jj) >= x(jj) - x(ii) + dev(jj) - dev(ii)];
            Qs = [Qs cancelledQ1, cancelledQ2];
        end
    end
    Constraints = [Constraints, Qs];
end



if settings.disruption.noCancel_dir0
    dir0 = find(direction == 02 | direction == 03);
    noCancel_dir0 = [x(dir0) == 0];
    Constraints = [Constraints, noCancel_dir0];
end
if settings.disruption.noCancel_dir1
    dir1 = find(direction == 12 | direction == 13);
    noCancel_dir1 = [x(dir1) == 0];
    Constraints = [Constraints, noCancel_dir1];
end


%% Build the objective function
Objective = 0;

w_dir0 = settings.weights.pax_dir0 / settings.weights.pax_dir1;
w_dir1 = 1;


% Penalize the delays
for tt = 1:Ntrains
%     Objective = Objective + C(tt) - deptime(tt);
    if mod(trains(tt),100) == 12 || mod(trains(tt),100) == 13
        Objective = Objective + w_dir1 * w_delay * d(tt);
        Objective = Objective + w_dir1 * w_late_entry * e(tt);     % Very small penalty for an entrance delay.
    else
        Objective = Objective + w_dir0 * w_delay * d(tt);
        Objective = Objective + w_dir0 * w_late_entry * e(tt);
    end
end

% Penalize the cancellation of trains
for tt = 1:Ntrains
    if mod(trains(tt),100) == 12 || mod(trains(tt),100) == 13
        Objective = Objective + w_dir1 * w_cancel * x(tt);
    else
        Objective = Objective + w_dir0 * w_cancel * x(tt);
    end
end

% Penalize deviation
w_dev_0 = settings.deviation.delay_dir0;
w_dev_1 = settings.deviation.delay_dir1;
for tt = 1:Ntrains
    if mod(trains(tt),100) == 12 || mod(trains(tt),100) == 13
        Objective = Objective + w_dir1 * w_dev_1 * dev(tt);
    else
        Objective = Objective + w_dir0 * w_dev_0 * dev(tt);
    end
end

%% Build and solve the model
if settings.solver.mipfocus > 0
    ops = sdpsettings('solver','gurobi','usex0',warmstart);
    ops.savesolveroutput = 1;
    ops.showprogress = 1;
    ops.gurobi.MIPFocus = settings.solver.mipfocus;
    if settings.solver.timelimit > 0
        ops.gurobi.TimeLimit = settings.solver.timelimit;
    end
else
    ops = sdpsettings('solver','gurobi','usex0',warmstart);
    ops.savesolveroutput = 1;
    ops.showprogress = 1;
    if settings.solver.timelimit > 0
        ops.cplex.timelimit = settings.solver.timelimit;
    end
end

[solu] = optimize(Constraints,Objective,ops);
solu.nr_freeOrders_0 = freeOrders_0_counter;
solu.nr_freeOrders_1 = freeOrders_1_counter;
if settings.solver.timelimit > 0
    solu.timelimit = settings.solver.timelimit;
end

if strcmp(ops.solver,'cplex')
    status = solu.solveroutput.output.cplexstatus;
elseif strcmp(ops.solver,'gurobi')
    statusT = solu.solveroutput.result.status;
    status = 0;
end
    
disp(int2str(status));

if status == 103 
    disp('Infeasible!')
    pause();
else
    arrivals =  round(value(t));
    departures = round(value(C));
    cancelled = round(value(x));
    vgamma = round(value(gamma));
    deviated = round(value(dev));
    
    % Delays?
    orig_delays = departures - orig_deptime;
    extra_delays = departures - deptime;
    
    cancelled_ev = [];
    deviated_ev = [];
    for xx = 1:length(cancelled)
        if cancelled(xx)
            rows = find(timetable.train_id == trains(xx));
            cancelled_ev = [cancelled_ev, rows'];
        elseif deviated(xx)
            rows = find(timetable.train_id == trains(xx));
            deviated_ev = [deviated_ev, rows'];
        end
    end
    timetable.cancelled(1) = 0;
    timetable.cancelled(cancelled_ev) = 1;
    timetable.deviated(1) = 0;
    timetable.deviated(deviated_ev) = 1;
    
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
    
    % Identify changed orders
    changedOrders_0 = [];
    if ~isempty(freeOrders_0)
        for ff = 1:size(freeOrders_0,1)
            ii = freeOrders_0(ff,1);
            jj = freeOrders_0(ff,2);
            if ~cancelled(ii) && ~cancelled(jj)
                before = releasetime(ii) < releasetime(jj);
                after = arrivals(ii,end) < arrivals(jj,end);
                if before ~= after
                    changedOrders_0 = [changedOrders_0; trains(ii) trains(jj)];
                end
            end
        end
    end
    changedOrders_1 = [];
    if ~isempty(freeOrders_1)
        for ff = 1:size(freeOrders_1,1)
            ii = freeOrders_1(ff,1);
            jj = freeOrders_1(ff,2);
            if ~cancelled(ii) && ~cancelled(jj)
                before = releasetime(ii) < releasetime(jj);
                after = arrivals(ii,1) < arrivals(jj,1);
                if before ~= after
                    changedOrders_1 = [changedOrders_1; trains(ii) trains(jj)];
                end
            end
        end
    end  
    solu.changedOrders_0 = changedOrders_0;
    solu.changedOrders_1 = changedOrders_1;
   
    % Update all timings
	if settings.disruption.nr_hours_after > 0
		new_timetable = updateTimetable_v5(timetable, regular, vgamma, blocksections, settings, trains, cancelled);
    else
		new_timetable = updateTimetable_v4(timetable, blocksections, settings, trains, cancelled);
	end
	
    % What type of measures do we take?
    t_closed = arrivals(:,closedM);
%     measures = deriveMeasures(timetable, t_closed, cancelled, releasetime, direction, settings);
   % [measures, statistics] = deriveMeasures_v2(timetable, t_closed, cancelled, releasetime, direction, settings, orig_delays, extra_delays, solu)

    
    % Plot new timetable
%     type = 'hour';
    type = 'complete';
%     type = 'base';
    include_blocks = 1;
    directions = [02 03 12 13];
%     directions = 0;

    HH = 17;
    MM = 50;
    SS = 0;
    firstTime = HH * 3600 + MM * 60 + SS;
    [arrtime, arrtimeHHMMSS, deptime, deptimeHHMMSS] = retrieveArrivalAndDepartureTimes(timetable, firstTime);
    [line_new, blocks_new] = plotTT_full(new_timetable, blocksections, settings, type, include_blocks, directions,firstTime);
    
end

end
