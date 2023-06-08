function [new_timetable, solu, traininfo, statistics, measures] = modelCase3_parallelmachines_full(timetable, blocksections, traininfo, settings)

try
    maxD = settings.disruption.maxDelay;
catch
    maxD = 3600;
end
maxD_THA = settings.constraints.maxDelayTHA;
if isempty(maxD)
    maxD_THA = 3600;
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
try
    w_delay = settings.weights.exitDelay;
catch
    w_delay = 1;
end

%% This version also accounts for the increase in blocking times next to the closed section.
freeOrders_0_counter = 0;
freeOrders_1_counter = 0;
freeOrders_0 = [];
freeOrders_1 = [];

% Path GUROBI
pathGUROBI = '/Library/gurobi1000/macos_universal2/matlab';
addpath(genpath(pathGUROBI));
% Path YALMIP files
pathYALMIP = '/Users/thomasvandendorpe/Dropbox/Thesis/Code/YALMIP-master';
addpath(genpath(pathYALMIP));

% Create different a machine for each track. Keep the closed one as
% together.
Nmachines = size(settings.tracks,1);
Ntrains = size(traininfo,2);

% Extract needed information
direction = [traininfo.dir];
type = {traininfo.type};
releasetime = [traininfo.entry];
for tt = 1:Ntrains
    events = traininfo(tt).ev;
    arr_ev(tt) = events(1);
    completion_ev(tt) = events(end);
    for mm = 1:Nmachines
        % Process time of train tt on track/machine mm
        lab = ['track' int2str(mm)];
        processtime(tt,mm) = sum(timetable.(lab).running(events));
    end
    if settings.constraints.deptimeowntrack
        tr = traininfo(tt).track;
        deptime(tt) = min(releasetime(tt) + processtime(tt,tr));
    else
        cols = ~settings.tracks.closed;
        if any(cols)
            deptime(tt) = min(releasetime(tt) + processtime(tt,cols));
        else
            deptime(tt) = min(releasetime(tt) + processtime(tt,:));
        end
    end
    lab = ['track' int2str(traininfo(tt).track)];
    orig_deptime(tt) = timetable.(lab).departure(events(end));
end
% On which tracks can the train be operated? Boolean variable.
z = ones(Ntrains,Nmachines);
trains = unique(timetable.track1.train_id);
for mm = 1:Nmachines
    for tt = 1:Ntrains
        row = find([traininfo(:).id] == trains(tt));
        allowedtrack = ismember(mm,traininfo(row).allowedtracks);
        if settings.tracks.closed(mm) == 1 || ~allowedtrack
            z(tt,mm) = 0;
        end
    end
end
        

% DETERMINE SETUPTIMES ...
setuptimes = createSetupTimeMatrix_case3(timetable, blocksections, settings, traininfo);



%% Create the variables
t = sdpvar(Ntrains,1);          % Starting time of train on any machine
C = sdpvar(Ntrains,1);          % Completion time on the corridor
D = sdpvar(Ntrains,1);          % Delay upon completion
x = binvar(Ntrains,1);          % Cancellation of full train
dev = binvar(Ntrains,1);        % Deviation of train
m = binvar(Ntrains,Nmachines);  % Train is operated by one machine only.
q = binvar(Ntrains,Ntrains,Nmachines);  % Order variables on the machines.

%% Constraints
bigM = (settings.disruption.duration+1) * 3600;  
Constraints = [];

% Create a warm start: use the original track
for tt = 1:Ntrains
    for mm = 1:Nmachines
        row = find([traininfo(:).id] == trains(tt));
        if traininfo(row).track == mm
            assign(m(tt,mm),1);
        else
            assign(m(tt,mm),0);
        end
    end
    assign(x(tt),0);
end
    
    

% Start time >= release time on the first section of the train
for tt = 1:Ntrains
	minStartTime = [t(tt) - releasetime(tt) * (1-x(tt)-dev(tt)) >= 0] : ['start_' int2str(traininfo(tt).id)];
	Constraints = [Constraints, minStartTime];
end

% Is deviation possible?
if ~settings.deviation.available
	noDev = [dev == 0];
	Constraints = [Constraints, noDev];
end

% Assign train to track, or cancel it.
for tt = 1:Ntrains
    assignTrain = [sum(m(tt,:)) + x(tt) + dev(tt) == 1]: ['assign_' int2str(traininfo(tt).id)];
    Constraints = [Constraints, assignTrain];
    % But only if it is allowed to!
    allowedTrack = [m(tt,:) - z(tt,:) <= 0]: ['allowed_' int2str(traininfo(tt).id) '_on_' int2str(mm)];
    Constraints = [Constraints, allowedTrack];
end


% Completion time depends on the machine on which it is operated! Several
% formulations are possible, this one is in line with Manne's.
for mm = 1:Nmachines
    for tt = 1:Ntrains
        processTime = [C(tt) - t(tt) + bigM * (1 - m(tt,mm)) ...
            >= processtime(tt,mm)]: ['p_' int2str(traininfo(tt).id) '_on_' int2str(mm)];
        Constraints = [Constraints, processTime];
    end
end

% Calculate delay
for tt = 1:Ntrains
    duedate = deptime(tt);
    delay = [D(tt) - C(tt) + duedate * (1 - x(tt) - dev(tt)) >= 0]: ...
        ['delay_' int2str(traininfo(tt).id)];
    Constraints = [Constraints, delay];
    
    % Maximum delay
    if strcmp(traininfo(tt).type,'THA')
        maxDelay = [D(tt) <= maxD_THA]: ['maxDelay_' int2str(traininfo(tt).id)];
    else
        maxDelay = [D(tt) <= maxD]: ['maxDelay_' int2str(traininfo(tt).id)];
    end
    Constraints = [Constraints, delay];
end
delayPositive = [D >= 0];
Constraints = [Constraints, delayPositive];

% Sequencing variables: which comes after which?
% Note that some will have to be fixed to 0 or 1 for the same direction!
seqConstraints = [];
for ii = 1:Ntrains
    for jj = ii:Ntrains
        if jj>ii
            ref_event_ii = C(ii);
            ref_event_jj = C(jj);
            
            for mm = 1:Nmachines
                t_setup_ii_jj = setuptimes.disrupted(ii,jj,mm);
                t_setup_jj_ii = setuptimes.disrupted(jj,ii,mm);
                
                seqConstraint1 = [t(jj) - ref_event_ii - t_setup_ii_jj + bigM *...
                        (3 - q(ii,jj,mm) - m(ii,mm) - m(jj,mm)) >= 0];
                seqConstraint2 = [t(ii) - ref_event_jj - t_setup_jj_ii + bigM *...
                        (2 + q(ii,jj,mm) - m(ii,mm) - m(jj,mm)) >= 0];
                seqConstraints = [seqConstraints, seqConstraint1, seqConstraint2];
            
            
                if traininfo(ii).dir == traininfo(jj).dir
                    label = ['fix q(' int2str(traininfo(ii).id) ',' int2str(traininfo(jj).id) ')'];
                    if releasetime(ii) < releasetime(jj)
                        Q = [q(ii,jj,:) >= 1] : label;
                    else
                        Q = [q(ii,jj,:) <= 0] : label;
                    end
                    seqConstraints = [seqConstraints, Q];
                end
            end
        end
    end
end
Constraints = [Constraints, seqConstraints];



% %% Cancel all
% cancelAll = [x == 1];
% Constraints = [Constraints cancelAll];
% %% Cancel none
% noCancel = [x == 0];
% Constraints = [Constraints noCancel];

if settings.disruption.noCancel_dir0
    dir0 = find([traininfo.dir] == 0);
    noCancel_dir0 = [x(dir0) + dev(dir0) == 0];
    Constraints = [Constraints, noCancel_dir0];
end
if settings.disruption.noCancel_dir1
    dir1 = find([traininfo.dir] == 1);
    noCancel_dir1 = [x(dir1) + dev(dir1) == 0];
    Constraints = [Constraints, noCancel_dir1];
end

% Off-balance constraint OR minLOS constraints
% minLOS_dir0 = settings.disruption.minLOS_dir0;
% minLOS_dir1 = settings.disruption.minLOS_dir1;
% if settings.constraints.OB
% 	if settings.disruption.aggregateOB
% 		dir1 = find([traininfo.dir]==1);
% 		dir0 = setdiff([1:Ntrains],dir1);
% 		delta = settings.disruption.duration * settings.disruption.offbalance;
% 		LBbalance = [sum(x(dir1)) - sum(x(dir0)) >= -delta]: 'LB_balance';
% 		UBbalance = [sum(x(dir1)) - sum(x(dir0)) <= delta]: 'UB_balance';
% 		Constraints = [Constraints, LBbalance, UBbalance];
% 	else
% 		% If we want this per hour
% 		if settings.disruption.OBrollinghorizon
% 			balanceConstraints = generateOffBalanceConstraintsRollingHorizon(trains, releasetime, x, settings);
% 			Constraints = [Constraints balanceConstraints];
% 		else
% 			delta = settings.disruption.offbalance;
% 			for hh = 1:settings.disruption.duration
% 				dir1 = find((mod(trains,10)==1) & (floor(trains/1000)==(hh-1)));
% 				dir0 = find((mod(trains,10)==0) & (floor(trains/1000)==(hh-1)));
% 				LBbalance = [sum(x(dir1)) - sum(x(dir0)) >= -delta]: ['LB_balance_hour_' int2str(hh)];
% 				UBbalance = [sum(x(dir1)) - sum(x(dir0)) <= delta]: ['UB_balance_hour_' int2str(hh)];
% 				Constraints = [Constraints, LBbalance, UBbalance];
% 			end
% 		end
% 	end
% elseif settings.constraints.minLOS
% 	% Consider the level of service as the minimum number of trains that has to be
% 	% running in each direction, per hour.
% 	for hh = 1:settings.disruption.duration
% 		dir1 = find(([traininfo.dir]==1) & (floor([traininfo.id]/1000)==(hh-1)));
% 		dir0 = find(([traininfo.dir]==0) & (floor([traininfo.id]/1000)==(hh-1)));
% 		if minLOS_dir0 >= 0
% 			LBbalance = [sum(1 - x(dir0)) >= minLOS_dir0]: ['LB_dir0_hour_' int2str(hh)];
% 			Constraints = [Constraints, LBbalance];
% 		end
% 		if minLOS_dir1 >= 0
% 			LBbalance = [sum(1 - x(dir1)) >= minLOS_dir1]: ['LB_dir1_hour_' int2str(hh)];
% 			Constraints = [Constraints, LBbalance];
% 		end
% 	end
% end

% Max. number of deviations per hour
if settings.deviation.available
	for hh = 1:(settings.disruption.duration + settings.disruption.nr_hours_after)
		dir1 = find(([traininfo.dir]==1) & (floor([traininfo.id]/1000)==(hh-1)));
		dir0 = find(([traininfo.dir]==0) & (floor([traininfo.id]/1000)==(hh-1)));

		maxDev_dir0 = [sum(dev(dir0)) <= settings.deviation.cap_dir0]: ['dev_dir0_hour_' int2str(hh)];
		maxDev_dir1 = [sum(dev(dir1)) <= settings.deviation.cap_dir1]: ['dev_dir0_hour_' int2str(hh)];
		
		Constraints = [Constraints, maxDev_dir0, maxDev_dir1];
	end
end


%% Build the objective function
Objective = 0;

w_dir0 = settings.weights.pax_dir0 / settings.weights.pax_dir1;
w_dir1 = 1;


% Penalize the delays
for tt = 1:Ntrains
%     Objective = Objective + C(tt) - deptime(tt);
    if traininfo(tt).dir == 1
        Objective = Objective + w_dir1 * w_delay * D(tt);
    else
        Objective = Objective + w_dir0 * w_delay * D(tt);
    end
end

% Penalize the cancellation of trains
for tt = 1:Ntrains
    if traininfo(tt).dir == 1
        Objective = Objective + w_dir1 * w_cancel * x(tt);
    else
        Objective = Objective + w_dir0 * w_cancel * x(tt);
    end
end

% Penalize deviation
w_dev_0 = settings.deviation.delay_dir0;
w_dev_1 = settings.deviation.delay_dir1;
for tt = 1:Ntrains
    if traininfo(tt).dir == 1
        Objective = Objective + w_dir1 * w_dev_1 * dev(tt);
    else
        Objective = Objective + w_dir0 * w_dev_0 * dev(tt);
    end
end

w_changeTrack = settings.weights.trackchange;
if w_changeTrack ~= 0
    for tt = 1:Ntrains
        row = find([traininfo(:).id] == trains(tt));
        tr = traininfo(row).track;
        Objective = Objective + w_changeTrack * (1 - m(tt,tr));
    end
end



%% Build and solve the model
if settings.solver.mipfocus > 0
    ops = sdpsettings('solver','gurobi','usex0',1);
    ops.savesolveroutput = 1;
    ops.showprogress = 1;
    ops.gurobi.MIPFocus = settings.solver.mipfocus;
    if settings.solver.timelimit > 0
        ops.gurobi.TimeLimit = settings.solver.timelimit;
    end
else
    ops = sdpsettings('solver','cplex','usex0',1);
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
    deviated = round(value(dev));
    massigned = round(value(m));
    
    % Delays?
    orig_delays = departures - orig_deptime';
    extra_delays = departures - deptime';
    
    for tt = 1:Ntrains
        row = find([traininfo(:).id] == trains(tt));
        traininfo(row).cancelled = cancelled(tt);
        traininfo(row).deviated = deviated(tt);
        if ~(traininfo(row).cancelled || traininfo(row).deviated)
            traininfo(row).delay = orig_delays(tt);
            traininfo(row).extra_delay = extra_delays(tt);
        else
            traininfo(row).delay = 0;
            traininfo(row).extra_delay = 0;
        end
        traininfo(row).adjusted_entry = arrivals(tt);
        traininfo(row).adjusted_exit = departures(tt);
        mm = find(massigned(tt,:) == 1);
        if ~isempty(mm)
            traininfo(row).newtrack = mm;
        else
            traininfo(row).newtrack = 0;
        end
    end
    

    % Update all timings
    new_timetable = updateTimetable_case3(timetable, settings, traininfo);
    
    
	
    % What type of measures do we take?
    [measures, statistics] = deriveMeasures_case3(new_timetable, traininfo, settings, solu);

	statistics.OFvalue = value(Objective);
    

end

end
