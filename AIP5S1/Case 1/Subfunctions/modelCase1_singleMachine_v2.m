function [timetable, solu, statistics] = modelCase1_singleMachine_v2(timetable, blocksections, trains, minHW, settings)

% Path GUROBI
pathGUROBI = '/Library/gurobi1000/macos_universal2/matlab';
addpath(genpath(pathGUROBI));
% Path YALMIP files
pathYALMIP = '/Users/thomasvandendorpe/Dropbox/Thesis/Code/YALMIP-master';
addpath(genpath(pathYALMIP));

%% Process until we have a machine scheduling problem
Ntrains = length(trains);
direction = zeros(Ntrains,1);
% Process times = total running time over the closed section.
processtime = zeros(Ntrains,1);
% Release times = arrival at the start of the closed section.
releasetime = zeros(Ntrains,1);
% Setup times are given by the minimum headways
setuptime = minHW;
% The original departure from the closed part
deptime = zeros(Ntrains,1);
orig_deptime = zeros(Ntrains,1);


closed_blocks = find(blocksections.closed);
for tt = 1:Ntrains
    id = trains(tt);
    ev = timetable(find(timetable.train_id == id),:);
    % Events on the closed section
    ev = ev(find(ismember(ev.blocksection,closed_blocks)),:);
    % Get the needed values
    processtime(tt) = sum(ev.running);
    release_ev(tt) = ev.event_id(1);        % First event of this train on the closed section.
    releasetime(tt) = ev.arrival(1);
    direction(tt) = ev.direction(1);
%     deptime(tt) = ev.departure(end);
    deptime(tt) = ev.arrival(1) + processtime(tt);
    orig_deptime(tt) = ev.departure(end);
    completion_ev(tt) = ev.event_id(end);   % Last event of this train on the closed section.
end


%% Create the variables
Ntrains = length(trains);
% Starting time
t = sdpvar(Ntrains,1);
% Completion time
C = sdpvar(Ntrains,1);
% Cancellation
x = binvar(Ntrains,1);
% Order variables
q = binvar(Ntrains,Ntrains);
% Delay upon completion
d = sdpvar(Ntrains,1);

%% Constraints
bigM = settings.disruption.duration * 3600*20;     % Just make it very big!
Constraints = [];
% Start time >= release time
for tt = 1:Ntrains
    minStartTime = [t(tt) - releasetime(tt)*(1 - x(tt)) >= 0] : 'start';
    Constraints = [Constraints, minStartTime];

    % Completion time > start time + process time
    processTime = [C(tt) - t(tt) - processtime(tt)*(1 - x(tt)) >= 0]: 'processTime';
    Constraints = [Constraints, processTime];
    
    delay = [d(tt) + deptime(tt) * (1 - x(tt)) - C(tt) >= 0]: 'delay';
    Constraints = [Constraints, delay];
    
    % NoCancel
%     noCancel = [x(tt)==0];
%     Constraints = [Constraints, noCancel];
    
end

% Sequencing variables: which comes after which?
% Note that some will have to be fixed to 0 or 1 for the same direction!
seqConstraints = [];
for ii = 1:Ntrains
    for jj = ii:Ntrains
        if jj>ii
            label = ['seq_' int2str(ii) '_' int2str(jj)];
            seq1 = [t(jj) + bigM * (1 - q(ii,jj) + x(jj) + x(ii)) - t(ii) ...
                    - setuptime(ii,jj) >= 0]: label;
            seq2 = [t(ii) + bigM * (q(ii,jj) + x(ii) + x(jj)) - t(jj) ...
                    - setuptime(jj,ii) >= 0]: label;
            seqConstraints = [seqConstraints, seq1, seq2];           
            if releasetime(ii) >= releasetime(jj) && direction(ii)==direction(jj)
                label = ['fix q(' int2str(ii) ',' int2str(jj) ')'];
                Q = [q(ii,jj) == 0] : label;
                seqConstraints = [seqConstraints, Q];
            end
        end
    end
end
Constraints = [Constraints, seqConstraints];

%% Build the objective function
Objective = 0;

% Penalize the delays
for tt = 1:Ntrains
%     Objective = Objective + C(tt) - deptime(tt);
    if tt == 4 || tt == 6 || tt == 14 || tt == 12 || tt==19 || tt==22 || tt == 28 || tt == 29 || tt == 36 || tt == 37
        Objective = Objective + 10000*d(tt);
    else 
        Objective = Objective + d(tt);
    end
end

% Penalize the cancellation of trains
w_cancel = 3600;
for tt = 1:Ntrains
    if tt == 4 || tt == 6 || tt==14 || tt == 12 || tt==19 || tt==22 || tt == 28 || tt == 29 || tt == 36 || tt == 37
        Objective = Objective + w_cancel * 10000* x(tt);
    else
        Objective = Objective + w_cancel * x(tt);
    end
end

%% Build and solve the model
ops = sdpsettings('solver','gurobi');
ops.savesolveroutput = 1;

[solu] = optimize(Constraints,Objective,ops);
status = solu.solveroutput.result.status;
disp(int2str(status));

if status == 103
    disp('Infeasible!')
    pause();
else
    arrivals =  round(value(t));
    departures = round(value(C));
    cancelled = round(value(x));
    orig_delays = departures - orig_deptime;
    extra_delays = departures - deptime;

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
    timetable.adjusted_arrival(release_ev) = arrivals';
    timetable.adjusted_departure(1) = 0;
    timetable.adjusted_departure(completion_ev) = departures';
    
    % Update all timings
    new_timetable = updateTimetable(timetable, blocksections, settings, trains, cancelled);
    
    
    % Plot new timetable
%     type = 'hour';
    type = 'complete';
%     type = 'base';
    include_blocks = 1;
    
    [measures, statistics] = deriveMeasures_v2(new_timetable, cancelled, releasetime, direction, settings, orig_delays, extra_delays, solu);
    [line_new, blocks_new] = plotTT_full(new_timetable, blocksections, settings, type, include_blocks);
    
 end

end
