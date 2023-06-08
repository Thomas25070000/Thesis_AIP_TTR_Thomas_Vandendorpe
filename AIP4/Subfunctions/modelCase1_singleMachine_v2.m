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
    releasetime(tt) = ev.arrival(1)
    direction(tt) = ev.direction(1);
    type(tt) = ev.train_type(1);
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
% Order respected variables
l = binvar(Ntrains,Ntrains);
% Delay upon completion
d = sdpvar(Ntrains,1);

%% Constraints
bigM = settings.disruption.duration * 3600 * 20;     % Just make it very big!
Constraints = [];
% Start time >= release time

for tt = 1:Ntrains
    minStartTime = [t(tt) - releasetime(tt) * (1-x(tt)) >= 0] : 'start';
    Constraints = [Constraints, minStartTime];

    % Completion time > start time + process time
    processTime = [C(tt) - t(tt) - processtime(tt) * (1-x(tt)) >= 0]: 'processTime';
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
        label = ['seq_' int2str(ii) '_' int2str(jj)];
        seq1 = [t(jj) + bigM * (2 - q(ii,jj) - l(ii,jj)) - t(ii) ...
                - setuptime(ii,jj,2) * (1 - x(ii)) >= 0]: label;
        seq2 = [t(jj) - bigM * (1 - q(ii,jj) + l(ii,jj)) - t(ii) ...
                 - setuptime(ii,jj,1) * (1 - x(ii)) <= 0]: label;
        seq3 = [t(jj) + bigM * (1 - q(ii,jj) + l(ii,jj)) - t(ii) ...
                 - setuptime(ii,jj,3) * (1 - x(ii)) >= 0]: label;
        seq4 = [t(ii) + bigM * (q(ii,jj)+1-l(ii,jj)) - t(jj) ...
                - setuptime(jj,ii,2) * (1 - x(jj)) >= 0]: label;
        seq5 = [t(ii) - bigM * (q(ii,jj)+l(ii,jj)) - t(jj) ...
                 - setuptime(jj,ii,1) * (1 - x(jj)) <= 0]: label;
        seq6 = [t(ii) + bigM * (q(ii,jj)+l(ii,jj)) - t(jj) ...
                 - setuptime(jj,ii,3) * (1 - x(jj)) >= 0]: label;
        seqConstraints = [seqConstraints, seq1, seq2, seq3, seq4, seq5, seq6];
        %seqConstraints = [seqConstraints, seq6];

        if direction(ii)==direction(jj)
            label = ['fix l(' int2str(ii) ',' int2str(jj) ')'];
            L = [l(ii,jj) == 1] : label;
            % We have to fix the l-value!
            label = ['fix q(' int2str(ii) ',' int2str(jj) ')'];
            seqConstraints = [seqConstraints, L];
        end
%         label = ['fix l(' int2str(ii) ',' int2str(jj) ')'];
%         L = [l(ii,jj) == 1] : label;
%         % We have to fix the l-value!
%         label = ['fix q(' int2str(ii) ',' int2str(jj) ')'];
%         seqConstraints = [seqConstraints, L];

%         if releasetime(ii) < releasetime(jj)
%             Q = [q(ii,jj) >= 1 - x(ii) - x(jj)] : label;
%         else
% %                     Q = [q(ii,jj) == 0] : label;
%             Q = [q(ii,jj) <= x(ii) + x(jj)] : label;
%         end
%         seqConstraints = [seqConstraints, Q];
        if direction(ii) == direction(jj)
            % We have to fix the q-value!
            label = ['fix q(' int2str(ii) ',' int2str(jj) ')'];
            if releasetime(ii) < releasetime(jj)
                Q = [q(ii,jj) >= 1 - x(ii) - x(jj)] : label;
            else
%                     Q = [q(ii,jj) == 0] : label;
                Q = [q(ii,jj) <= x(ii) + x(jj)] : label;
            end
            seqConstraints = [seqConstraints, Q];
        end

        
    end
end
Constraints = [Constraints, seqConstraints];

%% Build the objective function
Objective = 0;

% Penalize the delays
delay_IC = 1;
delay_P = 0.25;
for tt = 1:Ntrains
%     Objective = Objective + C(tt) - deptime(tt);
    if strcmp(type(tt),'IC')
        Objective = Objective + delay_IC*d(tt);
    else
        Objective = Objective + delay_P*d(tt);
    end
end

cancel_IC = 3600;
cancel_P = 900;
% Penalize the cancellation
% of trains
for tt = 1:Ntrains
    if strcmp(type(tt),'IC')
        Objective = Objective + cancel_IC * x(tt);
    else
        Objective = Objective + cancel_P * x(tt);
    end
end

% % Penalize the delays
% for tt = 1:Ntrains
% %     Objective = Objective + C(tt) - deptime(tt);
%     Objective = Objective + d(tt);
% end
% 
% % Penalize the cancellation
% % of trains
% w_cancel = 3600;
% for tt = 1:Ntrains
%     Objective = Objective + w_cancel * x(tt);
% end
% 
% %% Build and solve the model
ops = sdpsettings('solver','gurobi');
ops.savesolveroutput = 1;
ops.gurobi.TimeLimit = 120;

[solu] = optimize(Constraints,Objective,ops);
status = solu.solveroutput.result.status;
disp(int2str(status));

if status == 103
    disp('Infeasible!')
    pause();
else
    arrivals =  round(value(t))
    departures = round(value(C));
    cancelled = round(value(x));
    order = value(q)
    order_switch = value(l)
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
    [line_new, blocks_new] = plotTT_new(new_timetable, blocksections, settings, type, include_blocks);
    [measures, statistics] = deriveMeasures_v2(new_timetable, cancelled, releasetime, direction, settings, orig_delays, extra_delays, solu);
end

end
