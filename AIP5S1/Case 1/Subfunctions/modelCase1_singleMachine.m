function [timetable, solu] = modelCase1_singleMachine(timetable, blocksections, trains, minHW, settings)

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
    deptime(tt) = ev.departure(end);
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


%% Constraints
bigM = settings.disruption.duration * 3600 * 2;     % Just make it very big!
Constraints = [];
% Start time >= release time
for tt = 1:Ntrains
    minStartTime = [t(tt) - releasetime(tt) * (1-x(tt)) >= 0] : 'start';
    Constraints = [Constraints, minStartTime];

    processTime = [C(tt) - t(tt) - processtime(tt) * (1-x(tt)) >= 0]: 'processTime';
    Constraints = [Constraints, processTime];
end

% Sequencing variables: which comes after which?
% Note that some will have to be fixed to 0 or 1 for the same direction!
seqConstraints = [];
for ii = 1:Ntrains-1
    for jj = (ii+1):Ntrains
        if ii ~= jj
            label = ['seq_' int2str(ii) '_' int2str(jj)];
            seq = [t(jj) - t(ii) + bigM * (1 - q(ii,jj)) - setuptime(ii,jj) ...
                            * (1 - x(ii)) >= 0] : label;
%             seqConstraints = [seqConstraints, seq];
%             seq = [t(jj) - t(ii) + bigM * (1 - q(ii,jj)) - setuptime(ii,jj)>=0];
%                             * (1 - x(ii) - x(jj)) >= 0] : label;
%             seq = [t(jj) - t(ii) + bigM * (1 - q(ii,jj)) - setuptime(ii,jj) ...
%                             * (1 - x(ii) - x(jj)) >= 0] : label;
            seq2 = [t(ii) - t(jj) + q(ii,jj) * bigM >= 0]: label;
            seqConstraints = [seqConstraints, seq, seq2];
%             seqConstraints = [seqConstraints, seq2];
%             
%             % Ensure that they are now on the same time
            if ii >= jj
                same = [q(ii,jj) + q(jj,ii) == 1]: label;
                seqConstraints = [seqConstraints, same];
            end

%             if direction(ii) == direction(jj)
%                 % We have to fix the q-value!
%                 label = ['fix q(' int2str(ii) ',' int2str(jj) ')'];
%                 if releasetime(ii) < releasetime(jj)
%                     Q = [q(ii,jj) >= 1 - x(ii) - x(jj)] : label;
%                 else
%                     Q = [q(ii,jj) == 0] : label;
%                 end
%                 seqConstraints = [seqConstraints, Q];
%             end
        end
    end
end
Constraints = [Constraints, seqConstraints];

%% Build the objective function
Objective = 0;

% Penalize the delays
for tt = 1:Ntrains
    Objective = Objective + C(tt) - deptime(tt);
end

% Penalize the cancellation of trains
w_cancel = 3600;
for tt = 1:Ntrains
    Objective = Objective + w_cancel * x(tt);
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
%     type = 'complete';
    type = 'base';
    include_blocks = 1;
    [line_new, blocks_new] = plotTT(new_timetable, blocksections, settings, type, include_blocks);
    
end

end
