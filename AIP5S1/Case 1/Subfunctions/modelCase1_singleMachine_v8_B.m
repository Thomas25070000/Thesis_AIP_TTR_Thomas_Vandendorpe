function [timetable, solu] = modelCase1_singleMachine_v8_B(timetable, regular, blocksections, trains, minHW, settings)
%% This version also accounts for the increase in blocking times next to the closed section.

regTT = regular.TT;
regHW = regular.HW;

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
[minHWstack, minHWstack_reg, otherinfo] = createHeadwayStackDisruption(timetable, regular, blocksections, blocktype, bref, settings);

arr_ev = otherinfo.arr_ev;
release_ev = otherinfo.release_ev;
direction = otherinfo.direction;
deptime = otherinfo.deptime;
completion_ev = otherinfo.completion_ev;
outside_disruption_period = otherinfo.outside_disruption_period;
releasetime = otherinfo.releasetime;

processtime = otherinfo.processtime;
regprocesstime = otherinfo.regprocesstime;
approachtime = otherinfo.approachtime;
regapproachtime = otherinfo.regapproachtime;
cleartime = otherinfo.cleartime;

t_end = settings.disruption.duration * 3600;

w_cancel = 3000;


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

%% Constraints
bigM = settings.disruption.duration * 3600 * 20;     % Just make it very big!
Constraints = [];

% Start time >= release time on the first section of the train + indicator
% whether it is after the disruption ending.
for tt = 1:Ntrains
    if direction(tt)
        % First block is the left one
        fb = 1;
    else
        % First block is at the end.
        fb = Nmachines;
    end
    minStartTime = [t(tt,fb) - releasetime(tt) * (1-x(tt)) >= 0] : ['start_' int2str(trains(tt))];
    disruptionend1 = [t(tt,fb) - t_end * gamma(tt) >= 0]: ['end_of_disruption_' int2str(trains(tt))];
    disruptionend2 = [t(tt,fb) - bigM * gamma(tt) <= t_end]: ['end_of_disruption_' int2str(trains(tt))];
    Constraints = [Constraints, minStartTime, disruptionend1, disruptionend2];
    
    % Fix gamma?
    if releasetime(tt) + w_cancel < t_end
        % Never able to cross the disruption border, better to cancel it
        % than delay it that much.
        fixGamma = [gamma(tt) == 0]: ['fix_gamma_' int2str(trains(tt))];
        Constraints = [Constraints, fixGamma];
    elseif releasetime(tt) >= t_end
        fixGamma = [gamma(tt) == 1]: ['fix_gamma_' int2str(trains(tt))];
        Constraints = [Constraints, fixGamma];
    end
    
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
t_before = settings.TT.blocktimes.setup;

for ii = 1:Ntrains
    for jj = ii:Ntrains
        if jj>ii
            
%             if (direction(ii) == settings.disruption.direction) && ~settings.disruption.signals
%                 % We have to compare with the actual exit of ii!
%                 flag_exit_ii = 1;
%             else
%                 flag_exit_ii = 0;
%             end
%             if (direction(jj) == settings.disruption.direction) && ~settings.disruption.signals
%                 % We have to compare with the actual exit of jj!
%                 flag_exit_jj = 1;
%             else
%                 flag_exit_jj = 0;
%             end
%             
%             
%             if direction(ii) == 1
%                 first_machine_after_ii = find(bref==0) + 1;
%             else
%                 first_machine_after_ii = find(bref==0) - 1;
%             end
%             if direction(jj) == 1
%                 first_machine_after_jj = find(bref==0) + 1;
%             else
%                 first_machine_after_jj = find(bref==0) - 1;
%             end
%             
%             if first_machine_after_ii > Nmachines
%                 % In the other direction, it will be 0 anyway!
%                 first_machine_after_ii = 0;
%             end
%             if first_machine_after_jj > Nmachines
%                 first_machine_after_jj = 0;
%             end
            
%             t_setup_ii_jj = cleartime(ii) + approachtime(jj,closedM) + t_before;
%             t_setup_jj_ii = cleartime(jj) + approachtime(ii,closedM) + t_before;
            label = ['seq_' int2str(ii) '_' int2str(jj) '_disrupted'];
            
            % What is the reference event?
%             if (flag_exit_ii || (direction(ii) ~= direction(jj))) && first_machine_after_ii==0
%                 ref_event_ii = C(ii);
%             elseif (flag_exit_ii || (direction(ii) ~= direction(jj)))
%                 ref_event_ii = t(ii,first_machine_after_ii);
%             else
%                 % Look at the time ii exits the first block section of the
%                 % closed part.
%                 closedev = timetable(find(timetable.train_id == trains(ii) ...
%                                         & blocksections.closed(timetable.blocksection)),:);
%                 ref_event_ii = t(ii,closedM) + max(closedev.running);
%             end
            
            ref_event_ii = t(ii,closedM);
            ref_event_jj = t(jj,closedM);
            
%             closedev = regTT(find(regTT.train_id == trains(ii) ...
%                                         & blocksections.closed(regTT.blocksection)),:);
%             ref_event_ii_end = t(ii,closedM) + closedev.running(1);
%             closedev = regTT(find(regTT.train_id == trains(jj) ...
%                                         & blocksections.closed(regTT.blocksection)),:);
%             ref_event_jj_end = t(jj,closedM) + closedev.running(1);
            
            
            % What is the reference event?
%             if (flag_exit_jj || (direction(ii) ~= direction(jj))) && first_machine_after_jj==0
%                 ref_event_jj = C(jj);
%             elseif (flag_exit_jj || (direction(ii) ~= direction(jj)))
%                 ref_event_jj = t(jj,first_machine_after_jj);
%             else
%                 % Look at the time ii exits the first block section of the
%                 % closed part.
%                 closedev = timetable(find(timetable.train_id == trains(jj) ...
%                                         & blocksections.closed(timetable.blocksection)),:);
%                 ref_event_jj = t(jj,closedM) + max(closedev.running);
%             end
            t_setup_ii_jj = minHWstack(ii,jj,closedM);
            t_setup_jj_ii = minHWstack(jj,ii,closedM);


            % If trains fall outside of the disruption period, it does not
            % matter anymore if they go in different directions!
            flag_using_same_track = 0;
            if direction(ii) == direction(jj)
                flag_using_same_track = 1;
            elseif direction(ii) ~= direction(jj) && ~outside_disruption_period(ii) ...
                        && ~outside_disruption_period(jj)
                flag_using_same_track = 1;
            end
                
%             if flag_using_same_track
%                 seq1 = [t(jj,closedM) + bigM * (1 - q(ii,jj) + gamma(ii)) ...
%                         - ref_event_ii - t_setup_ii_jj * (1 - x(ii)) >= 0]: label;
%                 seq2 = [t(ii,closedM) + bigM * (q(ii,jj) + gamma(jj)) ...
%                         - ref_event_jj - t_setup_jj_ii * (1 - x(jj)) >= 0]: label;        
%                 seqConstraints = [seqConstraints, seq1, seq2];

            if direction(ii) ~= direction(jj)
                if direction(ii) == settings.disruption.direction
                    seq1 = [t(jj,closedM) + bigM * (1 - q(ii,jj) + gamma(ii)) ...
                            - ref_event_ii - t_setup_ii_jj * (1 - x(ii)) >= 0]: label;
                    seq2 = [t(ii,closedM) + bigM * (q(ii,jj) + gamma(ii)) ...
                            - ref_event_jj - t_setup_jj_ii * (1 - x(jj)) >= 0]: label;    
                else
                    seq1 = [t(jj,closedM) + bigM * (1 - q(ii,jj) + gamma(jj)) ...
                            - ref_event_ii - t_setup_ii_jj * (1 - x(ii)) >= 0]: label;
                    seq2 = [t(ii,closedM) + bigM * (q(ii,jj) + gamma(jj)) ...
                            - ref_event_jj - t_setup_jj_ii * (1 - x(jj)) >= 0]: label;
                end
            else
                seq1 = [t(jj,closedM) + bigM * (1 - q(ii,jj) + gamma(ii) + gamma(jj)) ...
                            - ref_event_ii - t_setup_ii_jj * (1 - x(ii)) >= 0]: label;
                seq2 = [t(ii,closedM) + bigM * (q(ii,jj) + gamma(jj) + gamma(ii)) ...
                            - ref_event_jj - t_setup_jj_ii * (1 - x(jj)) >= 0]: label;  
            end
            
            seqConstraints = [seqConstraints, seq1, seq2];
%             end

            % In case of ending after the end of the disruption
%             if direction(ii) == direction(jj)
%                 seq_end1 = [t(jj,closedM) + bigM * (1 - q(ii,jj)) ...
%                             - t(ii, closedM) + regHW(ii,jj) * gamma(ii) >= 0];
%                 seq_end2 = [t(ii,closedM) + bigM * q(ii,jj) ...
%                             - t(jj, closedM) + regHW(jj,ii) * gamma(jj) >= 0];     
%                 seqConstraints = [seqConstraints, seq_end1, seq_end2];  
%             end
            

       
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
                        
                        for mm = 1:Nmachines%-1
                            if bref(mm) > 0
                                try
%                                     delta_approach = regapproachtime(jj,mm) - approachtime(jj,mm);
%                                     t_setup = cleartime(ii) + t_before + approachtime(jj,mm);
                                    delta_hw = minHWstack_reg(ii,jj,mm) - minHWstack(ii,jj,mm);
                                    t_setup = minHWstack(ii,jj,mm);
                                catch
                                    pause;
%                                     delta_approach = regapproachtime(jj,mm) - approachtime(jj,1);
%                                     t_setup = cleartime(ii) + t_before + approachtime(jj,1);
                                end
                               
%                                 seqConstraint = [t(jj,mm) - t(ii,mm+1) - t_setup ...
%                                         *(1-x(ii)) - delta_approach * gamma(jj) >= 0];   % GAMMA
                                seqConstraint = [t(jj,mm) - t(ii,mm) - t_setup ...
                                    *(1-x(ii)) - delta_hw * gamma(jj) >= 0];   % GAMMA
                                corridorSeq = [corridorSeq seqConstraint];
                            end
                        end
                        % Also for the completing one!
                        
%                         delta_approach = regapproachtime(jj,end) - approachtime(jj,end);
%                         t_setup = cleartime(ii) + t_before + approachtime(jj,end);
%                         seqConstraint = [t(jj,end) - C(ii) - t_setup*(1-x(ii))...
%                                     - delta_approach * gamma(jj)>= 0];
%                         corridorSeq = [corridorSeq seqConstraint];
                    else
  
                        for mm = 1:Nmachines%-1
                            if bref(mm) > 0         % Don't do it for the disrupted section!
                                try
%                                     delta_approach = regapproachtime(ii,mm) - approachtime(ii,mm);
%                                     t_setup = cleartime(jj) + t_before + approachtime(ii,mm);
                                    delta_hw = minHWstack_reg(jj,ii,mm) - minHWstack(jj,ii,mm);
                                    t_setup = minHWstack(jj,ii,mm);
                                catch
                                    pause;
%                                     delta_approach = regapproachtime(ii,1) - approachtime(ii,1);
%                                     t_setup = cleartime(jj) + t_before + approachtime(ii,1);
                                end
                                
%                                 seqConstraint = [t(ii,mm) - t(jj,mm+1) - t_setup*(1-x(jj))...
%                                             - delta_approach * gamma(ii)>= 0];
                                seqConstraint = [t(ii,mm) - t(jj,mm) - t_setup*(1-x(jj))...
                                            - delta_hw * gamma(ii)>= 0];
                                corridorSeq = [corridorSeq seqConstraint];
                            end
                        end
                        % Also for the completing one!
%                         delta_approach = regapproachtime(jj,end) - approachtime(jj,end);
%                         t_setup = cleartime(jj) + t_before + approachtime(ii,end);
%                         
%                         seqConstraint = [t(ii,end) - C(jj) - t_setup*(1-x(jj)) ...
%                                         - delta_approach * gamma(ii)>= 0];
%                         corridorSeq = [corridorSeq seqConstraint];
                    end
                else
                    % Trains go in the other direction!
                    if releasetime(jj) > releasetime(ii)
%                         t_setup = cleartime(ii) + t_before + processtime(jj,end);
                        for mm = Nmachines:-1:1 %2
                            if bref(mm)>0
                                try
%                                     delta_approach = regapproachtime(jj,mm) - approachtime(jj,mm);
%                                     t_setup = cleartime(ii) + t_before + approachtime(jj,mm);
                                    delta_hw = minHWstack_reg(ii,jj,mm) - minHWstack(ii,jj,mm);
                                    t_setup = minHWstack(ii,jj,mm);
                                catch
                                    pause;
%                                     delta_approach = regapproachtime(jj,end) - approachtime(jj,end);
%                                     t_setup = cleartime(ii) + t_before + approachtime(jj,end);
                                end
                                
%                                 seqConstraint = [t(jj,mm) - t(ii,mm-1) - t_setup*(1-x(ii))...
%                                              - delta_approach * gamma(jj)>= 0];
                                seqConstraint = [t(jj,mm) - t(ii,mm) - t_setup*(1-x(ii))...
                                             - delta_hw * gamma(jj)>= 0];
                                corridorSeq = [corridorSeq seqConstraint];
                            end
                        end
                        % Also for the completing one!
%                         delta_approach = regapproachtime(jj,1) - approachtime(jj,1);
%                         t_setup = cleartime(ii) + t_before + approachtime(jj,1);
%                         seqConstraint = [t(jj,1) - C(ii) - t_setup*(1-x(ii)) ...
%                                     -  delta_approach * gamma(ii)>= 0];
%                         corridorSeq = [corridorSeq seqConstraint];
                    else
%                         t_setup = cleartime(jj) + t_before + processtime(ii,end);
                        for mm = Nmachines:-1:1 %2
                            if bref(mm)>0
                                try
%                                     delta_approach = regapproachtime(ii,mm) - approachtime(ii,mm);
%                                     t_setup = cleartime(jj) + t_before + approachtime(ii,mm);
                                    delta_hw = minHWstack_reg(jj,ii,mm) - minHWstack(jj,ii,mm);
                                    t_setup = minHWstack(jj,ii,mm);
                                catch
                                    pause;
%                                     delta_approach = regapproachtime(ii,end) - approachtime(ii,end);
%                                     t_setup = cleartime(jj) + t_before + approachtime(ii,end);
                                end
                                
%                                 seqConstraint = [t(ii,mm) - t(jj,mm-1) - t_setup*(1-x(jj))...
%                                             - delta_approach * gamma(ii) >= 0];
                                seqConstraint = [t(ii,mm) - t(jj,mm) - t_setup*(1-x(jj))...
                                            - delta_hw * gamma(ii) >= 0];
                                corridorSeq = [corridorSeq seqConstraint];
                            end
                        end
                        % Also for the completing one!
%                         delta_approach = regapproachtime(ii,1) - approachtime(ii,1);
%                         t_setup = cleartime(jj) + t_before + approachtime(ii,1);
%                         seqConstraint = [t(ii,1) - C(jj) - t_setup*(1-x(jj)) ...
%                                         - delta_approach * gamma(ii) >= 0];
%                         corridorSeq = [corridorSeq seqConstraint];
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
    vgamma = round(value(gamma));
    
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
    new_timetable = updateTimetable_v5(timetable, regular, vgamma, blocksections, settings, trains, cancelled);
    
    
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
