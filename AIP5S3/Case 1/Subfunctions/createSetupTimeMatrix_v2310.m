function [setuptimes, proctimes] = createSetupTimeMatrix(timetable, regTT, blocksections, settings)

% Construct a timetable which common for both!
fullTT = [];
trains = unique(timetable.train_id);
row = 0;
t_before = settings.TT.blocktimes.setup;
t_after = settings.TT.blocktimes.afterIC;

for tt = 1:length(trains)
    train = trains(tt);
    TT_train = timetable(find(timetable.train_id == train),:);
    reg_TT_train = regTT(find(regTT.train_id == train),:);
    
    if floor(train/1000) >= settings.disruption.duration
        afterend = 1;
    else
        afterend = 0;
    end
    
    % First row starts at 0.
    for ee = 1:size(TT_train,1)
        row = row + 1;
        fullTT(row).train = train;
        fullTT(row).bs = TT_train.blocksection(ee);
        if ee == 1
            fullTT(row).arrival = 0;
            [app, reg_app] = returnFirstApproachTime(blocksections, settings, TT_train.train_type(1), mod(train,10));
        else
            fullTT(row).arrival = fullTT(row-1).departure;
            reg_app = fullTT(row-1).running;
            app = fullTT(row-1).new_running;
        end
        fullTT(row).running = reg_TT_train.running(ee);
        fullTT(row).departure = fullTT(row).arrival + fullTT(row).running;
        fullTT(row).start = fullTT(row).arrival - reg_app - t_before;
        fullTT(row).finish = fullTT(row).departure + t_after;
        if ee == 1 
            fullTT(row).new_arrival = 0;
        else
            fullTT(row).new_arrival = fullTT(row-1).new_departure;
        end
        fullTT(row).new_running = TT_train.running(ee);
        fullTT(row).new_departure = fullTT(row).new_arrival + fullTT(row).new_running;
        if afterend
            fullTT(row).new_start = fullTT(row).new_arrival - reg_app - t_before;
        else
            fullTT(row).new_start = fullTT(row).new_arrival - app - t_before;
        end
        fullTT(row).new_finish = fullTT(row).new_departure + t_after;
    end
end

fullTT = struct2table(fullTT);

if ~settings.disruption.signals
    for tt = 1:length(trains)
        if mod(trains(tt),10) == settings.disruption.direction && ...
                mod(trains(tt),1000) < settings.disruption.duration
            % Find the events on the closed blocks.
            closedB = find(blocksections.closed);
            ev = find(fullTT.train == trains(tt) & ismember(fullTT.bs,closedB));
            start = fullTT.new_start(ev(1));
            finish = fullTT.new_finish(ev(end));
            fullTT.new_start(ev) = start;
            fullTT.new_finish(ev) = finish;
        end
    end
end

bref = [];
bb = 0;
while ~blocksections.closed(bb+1)
    bb = bb + 1;
    bref = [bref blocksections.id(bb)];
end
bref = [bref 0];
closedone = length(bref);
closedB = find(blocksections.closed);
bb = closedB(end);
while bb+1 <= size(blocksections,1)
    bb = bb + 1;
    bref = [bref blocksections.id(bb)];
end
Nmachines = length(bref);

Ntrains = length(trains);
setuptimes.disrupted = zeros(Ntrains,Ntrains,Nmachines);
setuptimes.onegamma = zeros(Ntrains,Ntrains,Nmachines);
setuptimes.bothgammas = zeros(Ntrains,Ntrains,Nmachines);
proctimes.disrupted = zeros(Ntrains,Nmachines);
proctimes.regular = zeros(Ntrains,Nmachines);

for ii = 1:Ntrains
    for jj = 1:Ntrains
        train_ii = trains(ii);
        train_jj = trains(jj);
        for mm = 1:Nmachines
            if mm ~= closedone
                ev_ii = find(fullTT.bs == bref(mm) & fullTT.train == train_ii);
                ev_jj = find(fullTT.bs == bref(mm) & fullTT.train == train_jj);
            else
                ev_ii = find(ismember(fullTT.bs,closedB) & fullTT.train == train_ii);
                ev_jj = find(ismember(fullTT.bs,closedB) & fullTT.train == train_jj);
            end
            if (mm ~= closedone) || (mod(train_ii,10) == mod(train_jj,10))
                disrupted = max(fullTT.new_finish(ev_ii) - fullTT.new_arrival(ev_ii(1)) ...
                            + fullTT.new_arrival(ev_jj(1)) - fullTT.new_start(ev_jj));
                % jj is outside disruption
                onegamma = max(fullTT.new_finish(ev_ii) - fullTT.new_arrival(ev_ii(1)) ...
                                + fullTT.arrival(ev_jj(1)) - fullTT.start(ev_jj));
                bothgammas = max(fullTT.finish(ev_ii) - fullTT.arrival(ev_ii(1)) ...
                                + fullTT.arrival(ev_jj(1)) - fullTT.start(ev_jj));
                            
            elseif (mm == closedone) && (mod(train_ii,10) ~= mod(train_jj,10))
                disrupted = fullTT.new_finish(ev_ii(end)) - fullTT.new_arrival(ev_ii(1)) ...
                            + fullTT.new_arrival(ev_jj(1)) - fullTT.new_start(ev_jj(1));
                % jj is outside disruption
                if mod(train_jj,10) == settings.disruption.direction
                    % Train jj runs on the closed track, so no hinder
                    % anymore!
                    onegamma = 0;
                else
                    onegamma = fullTT.new_finish(ev_ii(end)) - fullTT.new_arrival(ev_ii(1)) ...
                                    + fullTT.arrival(ev_jj(1)) - fullTT.start(ev_jj(1));

                end
                
                bothgammas = 0;
            else
                error('Problem');
            end
            
            % Correct for the running time of ii!
            disrupted = disrupted - fullTT.new_running(ev_ii(1));
%             if mod(train_jj,10) ~= settings.disruption.direction
                onegamma = onegamma - fullTT.new_running(ev_ii(1));
%             end
            if onegamma<0
                    disp('');
            end
               
            bothgammas = bothgammas - fullTT.running(ev_ii(1));
            
            setuptimes.disrupted(ii,jj,mm) = disrupted;
            setuptimes.onegamma(ii,jj,mm) = onegamma;
            setuptimes.bothgammas(ii,jj,mm) = bothgammas;
            
            proctimes.disrupted(ii,mm) = fullTT.new_running(ev_ii(1));
            proctimes.regular(ii,mm) = fullTT.running(ev_ii(1));
        end
    end
end
                            


end