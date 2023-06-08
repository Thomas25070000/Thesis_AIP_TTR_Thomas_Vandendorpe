function balanceConstraints = generateOffBalanceConstraintsRollingHorizon(trains, releasetime, x, settings)

delta = settings.disruption.offbalance;

% In which direction do they go?
direction = mod(trains,10);

dir0 = find(direction == 0);
dir12 = find(direction == 2);
dir13 = find(direction == 3);

release0 = releasetime(dir0);
release12 = releasetime(dir12);
release13 = releasetime(dir13);


% Frequency of trains?
freq = length(trains) / (2*(settings.disruption.duration + settings.disruption.nr_hours_after));
if min(release13)>=min(release12)
    if min(release12) <= min(release0)
        C = [];
        for tt = 1:length(dir12)-freq+1
            % Which run in direction 0 in this interval?
            time = release12(tt);
            t0 = dir0(find(release0 >= time & release0 < (time + 3600 - 1)));
            t1 = dir12(find(release12 >= time & release12 < (time + 3600 - 1)));
    
            LB = [sum(x((t0))) - sum(x((t1))) >= -delta];
            UB = [sum(x((t0))) - sum(x((t1))) <= delta];
    
            C = [C LB UB];
        end
    else
        C = [];
        for tt = 1:length(dir0)-freq+1
            % Which run in direction 0 in this interval?
            time = release0(tt);
            t0 = dir0(find(release0 >= time & release0 < (time + 3600 - 1)));
            t1 = dir12(find(release12 >= time & release12 < (time + 3600 - 1)));
    
            LB = [sum(x((t0))) - sum(x((t1))) >= -delta];
            UB = [sum(x((t0))) - sum(x((t1))) <= delta];
    
            C = [C LB UB];
        end
    end
end
if min(release13)<min(release12)
    if min(release13) <= min(release0)
        C = [];
        for tt = 1:length(dir13)-freq+1
            % Which run in direction 0 in this interval?
            time = release1(tt);
            t0 = dir0(find(release0 >= time & release0 < (time + 3600 - 1)));
            t1 = dir13(find(release13 >= time & release13 < (time + 3600 - 1)));
    
            LB = [sum(x((t0))) - sum(x((t1))) >= -delta];
            UB = [sum(x((t0))) - sum(x((t1))) <= delta];
    
            C = [C LB UB];
        end
    else
        C = [];
        for tt = 1:length(dir0)-freq+1
            % Which run in direction 0 in this interval?
            time = release0(tt);
            t0 = dir0(find(release0 >= time & release0 < (time + 3600 - 1)));
            t1 = dir13(find(release13 >= time & release13 < (time + 3600 - 1)));
    
            LB = [sum(x((t0))) - sum(x((t1))) >= -delta];
            UB = [sum(x((t0))) - sum(x((t1))) <= delta];
    
            C = [C LB UB];
        end
    end
end

balanceConstraints = C;



end